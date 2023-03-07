
#ifndef MEDYAN_SysParamsValidate_hpp
#define MEDYAN_SysParamsValidate_hpp

#include <unordered_set>
#include <sstream>

#include "SysParams.h"

namespace medyan {

inline bool checkMembraneSettings(const SimulConfig& conf) {
    using namespace std;

    bool passed = true;
    std::ostringstream errMsg;

    auto& membraneSettings = conf.membraneSettings;
    auto& setupVec = membraneSettings.setupVec;
    auto& initVec  = membraneSettings.initVec;

    // Check membrane setup parameters.
    for(auto& setup : setupVec) {
        switch(setup.vertexSystem) {
            case MembraneMeshVertexSystem::material:
                if(setup.hasLipidReservoir) {
                    passed = false;
                    errMsg << "Currently, material coordinate system does not support lipid reservoir.\n";
                }
                break;
            case MembraneMeshVertexSystem::normal:
                passed = false;
                errMsg << "Currently, normal coordinate system is not supported.\n";
                break;
            case MembraneMeshVertexSystem::general:
                // Okay.
                break;
        }
    }

    // Check membrane initialization parameters.
    for(auto& init : initVec) {
        // Check if the membrane profile has been specified.
        auto it = find_if(
            setupVec.begin(),
            setupVec.end(),
            [&](const MembraneSetup& setup) {
                return setup.name == init.name;
            }
        );
        if(it == setupVec.end()) {
            passed = false;
            errMsg << "Membrane profile " << init.name << " has not been specified.\n";
        }

        // Check if all the species names are present in the chemistry data.
        for(auto& species : init.speciesInitVec) {
            auto optIndex = conf.chemistryData.findIndexSpeciesMembraneDiffusing(species.name);
            if(!optIndex.has_value()) {
                passed = false;
                errMsg << "In initializing " << init.name << ", species " << species.name << " is not present in the chemistry data.\n";
            }
        }
    }

    // Print error message.
    if(!passed) {
        log::error("Errors occurred in membrane settings.");
        log::error(errMsg.str());
    }
    return passed;
}

inline bool checkChemParameters(const SimulConfig& conf) {
    using namespace std;

    auto& chemParams = conf.chemParams;
    auto& chem = conf.chemistryData;
    auto& meshChem = conf.membraneMeshChemistryInfo;

    bool passed = true;
    std::ostringstream errMsg;

    // Check filament species.
    if(chemParams.numFilaments < 1) {
        errMsg << "Must specify at least one type of filament.\n";
        passed = false;
    }
    
    for(auto filType = 0; filType < chemParams.numFilaments; filType++) {
    
        if(chem.speciesFilament[filType].size() == 0) {
            
            errMsg << "At least one filament species is required for each filament type.\n";
            passed = false;
        }

        if(chem.speciesPlusEnd[filType].size() == 0) {
            
            errMsg << "At least one plus end species is required for each filament type.\n";
            passed = false;
        }

        if(chem.speciesMinusEnd[filType].size() == 0) {
            
            errMsg << "At least one minus end species is required for each filament type.\n";
            passed = false;
        }

        
        if(chem.speciesBound[filType].size() == 0) {
            
            errMsg << "At least one Bound species is required for each filament type.\n";
            passed = false;
        }

        
        //check if binding sites are valid
        if(chem.B_BINDING_INDEX[filType] == "" && chem.speciesBrancher[filType].size() != 0) {
            errMsg << "A brancher binding site must be set for every filament type.\n";
            passed = false;
        }
        
        if(chem.L_BINDING_INDEX[filType] == "" && chem.speciesLinker[filType].size() != 0) {
            errMsg << "A linker binding site must be set for every filament type.\n";
            passed = false;
        }
        
        if(chem.M_BINDING_INDEX[filType] == "" && chem.speciesMotor[filType].size() != 0) {
            errMsg << "A motor binding site must be set for every filament type.\n";
            passed = false;
        }
    }

    // Check species name uniqueness.
    {
        unordered_set<string> allSpeciesNames;

        const auto checkAndAdd = [&](const string& name, string_view cat) {
            if(allSpeciesNames.find(name) != allSpeciesNames.end()) {
                errMsg << "Species " << name << " has been defined, when trying to define a species for " << cat << ".\n";
                passed = false;
            } else {
                allSpeciesNames.insert(name);
            }
        };

        for(auto& s : chem.speciesBulk)              checkAndAdd(s.name, "bulk");
        for(auto& s : chem.speciesDiffusing)         checkAndAdd(s.name, "diffusing");
        for(auto& s : chem.speciesMembraneDiffusing) checkAndAdd(s.name, "membrane diffusing");
        for(auto& snamevec : chem.speciesFilament) for(auto& sname : snamevec) checkAndAdd(sname, "filament");
        for(auto& snamevec : chem.speciesPlusEnd)  for(auto& sname : snamevec) checkAndAdd(sname, "plus end");
        for(auto& snamevec : chem.speciesMinusEnd) for(auto& sname : snamevec) checkAndAdd(sname, "minus end");
        for(auto& snamevec : chem.speciesBound)    for(auto& sname : snamevec) checkAndAdd(sname, "bound");
        for(auto& snamevec : chem.speciesLinker)   for(auto& sname : snamevec) checkAndAdd(sname, "linker");
        for(auto& snamevec : chem.speciesMotor)    for(auto& sname : snamevec) checkAndAdd(sname, "motor");
        for(auto& snamevec : chem.speciesBrancher) for(auto& sname : snamevec) checkAndAdd(sname, "brancher");
    }

    // Check additional motor params.
    const auto totalNumMotors = chem.numMotorSpecies();
    
    if(totalNumMotors != chemParams.motorNumHeadsMin.size()) {
        
        errMsg << format("Number of minimum motor heads ({}) does not match the number of motor species ({}).\n", chemParams.motorNumHeadsMin.size(), totalNumMotors);
        passed = false;
    }
    if(totalNumMotors != chemParams.motorNumHeadsMax.size()) {
        
        errMsg << format("Number of maximum motor heads ({}) does not match the number of motor species ({}).\n", chemParams.motorNumHeadsMax.size(), totalNumMotors);
        passed = false;
    }
    if(totalNumMotors != chemParams.motorStepSize.size()) {
        
        errMsg << format("Number of motor step sizes ({}) does not match the number of motor species ({}).\n", chemParams.motorStepSize.size(), totalNumMotors);
        passed = false;
    }

    // Check usage of various types of species.
    //----------------------------------
    const auto checkMembraneDiffusingSpeciesName = [&](const std::string& name) {
        auto opIndexMembraneDiffusing = chem.findIndexSpeciesMembraneDiffusing(name);
        if(!opIndexMembraneDiffusing.has_value()) {
            errMsg << name << " used in a surface reaction is not registered as a membrane diffusing species.\n";
            passed = false;
        }
    };
    const auto checkGeneralSpeciesName = [&](const std::string& name) {
        auto opIndexGeneral = chem.findIndexSpeciesGeneral(name);
        if(!opIndexGeneral.has_value()) {
            errMsg << name << " used in a reaction is not registered as a general species.\n";
            passed = false;
        }
    };
    const auto checkBulkSpeciesName = [&](std::string_view name) {
        auto opIndexBulk = chem.findIndexSpeciesBulk(name);
        if(!opIndexBulk.has_value()) {
            errMsg << name << " used in a reaction is not registered as a bulk species.\n";
            passed = false;
        }
    };
    const auto checkDiffusingSpeciesName = [&](std::string_view name) {
        auto opIndexDiffusing = chem.findIndexSpeciesDiffusing(name);
        if(!opIndexDiffusing.has_value()) {
            errMsg << name << " used in a reaction is not registered as a diffusing species.\n";
            passed = false;
        }
    };
    const auto checkDiffusingOrBulkSpeciesName = [&](const std::string& name) {
        auto opIndexDiffusing = chem.findIndexSpeciesDiffusing(name);
        auto opIndexBulk = chem.findIndexSpeciesBulk(name);
        if(!opIndexDiffusing.has_value() && !opIndexBulk.has_value()) {
            errMsg << name << " used in a reaction is not registered as a diffusing or bulk species.\n";
            passed = false;
        }
    };

    // Check surface reactions.
    for(auto& r : chem.reactionsSurface) {
        if(r.reactants.empty() && r.products.empty()) {
            errMsg << "A surface reaction has no reactant or products.\n";
            passed = false;
        }

        for(auto& name : r.reactants) {
            checkMembraneDiffusingSpeciesName(name);
        }
        for(auto& name : r.products) {
            checkMembraneDiffusingSpeciesName(name);
        }
    }

    // Check adsorption/desorption reactions.
    for(auto& r : chem.reactionsAdsorptionDesorption) {
        checkDiffusingOrBulkSpeciesName(r.speciesName3D);
        checkMembraneDiffusingSpeciesName(r.speciesName2D);
    }

    // Check postprocessed membrane mesh chemistry info.
    {
        const auto ns = meshChem.diffusingSpeciesNames.size();

        // Check diffusion
        for(const auto& di : meshChem.diffusion) {
            if(di.speciesIndex >= ns) {
                errMsg << format("Diffusion species index {} is out of range.", di.speciesIndex);
                passed = false;
            }
            if(di.diffusionCoeff < 0) {
                errMsg << format("Diffusion coefficient {} is unacceptable.", di.diffusionCoeff);
                passed = false;
            }
        }

        // Check internal reaction
        for(const auto& iri : meshChem.internalReactions) {
            for(auto i : iri.reactantSpeciesIndices) {
                if(i >= ns) {
                    errMsg << "Internal reaction reactant species index " << i << " is out of range.";
                    passed = false;
                }
            }
            for(auto i : iri.productSpeciesIndices) {
                if(i >= ns) {
                    errMsg << "Internal reaction product species index " << i << " is out of range.";
                    passed = false;
                }
            }
            if(iri.rateConstant < 0) {
                errMsg << "Internal reaction rate constant " << iri.rateConstant << " is unacceptable.";
                passed = false;
            }
        }
    }

    // Check all emission-absorption reactions.
    const auto checkReactionEmissionAbsorption = [&](const ReactionEmissionAbsorptionSetup& r) {
        checkGeneralSpeciesName(r.speciesName1);
        checkDiffusingOrBulkSpeciesName(r.speciesName2);
        if(r.speciesType2 == "diffusing") {
            checkDiffusingSpeciesName(r.speciesName2);
        } else if(r.speciesType2 == "bulk") {
            checkBulkSpeciesName(r.speciesName2);
        } else {
            errMsg << "Species type " << r.speciesType2 << " is not supported.";
            passed = false;
        }
    };
    for(auto& mtocInit : conf.mtocSettings.initVec) {
        for(auto& r : mtocInit.vecRxnEmiAbs) {
            checkReactionEmissionAbsorption(r.setup);
        }
    }
    for(auto& afmInit : conf.afmSettings.initVec) {
        for(auto& r : afmInit.vecRxnEmiAbs) {
            checkReactionEmissionAbsorption(r.setup);
        }
    }

    // Print error message.
    if(!passed) {
        log::error("Error(s) occurred in chemistry settings in the input file.");
        std::cout << errMsg.str() << std::flush;
    }
    return passed;
}

inline bool checkMechParameters(const medyan::SimulConfig& config) {
    
    //check ff and associated parameters for consistency

    bool passed = true;
    std::ostringstream errMsg;

    const auto& mech = config.mechParams;
    const auto& chem = config.chemParams;
    const auto& chemData = config.chemistryData;
    const auto& geo  = config.geoParams;
    const auto& ff = mech.mechanicsFFType;

    // Auxiliary function to check usage of membrane diffusing species.
    const auto checkMembraneDiffusingSpeciesName = [&](std::string_view name) {
        auto it = std::find_if(
            chemData.speciesMembraneDiffusing.begin(), chemData.speciesMembraneDiffusing.end(),
            [name](const ChemistryData::SpeciesMembraneDiffusingInfo& sinfo) { return sinfo.name == name; }
        );
        if(it == chemData.speciesMembraneDiffusing.end()) {
            errMsg << name << " is not registered as a membrane diffusing species.\n";
            passed = false;
        }
    };


    //FILAMENT
    if(ff.FStretchingType != "") {
        if(mech.FStretchingK.size() != chem.numFilaments) {
            errMsg << "Must set a filament stretching constant for all filaments. Exiting.\n";
            passed = false;
        }
    }
    if(ff.FBendingType != "") {
        if(mech.FBendingK.size() != chem.numFilaments) {
            errMsg << "Must set a filament bending constant for all filaments.\n";
            passed = false;
        }
        if(mech.FBendingTheta.size() != chem.numFilaments) {
            errMsg << "Must set a filament equilibrium bending angle for all filaments.\n";
            passed = false;
        }
    }
    if(ff.FTwistingType != "" &&
       mech.FTwistingK.size() != chem.numFilaments) {
        errMsg << "Must set a filament twisting constant for all filaments.\n";
        passed = false;
    }
    
    //LINKER
    short totalNumLinkers = config.chemistryData.numLinkerSpecies();
    
    if(ff.LStretchingType != "" &&
       mech.LStretchingK.size() != totalNumLinkers && totalNumLinkers > 0) {
        errMsg << format("Number of linker stretching constants ({}) does not match the number of linker species in system ({}).\n", mech.LStretchingK.size(), totalNumLinkers);
        passed = false;
    }
    if(ff.LBendingType != "" &&
       mech.LBendingK.size() != totalNumLinkers && totalNumLinkers > 0) {
        errMsg << format("Number of linker bending constants ({}) does not match the number of linker species in system ({}).\n", mech.LBendingK.size(), totalNumLinkers);
        passed = false;
    }
    if(ff.LBendingType != "" &&
       mech.LBendingTheta.size() != totalNumLinkers && totalNumLinkers > 0) {
        errMsg << format("Number of linker bending angles ({}) does not match the number of linker species in system ({}).\n", mech.LBendingTheta.size(), totalNumLinkers);
        passed = false;
    }
    if(ff.LTwistingType != "" &&
       mech.LTwistingK.size() != totalNumLinkers && totalNumLinkers > 0) {
        errMsg << format("Number of linker twisting constants ({}) does not match the number of linker species in system ({}).\n", mech.LTwistingK.size(), totalNumLinkers);
        passed = false;
    }
    if(ff.LTwistingType != "" &&
       mech.LTwistingPhi.size() != totalNumLinkers && totalNumLinkers > 0) {
        errMsg << format("Number of linker twisting angles ({}) does not match the number of linker species in system ({}).\n", mech.LTwistingPhi.size(), totalNumLinkers);
        passed = false;
    }
    
    //MOTOR
    short totalNumMotors = config.chemistryData.numMotorSpecies();
    
    if(ff.MStretchingType != "" &&
       mech.MStretchingK.size() != totalNumMotors && totalNumMotors > 0) {
        errMsg << format("Number of motor stretching constants ({}) does not match the number of motor species in system ({}).\n", mech.MStretchingK.size(), totalNumMotors);
        passed = false;
    }
    if(ff.MBendingType != "" &&
       mech.MBendingK.size() != totalNumMotors && totalNumMotors > 0) {
        errMsg << format("Number of motor bending constants ({}) does not match the number of motor species in system ({}).\n", mech.MBendingK.size(), totalNumMotors);
        passed = false;
    }
    if(ff.MBendingType != "" &&
       mech.MBendingTheta.size() != totalNumMotors && totalNumMotors > 0) {
        errMsg << format("Number of motor bending angles ({}) does not match the number of motor species in system ({}).\n", mech.MBendingTheta.size(), totalNumMotors);
        passed = false;
    }
    if(ff.MTwistingType != "" &&
       mech.MTwistingK.size() != totalNumMotors && totalNumMotors > 0) {
        errMsg << format("Number of motor twisting constants ({}) does not match the number of motor species in system ({}).\n", mech.MTwistingK.size(), totalNumMotors);
        passed = false;
    }
    if(ff.MTwistingType != "" &&
       mech.MTwistingPhi.size() != totalNumMotors && totalNumMotors > 0) {
        errMsg << format("Number of motor twisting angles ({}) does not match the number of motor species in system ({}).\n", mech.MTwistingPhi.size(), totalNumMotors);
        passed = false;
    }

    
    //BRANCHINGPOINT
    short totalNumBranchers = 0;
    for(auto& brancherEachFilament : chemData.speciesBrancher) {
        totalNumBranchers += brancherEachFilament.size();
    }    

    if(ff.BrStretchingType != "" &&
       mech.BrStretchingK.size() != totalNumBranchers && totalNumBranchers > 0) {
        errMsg << format("Number of branching point stretching constants ({}) does not match the number of brancher species in system ({}).\n", mech.BrStretchingK.size(), totalNumBranchers);
        passed = false;
    }
    if(ff.BrStretchingType != "" &&
       mech.BrStretchingL.size() != totalNumBranchers && totalNumBranchers > 0) {
        errMsg << format("Number of branching point stretching length ({}) does not match the number of brancher species in system ({}).\n", mech.BrStretchingL.size(), totalNumBranchers);
        passed = false;
    }
    if(ff.BrBendingType != "" &&
       mech.BrBendingK.size() != totalNumBranchers && totalNumBranchers > 0) {
        errMsg << format("Number of branching point bending constants ({}) does not match the number of brancher species in system ({}).\n", mech.BrBendingK.size(), totalNumBranchers);
        passed = false;
    }
    if(ff.BrBendingType != "" &&
       mech.BrBendingTheta.size() != totalNumBranchers && totalNumBranchers > 0) {
        errMsg << format("Number of branching point bending angles ({}) does not match the number of brancher species in system ({}).\n", mech.BrBendingTheta.size(), totalNumBranchers);
        passed = false;
    }
    if(ff.BrDihedralType != "" &&
       mech.BrDihedralK.size() != totalNumBranchers && totalNumBranchers > 0) {
        errMsg << format("Number of branching point dihedral constants ({}) does not match the number of brancher species in system ({}).\n", mech.BrDihedralK.size(), totalNumBranchers);
        passed = false;
    }
    if(ff.BrPositionType != "" &&
       mech.BrPositionK.size() != totalNumBranchers && totalNumBranchers > 0) {
        errMsg << format("Number of branching point position constants ({}) does not match the number of brancher species in system ({}).\n", mech.BrPositionK.size(), totalNumBranchers);
        passed = false;
    }
    
    //VOLUME
    if(ff.cylinderVolumeExclusionFFType != medyan::CylinderVolumeExclusionFFType::none &&
       mech.VolumeK.size() != chem.numFilaments) {
        errMsg << "Must set a cylinder volume force constant for every filament type.\n";
        passed = false;
    }
    if(ff.cylinderVolumeExclusionFFType != medyan::CylinderVolumeExclusionFFType::none && areEqual(mech.VolumeCutoff, 0.0)) {
        errMsg << "Must set a cylinder volume cutoff for mechanical equilibration.\n";
        passed = false;
    }
    if(ff.cylinderVolumeExclusionFFType == medyan::CylinderVolumeExclusionFFType::monomer && mech.volumeExclusionMonomerInterval.size() != chem.numFilaments) {
        errMsg << "Must set a cylinder volume exclusion monomer interval for every filament type.\n";
        passed = false;
    }
    if(ff.triangleBeadVolumeFFType != "" && mech.triangleBeadVolume.cutoff == 0.0) {
        errMsg << "The membrane-bead volume cutoff for load force is invalid.\n";
        passed = false;
    }
    if(ff.triangleBeadVolumeFFType != "" && mech.triangleBeadVolume.cutoffMech == 0.0) {
        errMsg << "The membrane-bead volume cutoff for mechanical equilibration is invalid.\n";
        passed = false;
    }
    
    //Boundary
    if(ff.BoundaryFFType != "" && areEqual(config.boundParams.BoundaryK, 0.0)) {
        errMsg << "Must set a boundary force constant.\n";
        passed = false;
    }
    if(ff.BoundaryFFType != "" && areEqual(config.boundParams.BScreenLength, 0.0)) {
        errMsg << "Must set a boundary screen length.\n";
        passed = false;
    }
    if(ff.BoundaryFFType != "" && areEqual(config.boundParams.BoundaryCutoff, 0.0)) {
        errMsg << "Must set a boundary cutoff for mechanical equilibration.\n";
        passed = false;
    }
    
    //Bubbles
    if(ff.BubbleFFType != "" &&
      (mech.BubbleK.size() != mech.numBubbleTypes ||
       mech.BubbleRadius.size() != mech.numBubbleTypes ||
       mech.BubbleScreenLength.size() != mech.numBubbleTypes)) {
        errMsg << "Must set all bubble mechanical constants for every bubble type.\n";
        passed = false;
    }
    if(ff.BubbleFFType != "" && areEqual(mech.BubbleCutoff, 0.0)) {
        errMsg << "Must set a bubble cutoff for mechanical equilibration.\n";
        passed = false;
    }

    // Check protein curvature mismatch parameters.
    for(auto& setup : mech.proteinCurvatureMismatchSetups) {
        checkMembraneDiffusingSpeciesName(setup.speciesName);
    }

    ///Cylinder and monomer lengths specified
    if(geo.cylinderSize.size() != chem.numFilaments) {
        
        errMsg << "Must specify a cylinder size for every type of filament.\n";
        passed = false;
    }
    if(geo.monomerSize.size() != chem.numFilaments) {
        
        errMsg << "Must specify a monomer size for every type of filament.\n";
        passed = false;
    }


    // Print error message.
    if(!passed) {
        log::error("Error(s) occurred in mechanics settings in the input file.");
        std::cout << errMsg.str() << std::flush;
    }
    return passed;
}

inline bool checkGeoParameters(const GeoParams& geo) {

    bool passed = true;
    std::ostringstream errMsg;

    // Check grid and compartmentSize.
    if((geo.NX > 0 && geo.NY > 0 && geo.NZ > 0 &&
        geo.compartmentSizeX > 0 &&
        geo.compartmentSizeY > 0 &&
        geo.compartmentSizeZ > 0)){
    }
    else {
        errMsg << "Grid parameters are invalid.\n";
        passed = false;
    }

    // Check mesh adapter parameters.
    {
        auto& s = geo.meshAdapterSettings;
        if(std::cos(s.curvatureResolution) <= s.edgeFlipMinDotNormal) {
            errMsg << "MeshAdapterSettings: cos(curvatureResolution) must be greater than edgeFlipMinDotNormal.\n";
            passed = false;
        }
        if(s.maxSize <= 0) {
            errMsg << "MeshAdapterSettings: maxSize must be positive.\n";
            passed = false;
        }
    }

    // Print error message.
    if(!passed) {
        log::error("Error(s) occurred in geometry settings in the input file.");
        std::cout << errMsg.str() << std::flush;
    }
    return passed;
}

inline bool checkDyRateParameters(const medyan::SimulConfig& sc) {
    auto& chemParams = sc.chemParams;
    auto& dy = sc.dyRateParams.dynamicRateType;

    bool passed = true;
    std::ostringstream errMsg;

    //check types match number of species
    if(dy.dFPolymerizationType.size() != chemParams.numFilaments &&
       !dy.dFPolymerizationType.empty()) {
        errMsg << format("Number of filament dynamic rate polymerization forms ({}) must match the number of filaments ({}).\n", dy.dFPolymerizationType.size(), chemParams.numFilaments);
        passed = false;
    }

    if(dy.dLUnbindingType.size() != sc.chemistryData.numLinkerSpecies() &&
       !dy.dLUnbindingType.empty() && sc.chemistryData.numLinkerSpecies() > 0){
        errMsg << format("Number of linker dynamic rate unbinding forms ({}) must match the number of species ({}).\n", dy.dLUnbindingType.size(), sc.chemistryData.numLinkerSpecies());
        passed = false;
    }
    
    if(dy.dMUnbindingType.size() != sc.chemistryData.numMotorSpecies() &&
       !dy.dMUnbindingType.empty() && sc.chemistryData.numMotorSpecies() > 0) {
        errMsg << format("Number of motor dynamic rate unbinding forms ({}) must match the number of species ({}).\n", dy.dMUnbindingType.size(), sc.chemistryData.numMotorSpecies());
        passed = false;
    }
    if(dy.dMWalkingType.size() != sc.chemistryData.numMotorSpecies() &&
       !dy.dMWalkingType.empty() && sc.chemistryData.numMotorSpecies() > 0) {
        errMsg << format("Number of motor dynamic rate walking forms ({}) must match the number of species ({}).\n", dy.dMWalkingType.size(), sc.chemistryData.numMotorSpecies());
        passed = false;
    }
    
    //now check parameters
    if(dy.dFPolymerizationType.size() != sc.dyRateParams.dFilPolymerizationCharLength.size()) {
        errMsg << "Must set a dynamic rate polymerization length for all filaments.\n";
        passed = false;
    }
    
    auto numCharLengths = 0;
    auto numAmps = 0;
    
    for(auto &changer : dy.dLUnbindingType) {
        
        if(changer == "CATCHSLIP") {
            numCharLengths += 2;
            numAmps += 2;
        }
        else if(changer == "SLIP") {
            numCharLengths += 1;
        }
        
    }
    if(numCharLengths != sc.dyRateParams.dLinkerUnbindingCharLength.size()) {
        errMsg << "Number of characteristic lengths specified for chosen linker unbinding dynamic rate forms is " << sc.dyRateParams.dLinkerUnbindingCharLength.size()
            << ", but " << numCharLengths << " is required.\n";
        passed = false;
    }
    
    if(numAmps != sc.dyRateParams.dLinkerUnbindingAmplitude.size()) {
        
        
        cout << "Number of amplitudes specified for chosen "
             << "linker unbinding dynamic rate forms is not accurate. Exiting."
        << endl;
        return false;
    }
    
    auto numCharForces = 0;
    
    for(auto &changer : dy.dMUnbindingType) {
        
        if(changer == "LOWDUTYCATCHSLIP") {
            numCharForces += 2;
        }
        else if(changer == "LOWDUTYCATCH") {
            numCharForces += 1;
        }
        else if(changer == "HIGHDUTYCATCH") {
            numCharForces += 1;
        }
        
    }
    if(numCharForces != sc.dyRateParams.dMotorUnbindingCharForce.size()) {
        errMsg << "Number of characteristic forces specified for chosen motor unbinding dynamic rate forms is " << sc.dyRateParams.dMotorUnbindingCharForce.size()
            << ", but " << numCharForces << " is required.\n";
        passed = false;
    }
    
    if(dy.dMWalkingType.size() != sc.dyRateParams.dMotorWalkingCharForce.size()) {
        errMsg << "Number of characteristic forces specified for chosen motor walking dynamic rate forms is " << sc.dyRateParams.dMotorWalkingCharForce.size()
            << ", but " << dy.dMWalkingType.size() << " is required.\n";
        passed = false;
    }


    // Print error message.
    if(!passed) {
        log::error("Error(s) occurred in mechanics settings in the input file.");
        std::cout << errMsg.str() << std::flush;
    }
    return passed;
}


// Verify the correctness of the input parameters.
inline bool checkSimulConfig(const SimulConfig& conf) {
    bool passed = true;

    passed &= checkGeoParameters(conf.geoParams);
    passed &= checkChemParameters(conf);
    passed &= checkMechParameters(conf);
    passed &= checkDyRateParameters(conf);
    passed &= checkMembraneSettings(conf);

    return passed;
}

} // namespace medyan

#endif

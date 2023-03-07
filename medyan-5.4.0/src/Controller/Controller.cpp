
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include <algorithm>
#include <chrono>
#include <unordered_set>

#ifdef CUDAACCL
    #include "nvToolsExt.h"
#endif

#include "BubbleInitializer.h"
#include "FilamentInitializer.h"
#include "MathFunctions.h"
#include "Output.h"
#include "Parser.h"
#include "Rand.h"
#include "SysParams.h"
#include "SysParamsValidate.hpp"
#include "VisualSystemRawData.hpp"
#include "Chemistry/ChemManager.h"
#include "Controller/Controller.h"
#include "Mechanics/Minimizer/CGMethod.hpp"
#include "Structure/Boundary.h"
#include "Structure/BoundaryElementImpl.h"
#include "Structure/BranchingPoint.h"
#include "Structure/Bubble.h"
#include "Structure/CompartmentGrid.h"
#include "Structure/Cylinder.h"
#include "Structure/Filament.h"
#include "Structure/Linker.h"
#include "Structure/MotorGhost.h"
#include "Structure/Movables.hpp"
#include "Structure/SubSystem.h"
#include "Structure/SubSystemFunc.hpp"
#include "Structure/SurfaceMesh/FixedVertexAttachmentInit.hpp"
#include "Structure/SurfaceMesh/FuncMembraneMech.hpp"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Structure/SurfaceMesh/MembraneHierarchy.hpp"
#include "Structure/SurfaceMesh/MembraneRegion.hpp"
#include "Structure/SurfaceMesh/SurfaceMeshGeneratorPreset.hpp"
#include "Util/Io/Log.hpp"
#include "Util/Profiler.hpp"

namespace medyan {
using namespace mathfunc;


inline void remeshMembrane(SubSystem& sys, const adaptive_mesh::MembraneMeshAdapter& adapter, Membrane& membrane) {
    // Requires _meshAdapter to be already initialized
    adapter.adapt(sys, membrane.getMesh());

    // Update necessary geometry for the system
    medyan::updateGeometryValueForSystem(sys, membrane.getMesh());

    for(auto& t : membrane.getMesh().getTriangles()) {
        medyan::movable::updatePosition(sys, t.attr.triangle(sys));
    }
}


SimulConfig Controller::initialize(const CommandLineConfig& cmdConfig) {
    using namespace std;

    // Set floating point environment.
    trapInvalidFP(cmdConfig.trapInvalidFP);

    SysParams::INITIALIZEDSTATUS = false;

    //----------------------------------
    // Parse input, get parameters.
    //----------------------------------
    auto simulConfig = getSimulConfigFromInput(cmdConfig.inputFile, cmdConfig.inputDirectory);
    if(!checkSimulConfig(simulConfig)) {
        log::error("Simulation configuration is not valid.");
        throw runtime_error("Simulation configuration is not valid.");
    }
    SysParams::GParams = simulConfig.geoParams;
    SysParams::BParams = simulConfig.boundParams;
    SysParams::MParams = simulConfig.mechParams;
    SysParams::CParams = simulConfig.chemParams;
    SysParams::DRParams = simulConfig.dyRateParams;
    SysParams::filamentSetup = simulConfig.filamentSetup;

    auto& chemParams = simulConfig.chemParams;
    auto& ChemData = simulConfig.chemistryData;

    cout << endl;


    //----------------------------------
    // Initialize controllers.
    //----------------------------------
    //Initialize geometry controller
    cout << "---" << endl;
    log::info("Initializing geometry...");
    _gController.initializeGrid(simulConfig.geoParams);
    log::info("Done.");

    //Initialize boundary
    cout << "---" << endl;
    log::info("Initializing boundary...");
    _gController.initializeBoundary(simulConfig.boundParams.boundaryType);
    log::info("Done.");


    //Initialize Mechanical controller
    cout << "---" << endl;
    log::info("Initializing mechanics...");
    _mController.initialize(simulConfig);
    log::info("Done.");

    // Initialize special protocols.
    log::info("Initializing special protocols...");
    initializeSpecialProtocols(simulConfig);
    log::info("Done.");


    //Activate necessary compartments for diffusion
    _gController.setActiveCompartments();

    if(_subSystem.getBoundary()->getShape() == BoundaryShape::Cylinder){
        for(auto& C : _subSystem.getCompartmentGrid()->getCompartments()){
            C->computeSlicedVolumeArea(_subSystem, Compartment::SliceMethod::cylinderBoundary);
        }
    }
    else{
        for(auto& C : _subSystem.getCompartmentGrid()->getCompartments()){
            C->computeNonSlicedVolumeArea();
        }
    }
    //Calculate surface area and volume for reaction rate scaling


    //Initialize chemical controller
    cout << "---" << endl;
    log::info("Initializing chemistry...");
    //read algorithm
    auto& CAlgorithm = chemParams.chemistryAlgorithm;
    auto& CSetup = chemParams.chemistrySetup;
    //run time for sim
    _runTime = CAlgorithm.runTime;

    //freq of snapshots, minimizations, neighborlist updates
    _minimizationTime = CAlgorithm.minimizationTime;
    _initialSlowDownTime= CAlgorithm.initialSlowDownTime;

    //if run time was not set, look for runsteps parameters
    _runSteps = CAlgorithm.runSteps;
    _snapshotSteps = CAlgorithm.snapshotSteps;
    _minimizationSteps = CAlgorithm.minimizationSteps;
    _neighborListSteps = CAlgorithm.neighborListSteps;

    // create the dissiption tracking object
    _subSystem.pdt = make_unique<DissipationTracker>(simulConfig);
    _cController.initialize(CAlgorithm.algorithm, ChemData, _subSystem.pdt.get(), simulConfig);
    log::info("Done.");



    cout << "---" << endl;
    log::info("Initializing dynamic rates...");
    _drController.initialize(simulConfig);
    log::info("Done.");


    //----------------------------------
    // Initialize elements in the system.
    //----------------------------------

    // Initialize neighbor lists.
    _subSystem.opBoundaryBubbleNL.emplace(simulConfig.boundParams.BoundaryCutoff);
    _subSystem.opBubbleBubbleNL.emplace(simulConfig.mechParams.BubbleCutoff);
    _subSystem.opBubbleBeadNL.emplace(simulConfig.mechParams.BubbleCutoff);

    // Initialize the membrane mesh adapter.
    _meshAdapter = std::make_unique< adaptive_mesh::MembraneMeshAdapter >(simulConfig.geoParams.meshAdapterSettings);
    
    //setup initial network configuration
    setupInitialNetwork(cmdConfig, simulConfig);

    //setup special structures
    setupSpecialStructures(simulConfig);

    //----------------------------------
    // Prepare outputs.
    //----------------------------------
    visual::copySystemMetaData(_subSystem, _mController.ffm, simulConfig, cmdConfig);
    output_.init(cmdConfig.outputDirectory / "traj.h5", simulConfig, _mController.ffm);
    _outputs.push_back(make_unique<BasicSnapshot>((cmdConfig.outputDirectory / "snapshot.traj").string(), &_subSystem));
    _outputs.push_back(make_unique<BirthTimes>((cmdConfig.outputDirectory / "birthtimes.traj").string(), &_subSystem));
    _outputs.push_back(make_unique<Forces>((cmdConfig.outputDirectory / "forces.traj").string(), &_subSystem));
    _outputs.push_back(make_unique<Tensions>((cmdConfig.outputDirectory / "tensions.traj").string(), &_subSystem));

    _outputs.push_back(make_unique<PlusEnd>((cmdConfig.outputDirectory / "plusend.traj").string(), &_subSystem));
    //ReactionOut should be the last one in the output list
    //Otherwise incorrect deltaMinusEnd or deltaPlusEnd values may be genetrated.
    _outputs.push_back(make_unique<ReactionOut>((cmdConfig.outputDirectory / "monomers.traj").string(), &_subSystem));
    //add br force out and local diffussing species concentration
    _outputs.push_back(make_unique<BRForces>((cmdConfig.outputDirectory / "repulsion.traj").string(), &_subSystem));

    //Set up chemistry output if any
    string chemsnapname = (cmdConfig.outputDirectory / "chemistry.traj").string();
    _outputs.push_back(make_unique<Chemistry>(chemsnapname, &_subSystem, ChemData,
                                     _subSystem.getCompartmentGrid()));

    ChemSim* _cs = _cController.getCS();
	ForceFieldManager* _ffm = _mController.getForceFieldManager();

    string concenname = (cmdConfig.outputDirectory / "concentration.traj").string();
    _outputs.push_back(make_unique<Concentrations>(concenname, &_subSystem, ChemData));

    if(chemParams.dissTracking){
        //Set up dissipation output if dissipation tracking is enabled
        string disssnapname = (cmdConfig.outputDirectory / "dissipation.traj").string();
        _outputs.push_back(make_unique<Dissipation>(disssnapname, &_subSystem, _cs));

        //Set up HRCD output if dissipation tracking is enabled
        string hrcdsnapname = (cmdConfig.outputDirectory / "HRCD.traj").string();

        _outputs.push_back(make_unique<HRCD>(hrcdsnapname, &_subSystem, _cs));

        //Set up HRMD output if dissipation tracking is enabled
        string hrmdsnapname = (cmdConfig.outputDirectory / "HRMD.traj").string();
        _outputs.push_back(make_unique<HRMD>(hrmdsnapname, &_subSystem, _cs));

        //Set up MotorWalkingEvents if event tracking is enabled
        string motorwalkingevents = (cmdConfig.outputDirectory / "motorwalkingevents.traj").string();
        _outputs.push_back(make_unique<MotorWalkingEvents>(motorwalkingevents, &_subSystem, _cs));

        //Set up motorunbindingevents if event tracking is enabled
        string motorunbindingevents = (cmdConfig.outputDirectory / "motorunbindingevents.traj").string();
        _outputs.push_back(make_unique<MotorUnbindingEvents>(motorunbindingevents, &_subSystem, _cs));

        //Set up LinkerUnbindingEvents if event tracking is enabled
        string linkerunbindingevents = (cmdConfig.outputDirectory / "linkerunbindingevents.traj").string();
        _outputs.push_back(make_unique<LinkerUnbindingEvents>(linkerunbindingevents, &_subSystem, _cs));

        //Set up LinkerBindingEvents if event tracking is enabled
        string linkerbindingevents = (cmdConfig.outputDirectory / "linkerbindingevents.traj").string();
        _outputs.push_back(make_unique<LinkerBindingEvents>(linkerbindingevents, &_subSystem, _cs));
    }

    if(simulConfig.mechParams.hessTracking){
        if(simulConfig.mechParams.hessMatrixPrintBool){
            //Set up HessianMatrix if hessiantracking is enabled
            string hessianmatrix = (cmdConfig.outputDirectory / "hessianmatrix.traj").string();
            _outputs.push_back(make_unique<HessianMatrix>(hessianmatrix, &_subSystem, _ffm));
        }
        //Set up HessianSpectra if hessiantracking is enabled
        string hessianspectra = (cmdConfig.outputDirectory / "hessianspectra.traj").string();
        _outputs.push_back(make_unique<HessianSpectra>(hessianspectra, &_subSystem, _ffm));

        //Set up Projections if hessiantracking is enabled
        string projections = (cmdConfig.outputDirectory / "projections.traj").string();
        _outputs.push_back(make_unique<Projections>(projections, &_subSystem, _ffm));
    }

    //Set up CMGraph output
    string cmgraphsnapname = (cmdConfig.outputDirectory / "CMGraph.traj").string();
    _outputs.push_back(make_unique<CMGraph>(cmgraphsnapname, &_subSystem));
    
    //Set up TMGraph output
    string tmgraphsnapname = (cmdConfig.outputDirectory / "TMGraph.traj").string();
    _outputs.push_back(make_unique<TMGraph>(tmgraphsnapname, &_subSystem));


    //Set up datadump output if any
    string datadumpname = (cmdConfig.outputDirectory / "datadump.traj").string();
    _outputdump.push_back(make_unique<Datadump>(datadumpname, &_subSystem, ChemData));

    //----------------------------------
    // Finishing initialization.
    //----------------------------------
    SysParams::INITIALIZEDSTATUS = true;
    return simulConfig;
}

void Controller::setupInitialNetwork(const CommandLineConfig& cmdConfig, SimulConfig& simulConfig) {
    using namespace std;

    cout << "---" << endl;
    log::info("Initializing bubbles...");

    //Read bubble setup, parse bubble input file if needed
    auto& BSetup = simulConfig.bubbleSetup;
    auto bubbles = simulConfig.bubbleData;

    //add other bubbles if specified
    auto bubblesGen = createBubblesRandomDist(_subSystem, BSetup.numBubbles, BSetup.bubbleType, simulConfig.mechParams);
    bubbles.bubbles.insert(bubbles.bubbles.end(), bubblesGen.bubbles.begin(), bubblesGen.bubbles.end());

    //add bubbles
    for (auto it: bubbles.bubbles) {

        if(it.type >= simulConfig.mechParams.numBubbleTypes) {
            log::error("Bubble data specified contains an invalid bubble type. Exiting.");
            throw std::runtime_error("Invalid bubble type.");
        }
        SubSystemFunc{}.emplaceTrackable<Bubble>(_subSystem, it.coord, it.type, simulConfig.mechParams);
    }
    log::info("Done. {} bubbles created", bubbles.bubbles.size());


    //--------------------------------------------------------------------------
    // Add MTOC and AFM bubbles.
    //--------------------------------------------------------------------------
    {
        log::info("Initializing MTOC...");
        for(auto& mtocInit : simulConfig.mtocSettings.initVec) {
            auto mtocIndex = SubSystemFunc{}.emplaceTrackable<MTOC>(_subSystem);
            auto& mtoc = _subSystem.mtocs[mtocIndex];

            // Set mechanical parameters.
            mtoc.attachmentStretchingK = mtocInit.attachmentStretchingK;
            mtoc.emiForce1 = mtocInit.emiForce1;
            mtoc.emiForce2 = mtocInit.emiForce2;

            // Create the bubble.
            auto bubbleIndex = SubSystemFunc{}.emplaceTrackable<Bubble>(_subSystem, mtocInit.bubbleCoord, mtocInit.bubbleType, simulConfig.mechParams);
            mtoc.setBubbleSysIndex(_subSystem, bubbleIndex);
            _subSystem.bubbles[bubbleIndex].fixed = mtocInit.bubbleFixed;

            // Create filament attachments.
            auto fils = createFilamentsMTOCDist(
                _subSystem,
                mtocInit,
                simulConfig.mechParams.BubbleRadius[mtocInit.bubbleType],
                simulConfig.geoParams
            ).filaments;

            // Add the filaments.
            for(auto& fil : fils) {
                const auto dist = distance(fil.coords[0], fil.coords[1]);
                const int numSegments = dist / simulConfig.geoParams.cylinderSize[mtocInit.filamentType];

                Filament* pf = _subSystem.addTrackable<Filament>(&_subSystem, mtocInit.filamentType, fil.coords, numSegments + 1, simulConfig, "ARC");
                mtoc.addFilament(pf);
            }

            // Setup chemistry.
            mtoc.setChemistry(_subSystem, mtocInit.vecRxnEmiAbs, simulConfig.chemistryData);
        }
        log::info("{} MTOCs added.", simulConfig.mtocSettings.initVec.size());

        log::info("Initializing AFM...");
        for(auto& afmInit : simulConfig.afmSettings.initVec) {
            auto afmIndex = SubSystemFunc{}.emplaceTrackable<AFM>(_subSystem);
            auto& afm = _subSystem.afms[afmIndex];

            // Set mechanical parameters.
            afm.attachmentStretchingK = afmInit.attachmentStretchingK;
            afm.emiForce1 = afmInit.emiForce1;
            afm.emiForce2 = afmInit.emiForce2;

            // Create the bubble.
            auto bubbleIndex = SubSystemFunc{}.emplaceTrackable<Bubble>(_subSystem, afmInit.bubbleCoord, afmInit.bubbleType, simulConfig.mechParams);
            afm.setBubbleSysIndex(_subSystem, bubbleIndex);
            _subSystem.bubbles[bubbleIndex].fixed = afmInit.bubbleFixed;

            // Create the boundary.
            PlaneBoundaryElement* afmpbe = _subSystem.addTrackable<PlaneBoundaryElement>(
                vec2Vector(afmInit.bubbleCoord),
                vector<floatingpoint>{0, 0, -1},
                simulConfig.boundParams.BoundaryK,
                simulConfig.boundParams.BScreenLength
            );
            afm.setPlaneBoundaryElement(afmpbe);

            // Create filament attachments.
            auto fils = createFilamentsAFMDist(
                _subSystem,
                afmInit,
                simulConfig.mechParams.BubbleRadius[afmInit.bubbleType],
                simulConfig.geoParams
            ).filaments;

            // Add the filaments.
            for(auto& fil : fils) {
                const auto dist = distance(fil.coords[0], fil.coords[1]);
                const int numSegments = dist / simulConfig.geoParams.cylinderSize[afmInit.filamentType];

                Filament* pf = _subSystem.addTrackable<Filament>(&_subSystem, afmInit.filamentType, fil.coords, numSegments + 1, simulConfig, "ARC");
                afm.addFilament(pf);
            }

            // Setup chemistry.
            afm.setChemistry(_subSystem, afmInit.vecRxnEmiAbs, simulConfig.chemistryData);
        }
        log::info("{} AFMs added.", simulConfig.afmSettings.initVec.size());
    }

    /**************************************************************************
    Now starting to add the membrane into the network.
    **************************************************************************/
    cout << "---" << endl;
    log::info("Initializing membranes...");

    const auto& membraneSettings = simulConfig.membraneSettings;
    const auto memChemInfo = MembraneMeshChemistryInfo::fromChemistryData(simulConfig.chemistryData);

    int numMembranes = 0;
    const auto addMembrane = [this, &numMembranes, &memChemInfo](
        const MembraneSetup& memSetup,
        const MembraneInit& memInit,
        const MembraneParser::MembraneInfo& memData
    ) -> void {
        auto newMembraneIndex = SubSystemFunc{}.emplaceTrackable<Membrane>(
            _subSystem,
            // Starting membrane constructor.
            memSetup,
            memData.vertexCoordinateList,
            memData.triangleVertexIndexList
        );
        auto& newMembrane = _subSystem.membranes[newMembraneIndex];

        // Optimize the mesh for membrane
        remeshMembrane(_subSystem, *_meshAdapter, newMembrane);

        // Set up mechanics
        initMechanicParams(_subSystem, newMembrane, memSetup, memInit);

        // Set up chemistry.
        setChemistry(_subSystem, newMembrane.getMesh(), memChemInfo);
        initMeshSpecies(_subSystem, newMembrane.getMesh(), memInit);

        ++numMembranes;
    };

    for(auto& memInit : membraneSettings.initVec) {
        // Find the corresponding membrane setup profile.
        auto& setupVec = membraneSettings.setupVec;
        auto it = find_if(
            setupVec.begin(),
            setupVec.end(),
            [&memInit](const MembraneSetup& setup) { return setup.name == memInit.name; }
        );
        if(it == setupVec.end()) {
            log::error("Membrane setup {} not found.", memInit.name);
            throw std::runtime_error("Membrane setup not found.");
        }

        // Initialize membrane mesh.
        if(memInit.meshParams.size() == 2 && memInit.meshParams[0] == "file") {
            // The input looks like this: mesh file path/to/file
            // Read membrane mesh information from an external file.
            auto memPath = cmdConfig.inputDirectory / std::filesystem::path(memInit.meshParams[1]);
            std::ifstream ifs(memPath);
            if (!ifs.is_open()) {
                log::error("Cannot open membrane file {}", memPath);
                throw std::runtime_error("Cannot open membrane file.");
            }

            const auto memDataVec = MembraneParser::readMembranes(ifs);

            for(auto& memData : memDataVec) {
                addMembrane(*it, memInit, memData);
            }
        }
        else {
            // Forward the input to the membrane mesh initializer
            const auto newMesh = mesh_gen::generateMeshViaParams< floatingpoint >(memInit.meshParams);

            addMembrane(*it, memInit, {newMesh.vertexCoordinateList, newMesh.triangleList});
        }
    }

    log::info("Done. {} membranes created.", numMembranes);

    // Create a region inside the membrane
    log::info("Creating membrane regions...");
    _regionInMembrane = (
        numMembranes == 0 ?
        make_unique<MembraneRegion<Membrane>>(_subSystem.getBoundary()) :
        MembraneRegion<Membrane>::makeByChildren(_subSystem, MembraneHierarchy< Membrane >::root())
    );
    _subSystem.setRegionInMembrane(_regionInMembrane.get());


    log::info("Setting up surface chemistry...");
    {
        setAdsorptionDesorptionReactions(_subSystem, simulConfig.chemistryData.reactionsAdsorptionDesorption);
    }

    log::info("Adjusting compartments by membranes...");

    // Deactivate all the compartments outside membrane, and mark boundaries as interesting
    for(auto& c : _subSystem.getCompartmentGrid()->getCompartments()) {
        if(!c->getTriangles().empty()) {
            // Contains triangles, so this compartment is at the boundary.
            c->boundaryInteresting = true;

            // Update partial activate status
            c->computeSlicedVolumeArea(_subSystem, Compartment::SliceMethod::membrane);
            _cController.updateActivation(*_subSystem.getCompartmentGrid(), c->getId(), Compartment::ActivateReason::Membrane);

        } else if( ! _regionInMembrane->contains(_subSystem, c->coordinates())) {
            // Compartment is outside the membrane
            _cController.deactivate(*_subSystem.getCompartmentGrid(), c->getId(), true);
        }
    }

    // Transfer species from all the inactive compartments
    {
        vector<Compartment*> ac, ic;
        for(auto& c : _subSystem.getCompartmentGrid()->getCompartments()) {
            if(c->isActivated()) ac.push_back(c.get());
            else                 ic.push_back(c.get());
        }
        auto nac = ac.size();
        for(auto c : ic) {
            for(auto &sp : c->getSpeciesContainer().species()) {
                int copyNumber = sp->getN();
                unordered_set<Species*> sp_targets;
                if(sp->getType() == SpeciesType::DIFFUSING) {
                    while(copyNumber > 0) {
                        sp->down();
                        auto tc = ac[Rand::randInteger(0, nac-1)];
                        auto sp_target = tc->findSpeciesByName(sp->getName());
                        sp_targets.insert(sp_target);
                        sp_target->up();
                        --copyNumber;
                    }
                }
                for(auto sp_target : sp_targets)
                    sp_target->updateReactantPropensities();
                sp->updateReactantPropensities();
            }
        }
    }

    log::info("Membrane initialization complete.");

    //----------------------------------
    // Add fixed vertex attachments.
    //----------------------------------
    log::info("Adding fixed vertex attachments...");
    {
        const auto numAdded = createFixedVertexAttachments<SubSystemFunc>(_subSystem, membraneSettings.fixedVertexAttachmentInits);
        log::info("Added {} fixed vertex attachments.", numAdded);
    }

    /**************************************************************************
    Now starting to add the filaments into the network.
    **************************************************************************/
    //Read filament setup, parse filament input file if needed
    auto& FSetup = simulConfig.filamentSetup;
    
    cout << "---" << endl;
    log::info("Initializing filaments...");

    if (SysParams::RUNSTATE == true) {
        auto fil = simulConfig.filamentData.filaments;
        //add other filaments if specified

        auto filGen = createFilamentsRandomDist(
            _subSystem,
            *_regionInMembrane,
            FSetup.numFilaments,
            FSetup.filamentType,
            FSetup.filamentLength,
            FSetup
        ).filaments;
        fil.insert(fil.end(), filGen.begin(), filGen.end());

        //add filaments

        for (auto& it: fil) {

            auto type = it.type;

            if (type >= simulConfig.chemParams.numFilaments) {
                log::error("Filament data specified contains an invalid filament type {}.", type);
                throw std::runtime_error("Filament data specified contains an invalid filament type.");
            }

            if(it.coords.size() < 2) {
                log::error("Filament data specified contains insufficient coordinates.");
                throw std::runtime_error("Filament data specified contains insufficient coordinates.");
            }
            else if (it.coords.size() == 2) {

                const FP d = distance(it.coords[0], it.coords[1]);
                const int numSegment = std::max<int>(
                    1,
                    static_cast<int>(std::round(d / simulConfig.geoParams.cylinderSize[type]))
                );

                // check how many segments can fit between end-to-end of the filament
                if (numSegment == 0)
                    _subSystem.addTrackable<Filament>(&_subSystem, type, it.coords, 2, simulConfig,
                                                        FSetup.projectionType);
                else
                    _subSystem.addTrackable<Filament>(&_subSystem, type, it.coords, numSegment + 1, simulConfig,
                                                        FSetup.projectionType);
            }
            else {
                const int numSegment = it.coords.size() - 1;

                _subSystem.addTrackable<Filament>(&_subSystem, type, it.coords, numSegment + 1, simulConfig,
                                                        FSetup.projectionType);
            }
        }
        cout << "Done. " << fil.size() << " filaments created." << endl;
        cout << "Total cylinders " << Cylinder::getCylinders().size() << endl;
    }
    else{
        cout<<endl;
	    cout<<"RESTART PHASE BEINGS."<<endl;
        //Create the restart pointer
        const string inputfileName = (cmdConfig.inputDirectory / FSetup.inputFile).string();
        _restart = new Restart(&_subSystem, simulConfig.chemistryData, inputfileName);
        //read set up.
        _restart->readNetworkSetup();
        _restart->setupInitialNetwork();
    }
}

void Controller::setupSpecialStructures(SimulConfig& simulConfig) {

    cout << "---" << endl;
    cout << "Setting up special structures...";

    cout << "Done." << endl;
}

void Controller::activatedeactivateComp(){

    auto& grid = *_subSystem.getCompartmentGrid();
    if(SysParams::Boundaries().transfershareaxis>=0){
        fCompmap.clear();
        bCompmap.clear();
        activatecompartments.clear();
        ControlfrontbackEndComp();
//            std::cout<<fCompmap.size()<<" "<<bCompmap.size()<<" "<<activatecompartments.size()<<endl;
        for(auto it=activatecompartments.begin();it!=activatecompartments.end();it++)
        {
            if(!(*it)->isActivated())
                _cController.activate(grid, (*it)->getId());
        }
        //deactivate compartments starting from the right extreme
        for (std::multimap<int,Compartment*>::reverse_iterator it=fCompmap.rbegin(); it!=fCompmap.rend(); ++it)
            _cController.deactivate(grid, it->second->getId());
        //deactivate compartments starting from the left extreme
        for (std::multimap<int,Compartment*>::iterator it=bCompmap.begin(); it!=bCompmap.end(); ++it)
            _cController.deactivate(grid, it->second->getId());
        fCompmap.clear();
        bCompmap.clear();

    }
}

void Controller::ControlfrontbackEndComp(){
    auto& grid = *_subSystem.getCompartmentGrid();

    Compartment* maxcomp=NULL;
    Compartment* mincomp=NULL;
    int tsaxis = SysParams::Boundaries().transfershareaxis;
    int planestomove = SysParams::Boundaries().planestomove;
    bool maxcompstate = false;
    bool mincompstate = false;
    if(planestomove == 2 || planestomove == 0) maxcompstate = true;
    if(planestomove == 2 || planestomove == 1) mincompstate = true;
    floatingpoint systemspan = 0.0;
    floatingpoint cmpsize = 0.0;
    if(tsaxis == 0)
    {systemspan = SysParams::Geometry().NX * SysParams::Geometry()
                                                              .compartmentSizeX;
    cmpsize = SysParams::Geometry().compartmentSizeX;}
    else if(tsaxis == 1) {
        systemspan = SysParams::Geometry().NY * SysParams::Geometry()
                .compartmentSizeY;
        cmpsize = SysParams::Geometry().compartmentSizeY;
    }
    else if(tsaxis == 2) {
        systemspan = SysParams::Geometry().NZ * SysParams::Geometry().compartmentSizeZ;
        cmpsize = SysParams::Geometry().compartmentSizeZ;
    }
    //copy vector to prevcopy
    bounds_prev[0] = bounds[0];bounds_prev[1] = bounds[1];
    bounds[0] = 0.0; bounds[1] =  systemspan;
    for(auto& C : _subSystem.getCompartmentGrid()->getCompartments()){
        auto cyls=C->getCylinders();
        if(cyls.size()>0){
            //maxcomp refers to the compartment on the right extreme of reaction volume
            // that represents the current chemical boundary of the system.
            if(maxcompstate) {
                if (maxcomp == NULL)
                    maxcomp = C.get();
                else {
                    //get current maxcomp coordinates
                    auto mcoord = maxcomp->coordinates();
                    //get compartment coorinates
                    auto ccord = C->coordinates();
                    //compare to see if the compartment is further to the right of maxcomp.
                    if (mcoord[SysParams::Boundaries().transfershareaxis] <
                        ccord[SysParams::Boundaries()
                                .transfershareaxis])
                        maxcomp = C.get();
                }
            }
            //mincomp refers to the compartment on the left extreme of reaction volume
            // that represents the current chemical boundary of the system.
            if(mincompstate) {
                if (mincomp == NULL)
                    mincomp = C.get();
                else {
                    auto mcoord = mincomp->coordinates();
                    auto ccord = C->coordinates();
                    //compare to see if the compartment is further to the left of mincomp.
                    if (mcoord[SysParams::Boundaries().transfershareaxis] >
                        ccord[SysParams::Boundaries().transfershareaxis])
                        mincomp = C.get();
                }
            }
        }
    }

    if(maxcompstate) {
        std::cout<<"1 maxcomp "<<maxcomp->coordinates()[0]<<endl;
        // front end is defined two compartments away from the current maxcomp.
        auto cmaxcomp = maxcomp->coordinates();
        //get the neighbor who is to the right of maxcomp.
        for (auto cnindex : maxcomp->getNeighborIndices()) if(cnindex != -1) {
            auto C = &grid.getCompartment(cnindex);
            auto cC = C->coordinates();
            if (cmaxcomp[SysParams::Boundaries().transfershareaxis] <
                cC[SysParams::Boundaries().transfershareaxis])
                maxcomp = C;
        }
        std::cout<<"2 maxcomp "<<maxcomp->coordinates()[0]<<endl;
        cmaxcomp = maxcomp->coordinates();
        //get the neighbor who is to the right of maxcomp.
        for (auto cnindex : maxcomp->getNeighborIndices()) if(cnindex != -1) {
            auto C = &grid.getCompartment(cnindex);
            auto cC = C->coordinates();
            if (cmaxcomp[SysParams::Boundaries().transfershareaxis] <
                cC[SysParams::Boundaries().transfershareaxis])
                maxcomp = C;
        }
        std::cout<<"3 maxcomp "<<maxcomp->coordinates()[0]<<endl;
        cmaxcomp = maxcomp->coordinates();
        assert((maxcomp != NULL) && "Non existent maxcomp. Exiting.");
        //Loop through compartments
        for (auto& C : _subSystem.getCompartmentGrid()->getCompartments()) {
            auto cC = C->coordinates();
            //if compartment is to the right of maxcomp and activated, add to a vector to
            // deactivate later.
            if (cC[SysParams::Boundaries().transfershareaxis] >
                cmaxcomp[SysParams::Boundaries().transfershareaxis]) {
                if (C->isActivated())
                    fCompmap.insert(pair<int, Compartment *>(
                            cC[SysParams::Boundaries().transfershareaxis], C.get()));
            }
                //if compartment is to the left of maxcomp and not active, add to a vector to
                // activate later.
            else {
                if (!(C->isActivated()))
                    activatecompartments.push_back(C.get());
            }
        }
        bounds[1] = maxcomp->coordinates()[SysParams::Boundaries().transfershareaxis] +
                cmpsize/2;
    }
    //back end is defined as the compartment that is two compartments to the left of
    // mincomp.
    if(mincompstate) {
        auto cmincomp = mincomp->coordinates();
        //get the neighbor who is to the left of mincomp.
        for (auto cnindex : mincomp->getNeighborIndices()) if(cnindex != -1) {
            auto C = &grid.getCompartment(cnindex);
            auto cC = C->coordinates();
            if (cmincomp[SysParams::Boundaries().transfershareaxis] >
                cC[SysParams::Boundaries().transfershareaxis])
                mincomp = C;
        }
        cmincomp = mincomp->coordinates();
        //get the neighbor who is to the left of mincomp.
        for (auto cnindex : mincomp->getNeighborIndices()) if(cnindex != -1) {
            auto C = &grid.getCompartment(cnindex);
            auto cC = C->coordinates();
            if (cmincomp[SysParams::Boundaries().transfershareaxis] >
                cC[SysParams::Boundaries().transfershareaxis])
                mincomp = C;
        }
        cmincomp = mincomp->coordinates();
        assert(mincomp != NULL && "Non existent mincomp. Exiting.");
        //Loop through compartments
        for (auto& C : _subSystem.getCompartmentGrid()->getCompartments()) {
            auto cC = C->coordinates();
            //if compartment C is to the left of mincomp and was added to
            // activatecompartments vector, remove. If it is already active, add to a vector
            // to deactivate later.
            if (cC[SysParams::Boundaries().transfershareaxis] <
                cmincomp[SysParams::Boundaries().transfershareaxis]) {
                auto it = std::find(activatecompartments.begin(),
                                    activatecompartments.end(), C.get());
                if (it != activatecompartments.end())
                    activatecompartments.erase(it);
                if (C->isActivated()) {
                    bCompmap.insert(pair<int, Compartment *>(
                            cC[SysParams::Boundaries().transfershareaxis], C.get()));
                }
            }
        }
        bounds[0] = mincomp->coordinates()[SysParams::Boundaries().transfershareaxis] -
                cmpsize/2;
    }
    //print the maximum (right boundary) and minimum (left boundary) compartment spans.
    std::cout<<"Maxbound "<<bounds[1]<<" Minbound "<<bounds[0]<<endl;
}

void Controller::moveBoundary(floatingpoint deltaTau) {
    auto& grid = *_subSystem.getCompartmentGrid();

    //calculate distance to move
    floatingpoint dist = SysParams::Boundaries().moveSpeed * deltaTau;

    if(SysParams::Boundaries().transfershareaxis>=0){
        vector<floatingpoint> distvec= {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        if(SysParams::Boundaries().transfershareaxis == 0){
            distvec[0] = bounds[0] - bounds_prev[0];
            distvec[1] = bounds[1] - bounds_prev[1];
        }
        else if(SysParams::Boundaries().transfershareaxis == 1){
            distvec[2] = bounds[0] - bounds_prev[0];
            distvec[3] = bounds[1] - bounds_prev[1];
        }
        else if(SysParams::Boundaries().transfershareaxis == 2){
            distvec[4] = bounds[0] - bounds_prev[0];
            distvec[5] = bounds[1] - bounds_prev[1];
        }
        _subSystem.getBoundary()->move(distvec);
    }
        //deprecated not good to use.
    else if(abs(dist)>0){
        vector<floatingpoint> distvec = {dist, -dist, dist, -dist, dist, -dist};
        //move it
        if(tau() >= SysParams::Boundaries().moveStartTime &&
           tau() <= SysParams::Boundaries().moveEndTime)
            _subSystem.getBoundary()->move(distvec);

        //activate, deactivate necessary compartments
        for(auto& C : grid.getCompartments()) {

            if(_subSystem.getBoundary()->within(C.get())) {

                if(C->isActivated()) continue;
                else _cController.activate(grid, C->getId());
            }
            else {
                if(!C->isActivated()) continue;
                else _cController.deactivate(grid, C->getId());
            }
        }
    }
    //calculate system volume.
    _subSystem.getBoundary()->volume();
}

void Controller::updateActiveCompartments() {
    // For this function to work, we must assume that each minimization step
    // will push the membrane boundary no more than 1 compartment, so that
    // changes will only happen at the neighborhood of the previous boundary.
    auto& allMembranes = _subSystem.membranes;
    auto& grid = *_subSystem.getCompartmentGrid();

    // Currently only the 0th membrane will be considered
    if(allMembranes.size()) {
        auto& theMembrane = *allMembranes.begin();
        // For non empty compartments, we mark them as interesting and update their status
        // For the "interesting" compartments last round but now empty, we fully activate or deactivate them
        // For the rest we do nothing, assuming that the membranes will NOT move across a whole compartment
        for(auto& c: grid.getCompartments()) {
            const auto& ts = c->getTriangles();
            if(!ts.empty()) {
                // Update partial activate status
                c->computeSlicedVolumeArea(_subSystem, Compartment::SliceMethod::membrane);
                _cController.updateActivation(grid, c->getId(), Compartment::ActivateReason::Membrane);

                // No matter whether the compartment is interesting before, mark it as interesting
                c->boundaryInteresting = true;
            } else if(c->boundaryInteresting) { // Interesting last round but now empty
                bool inMembrane = (
                    (!theMembrane.isClosed()) ||
                    (medyan::contains(_subSystem, theMembrane.getMesh(), c->coordinates()))
                );
                if(inMembrane) {
                    // Fully activate the compartment
                    c->resetVolumeFrac();
					const auto& fullArea = grid.compartmentAreas;
                    c->setPartialArea({{
                        fullArea[0], fullArea[0],
						fullArea[1], fullArea[1],
						fullArea[2], fullArea[2]
                    }});
                    _cController.updateActivation(grid, c->getId(), Compartment::ActivateReason::Membrane);
                } else {
                    // Deactivate the compartment
                    _cController.deactivate(grid, c->getId());
                }

                // Mark the compartment as not interesting
                c->boundaryInteresting = false;
            }
        }
    } // Otherwise, no membrane exists. Do nothing
}

void Controller::initializeSpecialProtocols(const SimulConfig& conf) {
    // Making filaments static.
    if(conf.chemParams.makeFilamentsStatic) {
        specialProtocolManager_.addProtocol(
            SpecialProtocol::Schedule::after,
            conf.chemParams.makeFilamentsStaticTime,
            [](SubSystem& sys, const SimulConfig& conf) {
                for(auto pc : Cylinder::getCylinders()) {
                    pc->getCCylinder()->passivatefilreactions();
                }
            }
        );
    }

    // Making linkers static.
    if(conf.chemParams.makeLinkersStatic) {
        specialProtocolManager_.addProtocol(
            SpecialProtocol::Schedule::after,
            conf.chemParams.makeLinkersStaticTime,
            [](SubSystem& sys, const SimulConfig& conf) {
                for(auto pl : Linker::getLinkers()) {
                    pl->getCLinker()->getOffReaction()->passivateReaction();
                }
            }
        );
    }

    // Pinning boundary filaments.
    if(conf.mechParams.pinBoundaryFilaments) {
        specialProtocolManager_.addProtocol(
            SpecialProtocol::Schedule::after,
            conf.mechParams.pinTime,
            [](SubSystem& sys, const SimulConfig& conf) {
                //if we've already added pinned filaments, return
                if(Bead::getPinnedBeads().size() != 0)
                    return;

                //loop through beads, check if within pindistance
                for(auto b : Bead::getBeads()) {

                    //pin only beads who are at the front of a plus end cylinder or back of a minus end cylinder
                    Filament* f = (Filament*) b->getParent();
                    Cylinder* plusEndC = f->getPlusEndCylinder();
                    Cylinder* minusEndC = f->getMinusEndCylinder();

                    if((plusEndC->getSecondBead() == b) ||
                    (minusEndC->getFirstBead() == b)) {

                        cout << sys.getBoundary()->distance(b->vcoordinate()) << endl;
                        cout << conf.mechParams.pinDistance << endl;


                        //if within dist to boundary, add
                        if(sys.getBoundary()->distance(b->vcoordinate()) < conf.mechParams.pinDistance) {

                            b->pinnedPosition = b->vcoordinate();
                            b->addAsPinned();
                        }
                    }
                }
            }
        );
    }

    // Pinning lower boundary filaments.
    if(conf.mechParams.pinLowerBoundaryFilaments) {
        specialProtocolManager_.addProtocol(
            SpecialProtocol::Schedule::after,
            conf.mechParams.pinTime,
            [](SubSystem& sys, const SimulConfig& conf) {
                //renew pinned filament list everytime

                //loop through beads, check if within pindistance
                for(auto b : Bead::getBeads()) {

                    //pin all beads besides plus end and minus end cylinder
                    Filament* f = (Filament*) b->getParent();
                    Cylinder* plusEndC = f->getPlusEndCylinder();
                    Cylinder* minusEndC = f->getMinusEndCylinder();

                    if((plusEndC->getSecondBead() != b) ||
                    (minusEndC->getFirstBead() != b)) {

                        auto index = Rand::randfloatingpoint(0,1);
                        //cout << index <<endl;
                        //if within dist to boundary and index > 0.5, add
                        if(sys.getBoundary()->lowerdistance(b->vcoordinate()) < conf.mechParams.pinDistance
                        && index < conf.mechParams.pinFraction && b->isPinned() == false) {
                            //cout << index << endl;
                            b->pinnedPosition = b->vcoordinate();
                            b->addAsPinned();
                        }
                    }
                }
            }
        );
    }

    // Auxiliary function to pin filaments in a certain region.
    // inRegion: takes 3D coordinates and returns whether the point is in the region.
    const auto pinInitialFilamentWith = [](SubSystem& sys, auto&& inRegion) {
        for(auto b : Bead::getBeads()) {
            if(inRegion(b->coordinate())) {
                b->pinnedPosition = b->vcoordinate();
                b->addAsPinned();
            }
        }
    };

    // Pinning filaments below Z.
    if(conf.mechParams.pinInitialFilamentBelowZ) {
        specialProtocolManager_.addProtocol(
            SpecialProtocol::Schedule::once,
            0.0,
            [pinInitialFilamentWith](SubSystem& sys, const SimulConfig& conf) {
                pinInitialFilamentWith(sys, [&conf](auto&& c) { return c[2] < conf.mechParams.pinInitialFilamentBelowZValue; });
            }
        );
    }

    // Membrane equilibrium area change.
    if(conf.mechParams.membraneEqAreaChange.has_value()) {
        specialProtocolManager_.addProtocol(
            SpecialProtocol::Schedule::after,
            0.0,
            [](SubSystem& sys, const SimulConfig& conf) {
                auto& change = conf.mechParams.membraneEqAreaChange.value();
                const auto inc = change.rate * conf.chemParams.chemistryAlgorithm.minimizationTime;
                for(auto& m : sys.membranes) {
                    auto& mesh = m.getMesh();
                    if(!mesh.metaAttribute().hasLipidReservoir) {
                        // Find old equilibrium area.
                        floatingpoint oldEqArea = 0;
                        for(auto& v : mesh.getVertices()) {
                            auto& objv = v.attr.vertex(sys);
                            oldEqArea += objv.mVertex.eqArea;
                        }
                        // Increase equilibrium area.
                        const auto newEqArea = std::max<FP>(change.minEqArea, oldEqArea + inc);
                        const auto ratio = newEqArea / oldEqArea;
                        // Redistribute equilibrium area.
                        for(auto& v : mesh.getVertices()) {
                            auto& objv = v.attr.vertex(sys);
                            objv.mVertex.eqArea *= ratio;
                        }
                    }
                }
            }
        );
    }
    // Membrane equilibrium volume change.
    if(conf.mechParams.membraneEqVolumeChange.has_value()) {
        specialProtocolManager_.addProtocol(
            SpecialProtocol::Schedule::after,
            0.0,
            [](SubSystem& sys, const SimulConfig& conf) {
                auto& change = conf.mechParams.membraneEqVolumeChange.value();
                const auto inc = change.rate * conf.chemParams.chemistryAlgorithm.minimizationTime;
                for(auto& m : sys.membranes) {
                    m.mMembrane.eqVolume = std::max<FP>(
                        change.minEqVolume,
                        m.mMembrane.eqVolume + inc
                    );
                    log::info("Membrane eqVolume changed to {}.", m.mMembrane.eqVolume);
                }
            }
        );
    }
}

void Controller::updatePositions(double& tp) {
    //NEED TO UPDATE CYLINDERS FIRST
    //Reset Cylinder update position state
    Cylinder::setpositionupdatedstate = false;
    for(auto c : Cylinder::getCylinders()) {
        c->updatePosition();
    }

    //Reset state to updated state
    Cylinder::setpositionupdatedstate = true;

    // Update all other movables.
    //---------------------------------
    SubSystemFunc{}.forEachMovable(_subSystem, movable::updatePosition);

    // Update bubble positions.
    const auto updateBubblePositions = [this, &tp] {
        for(auto& b : _subSystem.bubbles) {
            if(b.isAFM()) b.updatePositionManually(_subSystem);
        }
    
        if(SysParams::Chemistry().makeRateDepend && tau() - tp > 1) {
            tp+=1;
            
            for(auto &filament : Filament::getFilaments()) {
                double deltaL;
                double numCyl = 0;
                for (auto cylinder : filament->getCylinderVector()){
                    
                    deltaL += cylinder->getMCylinder()->getLength() -
                    cylinder->getMCylinder()->getEqLength();
                    numCyl += 1;
                }
                
                //print last
                Cylinder* cylinder = filament->getCylinderVector().back();
                deltaL += cylinder->getMCylinder()->getLength() -
                cylinder->getMCylinder()->getEqLength();
                numCyl += 1;
                
                double k = cylinder->getMCylinder()->getStretchingConst();
                
                //if the filament tension is higher than threshold, regardless of sign
                if(k*deltaL/numCyl > SysParams::Chemistry().makeRateDependForce ||
                -k*deltaL/numCyl > SysParams::Chemistry().makeRateDependForce ){
                    
                    Cylinder* pCyl = filament->getCylinderVector().back();
                    for(auto &r : pCyl->getCCylinder()->getInternalReactions()) {
                        if(r->getReactionType() == ReactionType::POLYMERIZATIONPLUSEND) {
                            float newrate = 5 * SysParams::Chemistry().originalPolyPlusRate;
                            r->setBareRate(newrate);
                            r->recalcRateVolumeFactor();
                            r->updatePropensity();
                        }
                    }
                }
                //else, set it back to orginal rate
                else{
                    Cylinder* pCyl = filament->getCylinderVector().back();
                    for(auto &r : pCyl->getCCylinder()->getInternalReactions()) {
                        if(r->getReactionType() == ReactionType::POLYMERIZATIONPLUSEND) {
                            float newrate = SysParams::Chemistry().originalPolyPlusRate;
                            r->setBareRate(newrate);
                            r->recalcRateVolumeFactor();
                            r->updatePropensity();
                        }
                    }
                }
            }
        }
    };
    if(SysParams::Chemistry().makeAFM) updateBubblePositions();

}


void Controller::updateReactionRates(const SimulConfig& conf) {
    SubSystemFunc{}.forEachReactable(_subSystem, reactable::updateReactionRates);

    // Update adsorption/desorption reaction rates.
    setAdsorptionDesorptionReactionRates(_subSystem, conf.chemistryData);
}

void Controller::updateNeighborLists(const CommandLineConfig& cmdConfig, SimulConfig& conf) {
    #ifdef CROSSCHECK_CYLINDER
        if(HybridNeighborList::_crosscheckdumpFileNL.is_open())
            HybridNeighborList::_crosscheckdumpFileNL.close();
        string crosscheckNLname = (cmdConfig.outputDirectory / "crosscheckNL.traj").string();
        HybridNeighborList::_crosscheckdumpFileNL.open(crosscheckNLname);
    #endif
    chrono::high_resolution_clock::time_point mins, mine;

    mins = chrono::high_resolution_clock::now();
    //Full reset of neighbor lists
    _subSystem.resetNeighborLists();
    mine = chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_runnl2(mine - mins);
    nl2time += elapsed_runnl2.count();

    mins = chrono::high_resolution_clock::now();
    _subSystem.updateBindingManagers(conf);
    #ifdef OPTIMOUT
        cout<<"updated BindingManagers"<<endl;
    #endif
    mine = chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_runb(mine - mins);
    bmgrtime += elapsed_runb.count();
}

void Controller::resetCounters() {
    for(Filament* f : Filament::getFilaments()) f->resetCounters();
}


void Controller::membraneAdaptiveRemesh() {
    // Requires _meshAdapter to be already initialized
    for(auto& m : _subSystem.membranes) {
        _meshAdapter->adapt(_subSystem, m.getMesh());

        // Update necessary geometry for the system
        medyan::updateGeometryValueForSystem(_subSystem, m.getMesh());
    }

    for(auto& t : _subSystem.triangles) {
        movable::updatePosition(_subSystem, t);
    }
}

void Controller::run(const CommandLineConfig& cmdConfig, SimulConfig& conf) {

    // Auxiliary functions.
    const auto copyVisualData = [&, this] {
        using namespace visual;
        using namespace visual::raw_data_cat;
        copySystemData(
            _subSystem, conf, cmdConfig,
            beadPosition | beadConnection | compartment | concentration
        );
    };

    const auto cgInterruptFunc = [&, this] {
        if(_subSystem.membranes.size() > 0) {
            membraneAdaptiveRemesh();
            _subSystem.resetNeighborLists();
            copyVisualData();
        }
    };
    const int cgInterruptNumIter = 50;
    const int cgInterruptCallLimit = (conf.mechParams.mechanicsAlgorithm.adaptMeshDuringMinimization && _subSystem.membranes.size() > 0)
        ? 200 : 0;
    const CGInterruptSettings cgInterruptSettings { cgInterruptNumIter, cgInterruptCallLimit, cgInterruptFunc };

    long totalSteps = 0;

    chrono::high_resolution_clock::time_point chk1, chk2, mins, mine;
    chk1 = chrono::high_resolution_clock::now();
    chrono::high_resolution_clock::time_point minsR, mineR;
//RESTART PHASE BEGINS
    if(SysParams::RUNSTATE==false){

        minsR = chrono::high_resolution_clock::now();
//Step 2A. Turn off diffusion, passivate filament reactions and add reactions to heap.
        _restart->settorestartphase();
	    cout<<"Turned off Diffusion, and filament reactions."<<endl;
        cout<<"Bound species added to reaction heap."<<endl;

//Step 3. ############ RUN LINKER/MOTOR REACTIONS TO BIND BRANCHERS, LINKERS, MOTORS AT RESPECTIVE POSITIONS.#######
        cout<<"Number of reactions to be fired "<<_restart->getnumchemsteps()<<endl;
        _cController.runSteps(_restart->getnumchemsteps(), conf.chemistryData, conf);
        cout<<"Reactions fired! Displaying number of reactions that are NOT fired in each"
              " compartment"
              ""<<endl;
//Step 4. Display the number of reactions yet to be fired. Should be zero.
        bool exitstatus = 0;
        for(auto& C : _subSystem.getCompartmentGrid()->getCompartments()) {
            for(auto &Mgr:C->getFilamentBindingManagers()){
                int numsites = 0;
#ifdef NLORIGINAL
                numsites = Mgr->numBindingSites();
#else
                numsites = Mgr->numBindingSitesstencil();
#endif
                if(numsites == 0)
                    cout<< numsites<<" ";
                else{
                    cout<<endl;
                    log::error("Compartment COORDS {} {} {}", C->coordinates()[0], C->coordinates()[1], C->coordinates()[2]);
                    log::error("Num binding sites {}", numsites);
                    string mgrname ="";
                    if(dynamic_cast<BranchingManager*>(Mgr.get()))
                        mgrname = " BRANCHING ";
                    else if (dynamic_cast<LinkerBindingManager*>(Mgr.get()))
                        mgrname = " LINKER ";
                    else
                        mgrname = " MOTOR ";
                    log::error("Printing {} binding sites that were not chosen", mgrname);
                    #ifdef NLORIGINAL
                    Mgr->printbindingsites();
					#else
                	Mgr->printbindingsitesstencil();
					#endif
                	exitstatus = true;
                }
            }}
        cout<<endl;
        if(exitstatus) {
            cout << "Few reactions were not fired! Cannot restart this trajectory. "
                    "Exiting after printing diffusing species in each compartment..." <<
                    endl;

            cout<< "COMPARTMENT DATA: CMPID DIFFUSINGSPECIES COPYNUM"<<endl;
            for(auto& cmp:_subSystem.getCompartmentGrid()->getCompartments()){
                cout <<cmp->getId()<<" ";
                for(auto sd : conf.chemistryData.speciesDiffusing) {
                    string name = sd.name;
                    auto s = cmp->findSpeciesByName(name);
                    auto copyNum = s->getN();
                    cout <<name<<" "<<copyNum<<" ";
                }
                cout <<endl;
            }
            exit(EXIT_FAILURE);
        }
///STEP 5. Reset time to required restart time.
        _cController.initializerestart(_restart->getrestartime(),_minimizationTime, conf.chemistryData);
	    _initialSlowDownTime += _restart->getrestartime();
        cout<<"Tau reset to "<<tau()<<endl;
///STEP 6. Reinitialize CBound eqlen, numHeads and numBoundHeads values as required by
/// datadump
        //sets
        _restart->CBoundinitializerestart();
///STEP 7. Assign copynumbers based on Chemistry input file or the datadump file as
// required by the user.
        _restart->restartupdateCopyNumbers();
        _restart->crosscheck();
        cout<<"Diffusion rates restored, diffusing molecules redistributed."<<endl;


        //Step 8. re-add pin positions
        {
            auto simulConfig = getSimulConfigFromInput(cmdConfig.inputFile, cmdConfig.inputDirectory);
            auto& filSetup = simulConfig.filamentSetup;

            if(SysParams::Mechanics().pinBoundaryFilaments){
                PinRestartParser ppin((cmdConfig.inputDirectory / filSetup.pinRestartFile).string());
                ppin.resetPins();
            }
        }

//Step 9. run mcontroller, update system, turn off restart state.
        updatePositions(tp);
        updateNeighborLists(cmdConfig, conf);

        mins = chrono::high_resolution_clock::now();
        cout<<"Minimizing energy"<<endl;

        _subSystem.prevMinResult = _mController.run(conf, cgInterruptSettings);
#ifdef OPTIMOUT
        mine= chrono::high_resolution_clock::now();
        chrono::duration<floatingpoint> elapsed_runm(mine - mins);
        minimizationtime += elapsed_runm.count();
        std::cout<<"Time taken for minimization "<<elapsed_runm.count()<<endl;
#endif
        //DO NOT MOVE THIS LINE
        SysParams::RUNSTATE=true;

        //reupdate positions and neighbor lists
        updatePositions(tp);
        updateNeighborLists(cmdConfig, conf);

//Step 10. Set Off rates back to original value.
        for(auto LL : Linker::getLinkers())
        {
            LL->getCLinker()->setOffRate(LL->getCLinker()->getOffReaction()->getBareRate());
            LL->updateReactionRates();
            LL->getCLinker()->getOffReaction()->updatePropensity();
            /*cout<<"L "<<LL->getId()<<" "<<LL->getMLinker()->getEqLength()<<" "
                << LL->getCLinker()->getOffRate()<<" "
                <<LL->getCLinker()->getOffReaction()->getBareRate()<<" "
                    <<LL->getMLinker()->stretchForce<<endl;*/

        }
        for(auto MM : MotorGhost::getMotorGhosts())
        {
            MM->getCMotorGhost()->setOffRate(MM->getCMotorGhost()->getOffReaction()->getBareRate());
            MM->updateReactionRates();
            MM->getCMotorGhost()->getOffReaction()->updatePropensity();
            /*cout<<"M "<<MM->getId()<<" "<<MM->getMMotorGhost()->getEqLength()<<" "
                << MM->getCMotorGhost()->getOffRate()<<" "
                <<MM->getCMotorGhost()->getOffReaction()->getBareRate()<<" "
                <<MM->getMMotorGhost()->stretchForce<<endl;*/
        }
        int dummy=0;
        for (auto BB: BranchingPoint::getBranchingPoints()) {
            dummy++;
            BB->getCBranchingPoint()->setOffRate(BB->getCBranchingPoint()->getOffReaction()->getBareRate());
            BB->updateReactionRates();
            BB->getCBranchingPoint()->getOffReaction()->updatePropensity();
            /*cout<<"B "<<BB->getId()<<" "<<BB->getMBranchingPoint()->getEqLength()<<" "
                << BB->getCBranchingPoint()->getOffRate()<<" "
                <<BB->getCBranchingPoint()->getOffReaction()->getBareRate()<<endl;*/
        }
//STEP 11: Get cylinders, activate filament reactions.
        for(auto& C : _subSystem.getCompartmentGrid()->getCompartments()) {
            for(auto x : C->getCylinders()) {
                x->getCCylinder()->activatefilreactions();
                x->getCCylinder()->activatefilcrossreactions();
            }}

//Step 11b. Activate general reactions.
        for(auto& C : _subSystem.getCompartmentGrid()->getCompartments()) {
            for(auto& rxn : C->getInternalReactionContainer().reactions()) {
                if(rxn->getReactionType() == ReactionType::REGULAR)
                    rxn->activateReaction();
            }}
        cout<<"Unbinding rates of bound species restored. filament reactions activated"<<endl;
//@

        _subSystem.updateBindingManagers(conf);

        updateReactionRates(conf);

        delete _restart;
        cout<< "Restart procedures completed. Starting original Medyan framework"<<endl;
        cout << "---" << endl;
        log::info("Current simulation time = {}", tau());
        cout << endl;
        //restart phase ends

        //Crosscheck tau to make sure heap is ordered accurately.
        _cController.crosschecktau();
        mineR = chrono::high_resolution_clock::now();
    }
    chrono::duration<floatingpoint> elapsed_runRestart(mineR-minsR);

    cout<<"Minimizing energy"<<endl;
    mins = chrono::high_resolution_clock::now();
    _subSystem.resetNeighborLists(); // TODO: resolve workaround
    Bead::rearrange();
    Cylinder::updateAllData();
    // update neighorLists before and after minimization. Need excluded volume
    // interactions.
	_subSystem.resetNeighborLists();
    copyVisualData();

    // Initial special protocols need to be executed before energy minimization
    executeSpecialProtocols(conf);
    auto minimizationResult = _mController.run(conf, cgInterruptSettings);
    membraneAdaptiveRemesh();
    _subSystem.resetNeighborLists();
    _mController.updateMechanics(conf);
    copyVisualData();
    _subSystem.prevMinResult = minimizationResult;
    mine= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_runm2(mine - mins);
    minimizationtime += elapsed_runm2.count();
    #ifdef OPTIMOUT
        std::cout<<"Time taken for minimization "<<elapsed_runm2.count()<<endl;
    #endif

    //activate/deactivate compartments
    mins = chrono::high_resolution_clock::now();
    //set initial values of variables.
    int tsaxis = SysParams::Boundaries().transfershareaxis;
    floatingpoint systemspan = 0.0;
    if(tsaxis == 0)
        systemspan = SysParams::Geometry().NX * SysParams::Geometry().compartmentSizeX;
    else if(tsaxis == 1)
        systemspan = SysParams::Geometry().NY * SysParams::Geometry().compartmentSizeY;
    else if(tsaxis == 2)
        systemspan = SysParams::Geometry().NZ * SysParams::Geometry().compartmentSizeZ;
    //copy vector to prevcopy
    bounds_prev[1] = systemspan;bounds_prev[0] = 0.0;
    bounds[1] = systemspan; bounds[0] =  0.0;
    activatedeactivateComp();
    moveBoundary(0.0);
    mine= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_runspl(mine - mins);
    specialtime += elapsed_runspl.count();

    //reupdate positions and neighbor lists
    mins = chrono::high_resolution_clock::now();
    updatePositions(tp);
    #ifdef OPTIMOUT
        cout<<"Positions updated"<<endl;
    #endif
    updateNeighborLists(cmdConfig, conf);
    #ifdef OPTIMOUT
        mine= chrono::high_resolution_clock::now();
        chrono::duration<floatingpoint> elapsed_runnl(mine - mins);
        nltime += elapsed_runnl.count();
        std::cout<<"NL time "<<elapsed_runnl.count()<<endl;
        mins = chrono::high_resolution_clock::now();
    #endif

    updateReactionRates(conf);
    mine= chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_runrxn(mine - mins);
    rxnratetime += elapsed_runrxn.count();

    for(auto& o: _outputs) o->print(0, conf);
    for(auto& o: _outputdump) o->print(0);
    output_.append(_subSystem, conf);

    resetCounters();

    int snapshotCounter = 1;

    // Should the simulation proceed with certain amount of time or certain number of chemical reactions?
    // true:  run with fixed time interval, until the given total time.
    // false: run with fixed number of chemical reactions, until the given total number of reactions.
    const bool useTime = _runTime != 0;

    log::info("Starting simulation in {} mode", useTime ? "time" : "step");

    {
        // Preparations before the first chemical reaction.
        //---------------------------------------------------------------------

        //number of mechanical minimizations per Snapshot and 
        //neighbor list updates 
        const int minsPerSnapshot=
            useTime
            ? max(1L,std::lround(conf.chemParams.chemistryAlgorithm.snapshotTime/_minimizationTime))
            : _minimizationSteps / _snapshotSteps;
        const int minsPerNeighborList=
            useTime
            ? max(1L,std::lround(conf.chemParams.chemistryAlgorithm.neighborListTime/_minimizationTime))
            : _minimizationSteps / _neighborListSteps;
        const int minsPerDatadump=
            useTime
            ? max(1L,std::lround(conf.chemParams.chemistryAlgorithm.datadumpTime/_minimizationTime))
            : 200; // Temporary fixed value.
        int64_t minimizationCounter= 0;

        //activate/deactivate compartments
        mins = chrono::high_resolution_clock::now();
        activatedeactivateComp();
	    // set initial mechanical energy of system through a call to force field manager if dissipation tracking is enabled
	    if(conf.chemParams.dissTracking){
		    _subSystem.pdt->setG1(minimizationResult.energiesAfter);
	    }
        mine= chrono::high_resolution_clock::now();
        chrono::duration<floatingpoint> elapsed_runspl(mine - mins);
        specialtime += elapsed_runspl.count();

        while(
            (useTime  && tau() < _runTime) ||
            (!useTime && totalSteps < _runSteps)
        ) {
            // Record some variables before executing the loop.
            //-----------------------------------------------------------------
            const auto oldTau = tau();

            // Run chemical controller.
            //-----------------------------------------------------------------
            #ifdef OPTIMOUT
                cout<<"Starting chemistry"<<endl;
			#endif
            SysParams::DURINGCHEMISTRY = true;
            mins = chrono::high_resolution_clock::now();
            float factor = 1.0;
            if(tau() <_initialSlowDownTime)
            	factor = 10.0;

            #ifdef CROSSCHECK_CYLINDER
                string crosscheckchemname = (cmdConfig.outputDirectory / "crosscheckChem.traj").string();
                if(CController::_crosscheckdumpFilechem.is_open())
                    CController::_crosscheckdumpFilechem.close();
                CController::_crosscheckdumpFilechem.open(crosscheckchemname);
            #endif
            const auto chemSuccess =
                useTime
                ? _cController.run(_minimizationTime / factor, conf.chemistryData, conf)
                : _cController.runSteps(_minimizationSteps, conf.chemistryData, conf);
            mine= chrono::high_resolution_clock::now();
            chrono::duration<floatingpoint> elapsed_runchem(mine - mins);
            chemistrytime += elapsed_runchem.count();
            SysParams::DURINGCHEMISTRY = false;
            #ifdef OPTIMOUT
                auto mtimex = CUDAcommon::tmin;
                cout<<"motorbinding calls "<<mtimex.motorbindingcalls<<endl;
                cout<<"motorunbinding calls "<<mtimex.motorunbindingcalls<<endl;
                cout<<"motorwalking calls "<<mtimex.motorwalkingcalls<<endl;
                cout<<"linkerbinding calls "<<mtimex.linkerbindingcalls<<endl;
                cout<<"linkerunbinding calls "<<mtimex.linkerunbindingcalls<<endl;
                CUDAcommon::tmin.motorbindingcalls = 0;
                CUDAcommon::tmin.motorunbindingcalls = 0;
                CUDAcommon::tmin.motorwalkingcalls = 0;
                CUDAcommon::tmin.linkerbindingcalls = 0;
                CUDAcommon::tmin.linkerunbindingcalls = 0;
            #endif
            //print output if chemistry fails.
            mins = chrono::high_resolution_clock::now();
            if(!chemSuccess) {
                for(auto& o: _outputs) { o->print(snapshotCounter, conf); }
                output_.append(_subSystem, conf);
                resetCounters();
                break;
            }

            mine= chrono::high_resolution_clock::now();
            chrono::duration<floatingpoint> elapsed_runout(mine - mins);
            outputtime += elapsed_runout.count();
            // Update with the last step. Note that the time has been automatically updated.
            if(!useTime) {
                totalSteps += _minimizationSteps;
            }

            // Run mechanical controller, and update system.
            //-----------------------------------------------------------------
            {
	            //check before rearrange
                #ifdef CROSSCHECK_CYLINDER
                    string crosscheckname = (cmdConfig.outputDirectory / "crosscheckcyl.traj").string();
                    Cylinder::_crosscheckdumpFile.open(crosscheckname);
                    if(!Cylinder::_crosscheckdumpFile.is_open()) {
                        Cylinder::_crosscheckdumpFile << "There was an error opening file " << crosscheckname
                            << " for output. Exiting." << endl;
                        exit(EXIT_FAILURE);
                    }
                    Cylinder::_crosscheckdumpFile << "Opening file " << crosscheckname << endl;
                    Cylinder::_crosscheckdumpFile << "NCylinders " << Cylinder::getCylinders().size() << endl;
                #endif
                mins = chrono::high_resolution_clock::now();

                Bead::rearrange();
                Cylinder::updateAllData();
                // Update neighbor lists for mechanical interactions for energy minimization.
                _subSystem.resetNeighborLists();

                string crosscheckmechname = (cmdConfig.outputDirectory / "crosscheckmech.traj").string();
                CGMethod::_crosscheckdumpMechFile.open(crosscheckmechname);

                minimizationResult = _mController.run(conf, cgInterruptSettings);
                _subSystem.prevMinResult = minimizationResult;
                // Membrane remeshing after mechanical energy minimization.
                membraneAdaptiveRemesh();
                _subSystem.resetNeighborLists();
                _mController.updateMechanics(conf);
                copyVisualData();
                #ifdef CROSSCHECK_CYLINDER
                    CGMethod::_crosscheckdumpMechFile.close();
                #endif
                mine= chrono::high_resolution_clock::now();

                
                chrono::duration<floatingpoint> elapsed_runm3(mine - mins);
                minimizationtime += elapsed_runm3.count();
                #ifdef OPTIMOUT
                    std::cout<<"Time taken for minimization "<<elapsed_runm3.count()<<endl;
				#endif

                //update position
                mins = chrono::high_resolution_clock::now();

                updatePositions(tp);

                // Update activation of the compartments
                updateActiveCompartments();
                #ifdef CROSSCHECK_CYLINDER
                    Cylinder::_crosscheckdumpFile.close();
                #endif

                #ifdef OPTIMOUT
                    cout<<"Position updated"<<endl;
				#endif

                mine= chrono::high_resolution_clock::now();
                chrono::duration<floatingpoint> elapsed_rxn2(mine - mins);
                updateposition += elapsed_rxn2.count();

                // perform multiple functions to update cumulative energy counters and reset the mechanical energy variables
                if(conf.chemParams.dissTracking){
                    _subSystem.pdt->setGMid(minimizationResult.energiesBefore);
                    _subSystem.pdt->updateAfterMinimization(minimizationResult.energiesAfter);
                }

	            //update reaction rates
	            mins = chrono::high_resolution_clock::now();
	            updateReactionRates(conf);
                #ifdef OPTIMOUT
	                cout<<"updated Reaction Rates"<<endl;
                #endif
	            mine= chrono::high_resolution_clock::now();
	            chrono::duration<floatingpoint> elapsed_rxn3(mine - mins);
	            rxnratetime += elapsed_rxn3.count();
                minimizationCounter++;
            }
            //End of minimization

            // Output snapshot.
            //-----------------------------------------------------------------
            if(minimizationCounter%minsPerSnapshot == 0) {
                mins = chrono::high_resolution_clock::now();
                for(auto& o: _outputs) o->print(snapshotCounter, conf);
                output_.append(_subSystem, conf);
                resetCounters();
                mine= chrono::high_resolution_clock::now();
                chrono::duration<floatingpoint> elapsed_runout2(mine - mins);
                outputtime += elapsed_runout2.count();

                if(conf.outputParams.logBriefProfilingEachSnapshot) {
                    log::info("Total time spent until cycle {}:", minimizationCounter);
                    log::info("- Chemistry:            {}", chemistrytime);
                    log::info("- Minimization:         {}", minimizationtime);
                    log::info("- Neighbor list update: {}", nltime);
                    log::info("- Position update:      {}", updateposition);
                    log::info("- Reaction rate update: {}", rxnratetime);
                    log::info("- Output:               {}", outputtime);
                    log::info("- Special protocols:    {}", specialtime);
                }
                ++snapshotCounter;
            }
            if(minimizationCounter%minsPerDatadump == 0) {
                for (auto& o: _outputdump) o->print(0);
            }
            if(minimizationCounter%minsPerSnapshot == 0 || conf.outputParams.logSimulationTimeEachCycle) {
                log::info(
                    "Current simulation time = {}{}{}",
                    tau(),
                    (minimizationCounter % minsPerSnapshot == 0) ? " (snapshot)" : "",
                    (minimizationCounter % minsPerDatadump == 0) ? " (dump)" : ""
                );
            }

            // Other system updates, including neighbor list update.
            //-----------------------------------------------------------------
            //activate/deactivate compartments
            mins = chrono::high_resolution_clock::now();
            activatedeactivateComp();
            //move the boundary
            moveBoundary(tau() - oldTau);
            mine= chrono::high_resolution_clock::now();
            chrono::duration<floatingpoint> elapsed_runspl(mine - mins);
            specialtime += elapsed_runspl.count();

            // update neighbor lists & Binding Managers
            if(minimizationCounter%minsPerNeighborList == 0) {
                mins = chrono::high_resolution_clock::now();
                updateNeighborLists(cmdConfig, conf);
                mine= chrono::high_resolution_clock::now();
                chrono::duration<floatingpoint> elapsed_runnl2(mine - mins);
                nltime += elapsed_runnl2.count();
                #ifdef OPTIMOUT
                cout<<"update NeighborLists"<<endl;
                #endif
            }
            //Special protocols
            mins = chrono::high_resolution_clock::now();
            //special protocols
            executeSpecialProtocols(conf);
            mine= chrono::high_resolution_clock::now();
            chrono::duration<floatingpoint> elapsed_runspl2(mine - mins);
            specialtime += elapsed_runspl2.count();

        }
    }

    //print last snapshots
    for(auto& o: _outputs) o->print(snapshotCounter, conf);
    output_.append(_subSystem, conf);
    
    
    
    
    // rockingsnapshot with last snapshot
    if(conf.mechParams.rockSnapBool){
    //Set up RockingSnapshot if hessiantracking is enabled
        ForceFieldManager* _ffm =  _mController.getForceFieldManager();
        Eigen::VectorXcd evalues = _ffm->evalues;
        for(auto k = 0; k < evalues.size(); k++){

            string rockingsnaphot = (cmdConfig.outputDirectory / ("rockingsnapshot_" + to_string(evalues.real()[k]) +".traj")).string();
            _rSnapShot = new RockingSnapshot(rockingsnaphot, &_subSystem, _ffm, k);
            _rSnapShot->savePositions();
            _rSnapShot->print(snapshotCounter);
            _rSnapShot->resetPositions();
            _rSnapShot->~RockingSnapshot();
         
        }
           
        
    };

    output_.finish();
	resetCounters();
    chk2 = chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_run(chk2-chk1);
    cout << "Time elapsed for run: dt=" << elapsed_run.count() << endl;
	#ifdef OPTIMOUT
    cout<<"Restart time for run=" << elapsed_runRestart.count()<<endl;
    cout<< "Chemistry time for run=" << chemistrytime <<endl;
    cout << "Minimization time for run=" << minimizationtime <<endl;
    cout<< "Neighbor-list+Bmgr-time for run="<<nltime<<endl;
    cout<< "Neighbor-list time for run="<<nl2time<<endl;
    cout<< "Bmgr-vec time for run="<<bmgrvectime<<endl;
    cout<< "SIMD time for run="<<SubSystem::SIMDtime<<endl;
    cout<< "HYBD time for run="<<SubSystem::HYBDtime<<endl;
    cout<< "Bmgr time for run="<<bmgrtime<<endl;
    cout<<"update-position time for run="<<updateposition<<endl;
    cout<<"rxnrate time for run="<<rxnratetime<<endl;
    cout<<"Output time for run="<<outputtime<<endl;
    cout<<"Special time for run="<<specialtime<<endl;
    cout << "Time elapsed for run: dt=" << elapsed_run.count() << endl;
    cout << "Total simulation time: dt=" << tau() << endl;
    cout<<"-----------"<<endl;
    if(true){
        auto mtime = CUDAcommon::tmin;
        cout<<"update-position time for run="<<updateposition<<endl;
        cout<<"update-position-cylinder time for run="<<updatepositioncylinder<<endl;
        cout<<"update-position-movable time for run="<<updatepositionmovable<<endl;
        cout<<"move-compartment cylinder ="<<mtime.timecylinderupdate<<" calls "<<mtime
        .callscylinderupdate<<endl;
        cout<<"move-compartment linker ="<<mtime.timelinkerupdate<<" calls "<<mtime
                .callslinkerupdate<<endl;
        cout<<"move-compartment motor ="<<mtime.timemotorupdate<<" calls "<<mtime
                .callsmotorupdate<<endl;
        auto cdetails = CUDAcommon::cdetails;
        cout<<"Clone internal reactions="<<cdetails.ccylclonetimer[0]<<", calls="<<
                cdetails.ccylclonecounter[0]<<", rxncounter="<<cdetails.ccylclonerxncounter[0]<<endl;
        cout<<"Clone rxn alone="<<cdetails.internalrxnclone<<endl;
        cout<<"Add cloned reaction="<<cdetails.internalrxnadd<<endl;
        cout<<"Find species to clone="<<cdetails.clonefindspecies<<endl;
        cout<<"Get affected reactions="<<cdetails.getaffectedrxns<<endl;
        cout<<"Clone crossCylinder reactions="<<cdetails.ccylclonetimer[1]<<", calls="<<
            cdetails.ccylclonecounter[1]<<", rxncounter="<<cdetails.ccylclonerxncounter[1]<<endl;
        cout<<"Clone reactingCylinder reactions="<<cdetails.ccylclonetimer[2]<<", calls="<<
            cdetails.ccylclonecounter[2]<<", rxncounter="<<cdetails.ccylclonerxncounter[2]<<endl;
        cout<<"-----------"<<endl;

        cout << "Minimization time for run=" << minimizationtime <<endl;
        cout<<"Printing minimization times in seconds."<<endl;
        cout<<"starting minimization "<<mtime.vectorize<<endl;
        cout<<"Finding lambda "<<mtime.findlambda<<endl;
        cout<<"copy forces "<<mtime.copyforces<<endl;
        cout<<"other computations "<<mtime.tother<<endl;
        cout<<"end minimization "<<mtime.endminimization<<endl;
        cout<<"compute energies "<<mtime.computeenergy<<" calls "<<mtime
        .computeenergycalls<<endl;
        cout<<"compute energieszero "<<mtime.computeenergyzero<<" calls "<<mtime
        .computeenerycallszero<<endl;
	    cout<<"compute energiesnonzero "<<mtime.computeenergynonzero<<" calls "<<mtime
			    .computeenerycallsnonzero<<endl;
        cout<<"compute forces "<<mtime.computeforces<<" calls "<<mtime
                .computeforcescalls<<endl;
        cout<<"Time taken to compute energy in each forcefield "<<endl;
        cout<<"Filament, Linker, Motor, Branching, Excluded Volume, and Boundary"<<endl;
        for(auto x:mtime.individualenergies)
            cout<<x<<" ";
        cout<<endl;
	    cout<<"Time taken to compute energy in "<<endl;
	    cout<<"Filament Stretching "<<mtime.stretchingenergy<<endl;
	    cout<<"Filament Bending "<<mtime.bendingenergy<<endl;
        cout<<"Time taken to compute energyzero in each forcefield "<<endl;
        for(auto x:mtime.individualenergieszero)
            cout<<x<<" ";
        cout<<endl;
        cout<<"Time taken to compute energynonzero in each forcefield "<<endl;
        for(auto x:mtime.individualenergiesnonzero)
            cout<<x<<" ";
        cout<<endl;
        cout<<"Time taken to compute forces in each forcefield "<<endl;
        for(auto x:mtime.individualforces)
            cout<<x<<" ";
        cout<<endl;
	    cout<<"Time taken to compute forces in "<<endl;
	    cout<<"Filament Stretching "<<mtime.stretchingforces<<endl;
	    cout<<"Filament Bending "<<mtime.bendingforces<<endl;

		cout<<"Number of interactions considered in each force field"<<endl;

		cout<<"Filament Stretching "<<mtime.numinteractions[0]<<endl;
	    cout<<"Filament Bending "<<mtime.numinteractions[1]<<endl;
	    cout<<"Linker Stretching "<<mtime.numinteractions[2]<<endl;
	    cout<<"Motor Stretching "<<mtime.numinteractions[3]<<endl;
	    cout<<"Branching Stretching "<<mtime.numinteractions[4]<<endl;
	    cout<<"Branching Bending "<<mtime.numinteractions[5]<<endl;
	    cout<<"Branching Dihedral "<<mtime.numinteractions[6]<<endl;
	    cout<<"Branching Position "<<mtime.numinteractions[7]<<endl;
	    cout<<"Cylinder-Cylinder Repulsion "<<mtime.numinteractions[8]<<endl;
	    cout<<"Cylinder-Boundary Repulsion "<<mtime.numinteractions[9]<<endl;
    }
    if(true) {
        cout << "Printing callback times" << endl;
        auto ctime = CUDAcommon::ctime;
        auto ccount = CUDAcommon::ccount;
        cout << "UpdateBrancherBindingCallback " << ctime.tUpdateBrancherBindingCallback
             << " count "
             << ccount.cUpdateBrancherBindingCallback << endl;
        cout << "UpdateLinkerBindingCallback " << ctime.tUpdateLinkerBindingCallback
             << " count "
             << ccount.cUpdateLinkerBindingCallback << endl;
        cout << "UpdateMotorBindingCallback " << ctime.tUpdateMotorBindingCallback
             << " count "
             << ccount.cUpdateMotorBindingCallback << endl;
        cout << "UpdateMotorIDCallback " << ctime.tUpdateMotorIDCallback << " count "
             << ccount.cUpdateMotorIDCallback << endl;
        cout << "FilamentExtensionPlusEndCallback "
             << ctime.tFilamentExtensionPlusEndCallback << " count "
             << ccount.cFilamentExtensionPlusEndCallback << endl;
        cout << "FilamentExtensionMinusEndCallback "
             << ctime.tFilamentExtensionMinusEndCallback << " count "
             << ccount.cFilamentExtensionMinusEndCallback << endl;
        cout << "FilamentRetractionPlusEndCallback "
             << ctime.tFilamentRetractionPlusEndCallback << " count "
             << ccount.cFilamentRetractionPlusEndCallback << endl;
        cout << "FilamentRetractionMinusEndCallback "
             << ctime.tFilamentRetractionMinusEndCallback << " count "
             << ccount.cFilamentRetractionMinusEndCallback << endl;
        cout << "FilamentPolymerizationPlusEndCallback "
             << ctime.tFilamentPolymerizationPlusEndCallback << " count "
             << ccount.cFilamentPolymerizationPlusEndCallback << endl;
        cout << "FilamentPolymerizationMinusEndCallback "
             << ctime.tFilamentPolymerizationMinusEndCallback << " count "
             << ccount.cFilamentPolymerizationMinusEndCallback << endl;
        cout << "FilamentDepolymerizationPlusEndCallback "
             << ctime.tFilamentDepolymerizationPlusEndCallback << " count "
             << ccount.cFilamentDepolymerizationPlusEndCallback << endl;
        cout << "FilamentDepolymerizationMinusEndCallback "
             << ctime.tFilamentDepolymerizationMinusEndCallback << " count "
             << ccount.cFilamentDepolymerizationMinusEndCallback << endl;
        cout << "BranchingPointUnbindingCallback " << ctime.tBranchingPointUnbindingCallback
             << " count "
             << ccount.cBranchingPointUnbindingCallback << endl;
        cout << "BranchingCallback " << ctime.tBranchingCallback << " count "
             << ccount.cBranchingCallback << endl;
        cout << "LinkerUnbindingCallback " << ctime.tLinkerUnbindingCallback << " count "
             << ccount.cLinkerUnbindingCallback << endl;
        cout << "LinkerBindingCallback " << ctime.tLinkerBindingCallback << " count "
             << ccount.cLinkerBindingCallback << endl;
        cout << "MotorUnbindingCallback " << ctime.tMotorUnbindingCallback << " count "
             << ccount.cMotorUnbindingCallback << endl;
        cout << "MotorBindingCallback " << ctime.tMotorBindingCallback << " count "
             << ccount.cMotorBindingCallback << endl;
        cout << "MotorWalkingCallback " << ctime.tMotorWalkingCallback << " count "
             << ccount.cMotorWalkingCallback << endl;
        cout << "MotorMovingCylinderCallback " << ctime.tMotorMovingCylinderCallback
             << " count "
             << ccount.cMotorMovingCylinderCallback << endl;
        cout << "FilamentCreationCallback " << ctime.tFilamentCreationCallback << " count "
             << ccount.cFilamentCreationCallback << endl;
        cout << "FilamentSeveringCallback " << ctime.tFilamentSeveringCallback << " count "
             << ccount.cFilamentSeveringCallback << endl;
        cout << "FilamentDestructionCallback " << ctime.tFilamentDestructionCallback
             << " count "
             << ccount.cFilamentDestructionCallback << endl;
        cout << "------------" << endl;
        cout << "Printing neighbor times" << endl;
        cout << "Dynamic neighbor " << SubSystem::timedneighbor << endl;
        cout << "Neighbor " << SubSystem::timeneighbor << endl;
        cout << "Trackable " << SubSystem::timetrackable << endl;

        cout << "-------------" << endl;
        cout << "Filament extendPlusEnd 1 " << Filament::FilextendPlusendtimer1 << endl;
        cout << "Filament extendPlusEnd 2 " << Filament::FilextendPlusendtimer2 << endl;
        cout << "-------------" << endl;
        cout << "Cylinder constructor" << endl;
        cout << "part1 " << Cylinder::timecylinder1 << " part2 " << Cylinder::timecylinder2
             << " "
                "Ccylinder "
             << Cylinder::timecylinderchem << " mCylinder " << Cylinder::timecylindermech
             << endl;
        cout << "initializeCCylinder for loop " << ChemManager::tchemmanager1 << endl;
        cout << "extension Front/Back " << ChemManager::tchemmanager2 << endl;
        cout << "initialize " << ChemManager::tchemmanager3 << endl;
        cout << "last part " << ChemManager::tchemmanager4 << endl;
        cout << "------------" << endl;
        cout << "PolyPlusEndTemplate time" << endl;
        cout << "For loop " << CUDAcommon::ppendtime.rxntempate1 << " part2 (findspecies) "
             << CUDAcommon::ppendtime.rxntempate2 << " part3 (create rxn) "
             << CUDAcommon::ppendtime
                     .rxntempate3 << " part4 (Callback) "
             << CUDAcommon::ppendtime.rxntempate4 << endl;
        cout<<" Displaying chemistry times"<<endl;
        cout<<"Counts fired for each ReactionType"<<endl;
        for(auto i = 0; i<17;i++)
            cout<<CUDAcommon::cdetails.reactioncount[i]<<" ";
        cout<<endl;
        cout<<"Time taken to fire each ReactionType"<<endl;
        for(auto i = 0; i<17;i++)
            cout<<CUDAcommon::cdetails.totaltime[i]<<" ";
        cout<<endl;
        cout<<"Time taken to emitSignal for each ReactionType"<<endl;
        for(auto i = 0; i<17;i++)
            cout<<CUDAcommon::cdetails.emitsignal[i]<<" ";
        cout<<endl;
        cout<<"Time taken for dependency updates for each ReactionType"<<endl;
        for(auto i = 0; i<17;i++)
            cout<<CUDAcommon::cdetails.dependencytime[i]<<" ";
        cout<<endl;
        cout<<"Total number of dependencies for reactions fired based on ReactionType"<<endl;
        for(auto i = 0; i<17;i++)
            cout<<CUDAcommon::cdetails.dependentrxncount[i]<<" ";
        cout<<endl;
        cout<<"Diffusion passivate vs activate calls"<<endl;
        cout<<CUDAcommon::cdetails.diffusion_passivate_count<<" "<<CUDAcommon::cdetails.diffusion_activate_count<<endl;
    }
	#endif
    log::info("Done with simulation!");
}

} // namespace medyan

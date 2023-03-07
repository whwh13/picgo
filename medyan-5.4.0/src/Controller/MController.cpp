
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

#include "Controller/MController.h"

#include "SubSystem.h"

#include "LinkerFF.h"
#include "MotorGhostFF.h"
#include "BranchingFF.h"
#include "BubbleFF.h"
#include "Mechanics/ForceField/Boundary/BoundaryFF.h"
#include "Mechanics/ForceField/Filament/FilamentFF.h"
#include "Mechanics/ForceField/Membrane/MembraneFF.hpp"
#include "Mechanics/ForceField/Volume/CylinderExclVolume.h"
#include "Mechanics/ForceField/Volume/CylinderVolumeMon.hpp"
#include "Mechanics/ForceField/Volume/TriangleBeadVolumeFF.hpp"
#include "Mechanics/ForceField/VolumeConservation/VolConsrvFF.hpp"
#include "Mechanics/Minimizer/CGMethod.hpp"
#include "Util/Io/Log.hpp"

namespace medyan {

void MController::initializeMinAlgorithm (const MechParams::MechanicsAlgorithm& MAlgorithm) {


    cgParams = medyan::ConjugateGradientParams {
        [&] {
            if (MAlgorithm.ConjugateGradient == "FLETCHERRIEVES") {
                return medyan::ConjugateGradientDescentSearch::fletcherRieves;
            }
            else if (MAlgorithm.ConjugateGradient == "POLAKRIBIERE") {
                return medyan::ConjugateGradientDescentSearch::polakRibiere;
            }
            else if (MAlgorithm.ConjugateGradient == "STEEPESTDESCENT") {
                return medyan::ConjugateGradientDescentSearch::steepest;
            }
            else {
                log::error("Conjugate gradient method not recognized. Exiting.");
                throw std::runtime_error("Unrecognized conjugate gradient method.");
            }
        }(),
        MAlgorithm.maxDistance,
        MAlgorithm.lambdaMax,
        MAlgorithm.lambdarunningaverageprobability,
        MAlgorithm.linesearchalgorithm,
        MAlgorithm.gradientTolerance,
        MAlgorithm.energyChangeRelativeTolerance,
        0, 10000, true,
        MAlgorithm.tryToRecoverInLineSearchError,
    };

}

void MController::initializeFF (const SimulConfig& conf) {
    auto& mechParams = conf.mechParams;
    auto& ffTypes = mechParams.mechanicsFFType;

    // Initialize all force fields.
    {
        auto filamentFFs = createFilamentForceFields(
            ffTypes.FStretchingType,
            ffTypes.FBendingType,
            ffTypes.FTwistingType
        );
        for(auto& pff : filamentFFs) ffm.forceFields.push_back(std::move(pff));
    }

    {
        auto allLinkerFFs = createLinkerForceFields(
            ffTypes.LStretchingType,
            ffTypes.LBendingType,
            ffTypes.LTwistingType
        );
        for(auto& ff : allLinkerFFs) ffm.forceFields.push_back(std::move(ff));
    }
    {
        auto allMotorFFs = createMotorForceFields(
            ffTypes.MStretchingType,
            ffTypes.MBendingType,
            ffTypes.MTwistingType
        );
        for(auto& ff : allMotorFFs) ffm.forceFields.push_back(std::move(ff));
    }
    {
        auto allBranchingFFs = createBranchingForceFields(
            ffTypes.BrStretchingType,
            ffTypes.BrBendingType,
            ffTypes.BrDihedralType,
            ffTypes.BrPositionType
        );
        for(auto& ff : allBranchingFFs) ffm.forceFields.push_back(std::move(ff));
    }

    {
        auto membraneFFRes = createMembraneForceFields(
            ffTypes.memStretchingFFType,
            ffTypes.memTensionFFType,
            ffTypes.memBendingFFType,
            conf.geoParams.surfaceGeometrySettings
        );
        for(auto& ff : membraneFFRes) ffm.forceFields.push_back(std::move(ff));
        ffm.surfaceGeometrySettings = conf.geoParams.surfaceGeometrySettings;
    }

    {
        auto volConsrvFFRes = VolumeConservationFFFactory {}(
            ffTypes.volumeConservationFFType
        );
        for(auto& ff : volConsrvFFRes) ffm.forceFields.push_back(std::move(ff));
    }
    
    //These FF's have a neighbor list associated with them
    //add to the subsystem's database of neighbor lists.
    //--------------------------------------------------------------------------

    // Cylinder volume exclusion potentials.
    if(ffTypes.cylinderVolumeExclusionFFType == CylinderVolumeExclusionFFType::integral) {
        auto volumeFF = std::make_unique<CylinderExclVolume<CylinderExclVolRepulsion>>(conf);

        #if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
            //Original neighborlist is created even when HybridNeighborList is use as the
            // neighborLists are implemented as predominantly modifications in the back-end.
            // Front end functions remain the same and function calls are internally
            // redirected to HybridNeighborList functions.
            volumeFF->setHNeighborLists(_subSystem->getHNeighborList());
        #endif

        //Get the force field access to the HNLID
        for(auto nl : volumeFF->getNeighborLists()) {
            if(nl != nullptr)
                _subSystem->addNeighborList(nl);
        }
        ffm.forceFields.push_back(std::move(volumeFF));
    }
    else if(ffTypes.cylinderVolumeExclusionFFType == CylinderVolumeExclusionFFType::monomer) {
        auto pff = std::make_unique<CylinderVolumeMon>(_subSystem->getHNeighborList(), conf);

        for(auto pnl : pff->getNeighborLists()) {
            if(pnl != nullptr) {
                _subSystem->addNeighborList(pnl);
            }
        }
        ffm.forceFields.push_back(std::move(pff));
    }

    {
        auto allTriBeadVolFF = TriangleBeadVolumeFFFactory{}(
            *_subSystem,
            ffTypes.triangleBeadVolumeFFType
        );
        for(auto& ff : allTriBeadVolFF) {
            auto& ffp = ffm.forceFields.emplace_back(std::move(ff));
            for(auto nl : ffp->getNeighborLists()) {
                if(nl) _subSystem->addNeighborList(nl);
            }
        }
    }

    {
        auto allBoundaryFF = createBoundaryForceFields(ffTypes.BoundaryFFType, conf);
        for(auto& ff : allBoundaryFF) {
            auto& ffp = ffm.forceFields.emplace_back(std::move(ff));
            for(auto nl : ffp->getNeighborLists()) {
                if(nl) _subSystem->addNeighborList(nl);
            }
        }
    }

    {
        auto allBubbleFF = createBubbleForceFields(ffTypes.BubbleFFType, ffTypes.MTOCFFType, ffTypes.AFMFFType);
        for(auto& ff : allBubbleFF) {
            auto& ffp = ffm.forceFields.emplace_back(std::move(ff));
            for(auto nl : ffp->getNeighborLists()) {
                if(nl) _subSystem->addNeighborList(nl);
            }
        }
    }


    // Print all the force field names.
    log::debug("System is initialized with the following force fields:");
    for(auto& pff : ffm.forceFields) {
        log::debug("- {}", pff->getName());
    }
}

} // namespace medyan

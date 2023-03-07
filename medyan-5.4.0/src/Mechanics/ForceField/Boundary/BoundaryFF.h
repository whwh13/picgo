
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

#ifndef MEDYAN_BoundaryFF_h
#define MEDYAN_BoundaryFF_h

#include <vector>

#include "common.h"
#include "Mechanics/ForceField/Boundary/BoundaryCylinderRepulsion.h"
#include "Mechanics/ForceField/Boundary/BoundaryCylinderRepulsionExp.h"
#include "Mechanics/ForceField/Boundary/BoundaryCylinderRepulsionIn.h"
#include "Mechanics/ForceField/Boundary/BoundaryCylinderRepulsionExpIn.h"
#include "Mechanics/ForceField/Boundary/BoundaryBubbleRepulsion.h"
#include "Mechanics/ForceField/Boundary/BoundaryBubbleRepulsionExp.h"
#include "Mechanics/ForceField/Boundary/BoundaryCylinderAttachment.h"
#include "Mechanics/ForceField/Boundary/BoundaryCylinderAttachmentHarmonic.h"
#include "Mechanics/ForceField/ForceField.h"

namespace medyan {

inline auto createBoundaryForceFields(std::string_view type, const SimulConfig& conf) {
    std::vector<std::unique_ptr<ForceField>> forceFields;

    if (type == "REPULSIONEXP") {
        forceFields.push_back(
            std::make_unique<BoundaryCylinderRepulsion<BoundaryCylinderRepulsionExp>>(conf)
        );
        forceFields.push_back(
            std::make_unique<BoundaryBubbleRepulsion<BoundaryBubbleRepulsionExp>>()
        );
    }
    else if(type == "REPULSIONEXPIN") {
        forceFields.push_back(
            std::make_unique<BoundaryCylinderRepulsionIn<BoundaryCylinderRepulsionExpIn>>(conf)
        );
        forceFields.push_back(
            std::make_unique<BoundaryBubbleRepulsion<BoundaryBubbleRepulsionExp>>()
        );
    }
    else if(type == "") {
        log::warn("No boundary FF is specified.");
    }
    else {
        log::error("Boundary FF {} not recognized. Exiting.", type);
        throw std::runtime_error("Boundary FF not recognized.");
    }

    //if pinning to boundaries
    if(conf.mechParams.pinBoundaryFilaments ||
        conf.mechParams.pinInitialFilamentBelowZ ||
        conf.mechParams.pinLowerBoundaryFilaments
    ) {
        forceFields.push_back(
            std::make_unique<BoundaryCylinderAttachment<BoundaryCylinderAttachmentHarmonic>>()
        );
    }

    return forceFields;
}

} // namespace medyan

#endif

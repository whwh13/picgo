
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

#ifndef MEDYAN_BubbleFF_h
#define MEDYAN_BubbleFF_h

#include <string_view>
#include <vector>

#include "common.h"
#include "Mechanics/ForceField/Bubble/AFMAttachment.h"
#include "Mechanics/ForceField/Bubble/BubbleBubbleRepulsion.h"
#include "Mechanics/ForceField/Bubble/BubbleCylinderRepulsion.h"
#include "Mechanics/ForceField/Bubble/FixedBubbleCoordinates.hpp"
#include "Mechanics/ForceField/Bubble/MTOCAttachment.h"
#include "Mechanics/ForceField/Bubble/MTOCBending.h"

namespace medyan {

inline auto createBubbleForceFields(std::string_view type, std::string_view mtoc, std::string_view afm) {
    std::vector<std::unique_ptr<ForceField>> forceFields;

    // Fixed bead coordinates are always present.
    forceFields.push_back(std::make_unique<FixedBubbleCoordinates>());

    if(type == "REPULSIONEXP") {
        forceFields.push_back(std::make_unique<BubbleCylinderRepulsion>());
        forceFields.push_back(std::make_unique<BubbleBubbleRepulsion>());
    }
    else if(type == "") {}
    else {
        log::error("Bubble FF not recognized. Exiting.");
        throw std::runtime_error("Bubble FF not recognized.");
    }

    // MTOC specific.
    if(mtoc == "ATTACHMENTHARMONIC") {
        forceFields.push_back(std::make_unique<MTOCAttachment>());
    }
    else if(mtoc == "") {}
    else {
        log::error("MTOC FF not recognized. Exiting.");
        throw std::runtime_error("MTOC FF not recognized.");
    }

    // AFM specific.
    if(afm == "ATTACHMENTHARMONIC") {
        forceFields.push_back(std::make_unique<AFMAttachment>());
    }
    else if(afm == "") {}
    else {
        log::error("AFM FF not recognized. Exiting.");
        throw std::runtime_error("AFM FF not recognized.");
    }

    return forceFields;
}

} // namespace medyan

#endif

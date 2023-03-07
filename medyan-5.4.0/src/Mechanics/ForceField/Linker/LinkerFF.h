
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

#ifndef MEDYAN_LinkerFF_h
#define MEDYAN_LinkerFF_h

#include <vector>

#include "common.h"
#include "Mechanics/ForceField/Linker/LinkerStretching.h"
#include "Mechanics/ForceField/Linker/LinkerStretchingHarmonic.h"
#include "ForceField.h"

namespace medyan {

inline auto createLinkerForceFields(std::string_view stretching, std::string_view bending, std::string_view twisting) {
    std::vector<std::unique_ptr<ForceField>> forceFields;

    if (stretching == "HARMONIC")
        forceFields.push_back(
            std::make_unique<LinkerStretching<LinkerStretchingHarmonic>>()
        );
    else if(stretching == "") {}
    else {
        log::error("Linker stretching FF {} not recognized. Exiting.", stretching);
        throw std::runtime_error("Linker stretching FF not recognized.");
    }

    return forceFields;
}

} // namespace medyan

#endif

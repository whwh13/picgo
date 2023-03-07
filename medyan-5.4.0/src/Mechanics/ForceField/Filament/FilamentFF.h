
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

#ifndef MEDYAN_FilamentFF_h
#define MEDYAN_FilamentFF_h

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "common.h"
#include "Mechanics/ForceField/Filament/FilamentBending.h"
#include "Mechanics/ForceField/Filament/FilamentBendingCosine.h"
#include "Mechanics/ForceField/Filament/FilamentBendingHarmonic.h"
#include "Mechanics/ForceField/Filament/FilamentStretching.h"
#include "Mechanics/ForceField/Filament/FilamentStretchingHarmonic.h"
#include "Mechanics/ForceField/Filament/FilamentStretchingandBending.h"
#include "Mechanics/ForceField/Filament/FilamentStretchingHarmonicandBendingCosine.h"
#include "Mechanics/ForceField/Filament/FilamentStretchingHarmonicandBendingHarmonic.h"

namespace medyan {

inline auto createFilamentForceFields(std::string_view stretching, std::string_view bending, std::string_view twisting) {
    std::vector<std::unique_ptr<ForceField>> res;

    if(stretching == "HARMONIC") {
        res.push_back(std::make_unique<FilamentStretching<FilamentStretchingHarmonic>>());
    }
    else if(stretching == "") {}
    else {
        log::error("Filament stretching FF not recognized. Exiting.");
        throw std::runtime_error("Filament stretching FF not recognized.");
    }

    if(bending == "HARMONIC") {
        res.push_back(std::make_unique<FilamentBending<FilamentBendingHarmonic>>());
    }
    else if(bending == "COSINE") {
        res.push_back(std::make_unique<FilamentBending<FilamentBendingCosine>>());
    }
    else if(bending == "") {}
    else {
        log::error("Filament bending FF not recognized. Exiting.");
        throw std::runtime_error("Filament bending FF not recognized.");
    }

    return res;
}

} // namespace medyan

#endif

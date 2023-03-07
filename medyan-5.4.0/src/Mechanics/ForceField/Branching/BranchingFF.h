
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

#ifndef MEDYAN_BranchingFF_h
#define MEDYAN_BranchingFF_h

#include <vector>

#include "common.h"
#include "Mechanics/ForceField/Branching/BranchingStretching.h"
#include "Mechanics/ForceField/Branching/BranchingStretchingHarmonic.h"
#include "Mechanics/ForceField/Branching/BranchingBending.h"
#include "Mechanics/ForceField/Branching/BranchingBendingCosine.h"
#include "Mechanics/ForceField/Branching/BranchingDihedral.h"
#include "Mechanics/ForceField/Branching/BranchingDihedralCosine.h"
#include "Mechanics/ForceField/Branching/BranchingDihedralCosineV2.h"
#include "Mechanics/ForceField/Branching/BranchingDihedralQuadratic.hpp"
#include "Mechanics/ForceField/Branching/BranchingDihedralQuadraticV2.h"
#include "Mechanics/ForceField/Branching/BranchingPosition.h"
#include "Mechanics/ForceField/Branching/BranchingPositionCosine.h"

namespace medyan {

inline auto createBranchingForceFields(
    std::string_view stretching,
    std::string_view bending,
    std::string_view dihedral,
    std::string_view position
) {
    std::vector<std::unique_ptr<ForceField>> forceFields;

    if(stretching == "HARMONIC") {
        forceFields.push_back(
            std::make_unique<BranchingStretching<BranchingStretchingHarmonic>>()
        );
    }
    else if(stretching == "") {}
    else {
        log::error("Branching stretching FF {} not recognized. Exiting.", stretching);
        throw std::runtime_error("Branching stretching FF not recognized.");
    }

    if(bending == "COSINE") {
        forceFields.push_back(
            std::make_unique<BranchingBending<BranchingBendingCosine>>()
        );
    }
    else if(bending == "") {}
    else {
        log::error("Branching bending FF {} not recognized. Exiting.", bending);
        throw std::runtime_error("Branching bending FF not recognized.");
    }

    if(dihedral == "COSINE") {
        forceFields.push_back(
            std::make_unique<BranchingDihedral<BranchingDihedralCosine>>()
        );
    }
    else if(dihedral == "COSINEV2") {
        forceFields.push_back(
            std::make_unique<BranchingDihedral<BranchingDihedralCosineV2>>()
        );
    }
    else if(dihedral == "QUADRATIC") {
        forceFields.push_back(
            std::make_unique<BranchingDihedral< BranchingDihedralQuadratic>>()
        );
    }
    else if(dihedral == "QUADRATICV2") {
        forceFields.push_back(
            std::make_unique<BranchingDihedral< BranchingDihedralQuadraticV2>>()
        );
    }
    else if(dihedral == "") {}
    else {
        log::error("Branching dihedral FF {} not recognized.", dihedral);
        throw std::runtime_error("Unrecognized branching dihedral force field");
    }

    if(position == "COSINE") {
        forceFields.push_back(
            std::make_unique<BranchingPosition<BranchingPositionCosine>>()
        );
    }
    else if(position == "") {}
    else {
        log::error("Branching position FF {} not recognized. Exiting.", position);
        throw std::runtime_error("Branching position FF not recognized.");
    }

    return forceFields;
}

} // namespace medyan

#endif

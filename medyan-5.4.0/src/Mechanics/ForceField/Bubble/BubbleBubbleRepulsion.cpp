
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

#include "Mechanics/ForceField/Bubble/BubbleBubbleRepulsion.h"

#include "Bubble.h"
#include "Bead.h"
#include "Structure/DofSerializer.hpp"

namespace medyan {
void BubbleBubbleRepulsion::vectorize(const FFCoordinateStartingIndex& si) {
    ps_ = si.ps;

    auto& nl = ps_->opBubbleBubbleNL.value();

    pairInteractions_.clear();
    for(auto& bb : ps_->bubbles) {
        for(auto& bb2Index : nl.getNeighbors(bb.sysIndex)) {
            auto& bb2 = ps_->bubbles[bb2Index];
            auto& pair = pairInteractions_.emplace_back();
            pair.coordIndex1 = findBubbleCoordIndex(bb, si);
            pair.coordIndex2 = findBubbleCoordIndex(bb2, si);
            pair.krep = bb.getRepulsionConst();
            pair.slen = bb.getScreeningLength();
            pair.radius1 = bb.getRadius();
            pair.radius2 = bb2.getRadius();
        }
    }
}


floatingpoint BubbleBubbleRepulsion::computeEnergy(floatingpoint* coord) {
    
    floatingpoint energy = 0;
    
    for (auto& pair : pairInteractions_) {

        const auto u = _FFType.energy(
            coord,
            pair.coordIndex1, pair.coordIndex2,
            pair.radius1, pair.radius2,
            pair.krep, pair.slen
        );

        energy += u;
    }
    
    return energy;
}

void BubbleBubbleRepulsion::computeForces(floatingpoint *coord, floatingpoint *force) {
    
    for (auto& pair : pairInteractions_) {
        _FFType.forces(
            coord, force,
            pair.coordIndex1, pair.coordIndex2,
            pair.radius1, pair.radius2,
            pair.krep, pair.slen
        );
    }
}


} // namespace medyan


//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "MTOCBending.h"

#include "MTOCBendingCosine.h"

#include "MTOC.h"
#include "Bubble.h"
#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"
#include "Structure/DofSerializer.hpp"

namespace medyan {
void MTOCBending::vectorize(const FFCoordinateStartingIndex& si) {

    ps_ = si.ps;

    interactions_.clear();
    for(auto& mtoc : ps_->mtocs) {
        auto& bb = mtoc.getBubble(*ps_);
        for(auto pf : mtoc.getFilaments()) {
            interactions_.push_back({
                findBubbleCoordIndex(bb, si),
                findBeadCoordIndex(*pf->getMinusEndCylinder()->getFirstBead(), si),
                findBeadCoordIndex(*pf->getMinusEndCylinder()->getSecondBead(), si),
                bb.getMTOCBendingK(),
            });
        }
    }
}

floatingpoint MTOCBending::computeEnergy(floatingpoint* coord) {

    floatingpoint energy = 0;
    for(auto& eachInteraction : interactions_) {
        const auto u = _FFType.energy(
            coord,
            eachInteraction.bubbleCoordIndex,
            eachInteraction.beadCoordIndex1,
            eachInteraction.beadCoordIndex2,
            eachInteraction.kbend
        );
        energy += u;
    }
    
    return energy;
}

void MTOCBending::computeForces(floatingpoint *coord, floatingpoint *force) {
    for(auto& eachInteraction : interactions_) {
        _FFType.forces(
            coord, force,
            eachInteraction.bubbleCoordIndex,
            eachInteraction.beadCoordIndex1,
            eachInteraction.beadCoordIndex2,
            eachInteraction.kbend
        );
    }
    
}


} // namespace medyan

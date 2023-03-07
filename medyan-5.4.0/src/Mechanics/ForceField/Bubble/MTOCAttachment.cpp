
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

#include "MTOCAttachment.h"

#include "MTOCAttachmentHarmonic.h"

#include "MTOC.h"
#include "Bubble.h"
#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"
#include "Structure/DofSerializer.hpp"

namespace medyan {

void MTOCAttachment::vectorize(const FFCoordinateStartingIndex& si) {

    ps_ = si.ps;

    pairInteractions_.clear();
    for(auto& mtoc : ps_->mtocs) {
        auto& bb = mtoc.getBubble(*ps_);
        for(auto pf : mtoc.getFilaments()) {
            pairInteractions_.push_back({
                findBubbleCoordIndex(bb, si),
                findBeadCoordIndex(*pf->getMinusEndCylinder()->getFirstBead(), si),
                mtoc.attachmentStretchingK,
                bb.getRadius(),
            });
        }
    }
}


floatingpoint MTOCAttachment::computeEnergy(floatingpoint* coord) {

    floatingpoint energy = 0;
    for(auto& pair : pairInteractions_) {
        const auto u = impl.energy(
            coord,
            pair.bubbleCoordIndex, pair.beadCoordIndex,
            pair.kstr, pair.radius
        );
        energy += u;
    }
    return energy;
}

void MTOCAttachment::computeForces(floatingpoint *coord, floatingpoint *force) {
    for(auto& pair : pairInteractions_) {
        impl.forces(
            coord, force,
            pair.bubbleCoordIndex, pair.beadCoordIndex,
            pair.kstr, pair.radius
        );
    }
}


} // namespace medyan

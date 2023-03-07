
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

#include "AFMAttachment.h"

#include "AFMAttachmentHarmonic.h"

#include "AFM.h"
#include "Bubble.h"
#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"
#include "Structure/DofSerializer.hpp"

namespace medyan {
void AFMAttachment::vectorize(const FFCoordinateStartingIndex& si) {

    ps_ = si.ps;

    pairInteractions_.clear();
    for(auto& afm : ps_->afms) {
        auto& bb = afm.getBubble(*ps_);
        for(auto pf : afm.getFilaments()) {
            pairInteractions_.push_back({
                findBubbleCoordIndex(bb, si),
                findBeadCoordIndex(*pf->getMinusEndCylinder()->getFirstBead(), si),
                afm.attachmentStretchingK,
                bb.getRadius(),
            });
        }
    }
}

floatingpoint AFMAttachment::computeEnergy(floatingpoint *coord) {
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

void AFMAttachment::computeForces(floatingpoint *coord, floatingpoint *force) {
    for(auto& pair : pairInteractions_) {
        impl.forces(
            coord, force,
            pair.bubbleCoordIndex, pair.beadCoordIndex,
            pair.kstr, pair.radius
        );
    }
}

} // namespace medyan

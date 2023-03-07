
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

#include "BoundaryCylinderAttachment.h"

#include "BoundaryCylinderAttachmentHarmonic.h"
#include "BoundaryElement.h"

#include "Bubble.h"
#include "Bead.h"

namespace medyan {
template <class BAttachmentInteractionType>
void BoundaryCylinderAttachment<BAttachmentInteractionType>::vectorize(const FFCoordinateStartingIndex& si, const SimulConfig& conf) {
    
    //first coord in beadset is bead, then pin position
    beadSet.assign(Bead::getPinnedBeads().size(), 0);
    kattr.assign(Bead::getPinnedBeads().size(), 0);
    pins.resize(Bead::getPinnedBeads().size());
    
    int i = 0;
    for(auto b : Bead::getPinnedBeads()) {

        beadSet[n * i] = b->getIndex() * 3 + si.bead;
        kattr[n * i] = SysParams::Mechanics().pinK;

        pins[n * i] = mathfunc::vector2Vec< 3, floatingpoint >(b->getPinPosition());

        i++;
    }
}


template <class BAttachmentInteractionType>
floatingpoint BoundaryCylinderAttachment<BAttachmentInteractionType>::computeEnergy(floatingpoint *coord) {

    return _FFType.energy(coord, beadSet.data(), kattr.data(), pins);

}

template <class BAttachmentInteractionType>
void BoundaryCylinderAttachment<BAttachmentInteractionType>::computeForces(floatingpoint *coord, floatingpoint *f) {

    _FFType.forces(coord, f, beadSet.data(), kattr.data(), pins);
}

///Template specializations
template floatingpoint BoundaryCylinderAttachment<BoundaryCylinderAttachmentHarmonic>::computeEnergy(floatingpoint *coord);
template void BoundaryCylinderAttachment<BoundaryCylinderAttachmentHarmonic>::computeForces(floatingpoint *coord, floatingpoint *f);
template void BoundaryCylinderAttachment<BoundaryCylinderAttachmentHarmonic>::vectorize(const FFCoordinateStartingIndex&, const SimulConfig&);

} // namespace medyan

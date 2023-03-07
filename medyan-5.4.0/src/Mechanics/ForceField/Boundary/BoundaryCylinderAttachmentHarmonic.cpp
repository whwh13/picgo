
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

#include "BoundaryCylinderAttachmentHarmonic.h"
#include "BoundaryCylinderAttachment.h"

#include "Bead.h"

#include "MathFunctions.h"

namespace medyan {
using namespace mathfunc;

floatingpoint BoundaryCylinderAttachmentHarmonic::energy(
    floatingpoint *coord, int *beadSet,
    floatingpoint *kattr, const std::vector< Vec< 3, floatingpoint > >& pins
) const {
    

    int n = BoundaryCylinderAttachment<BoundaryCylinderAttachmentHarmonic>::n;
    int nint = Bead::getPinnedBeads().size();

    floatingpoint U_i;
    floatingpoint U = 0;

    for(int i = 0; i < nint; i += 1) {

        const auto coord1 = makeRefVec< 3 >(coord + beadSet[n * i]);

        const auto distsq = distance2(coord1, pins[i]);
        U_i = 0.5 * kattr[i] * distsq;

        U += U_i;
    }
    return U;
}


void BoundaryCylinderAttachmentHarmonic::forces(
    floatingpoint *coord, floatingpoint *f, int *beadSet,
    floatingpoint *kattr, const std::vector< Vec< 3, floatingpoint > >& pins
) const {
    
    int n = BoundaryCylinderAttachment<BoundaryCylinderAttachmentHarmonic>::n;
    int nint = Bead::getPinnedBeads().size();

    for(int i = 0; i < nint; i += 1) {

        const auto coord1 = makeRefVec< 3 >(coord + beadSet[n * i]);
        auto       force1 = makeRefVec< 3 >(f     + beadSet[n * i]);

        force1 += kattr[i] * (pins[i] - coord1);
    }
}

} // namespace medyan


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

#ifndef MEDYAN_MTOCBendingCosine_h
#define MEDYAN_MTOCBendingCosine_h

#include "common.h"
#include "Util/Math/Vec.hpp"

namespace medyan {

/// A harmonic potential used by the MTOCAttachment template.
struct MTOCBendingCosine {
    floatingpoint energy(
        const floatingpoint *coord,
        Index bubbleCoordIndex, Index beadCoordIndex1, Index beadCoordIndex2,
        floatingpoint kbend
    ) const {
        auto coord1 = makeRefVec<3>(coord + bubbleCoordIndex);
        auto coord2 = makeRefVec<3>(coord + beadCoordIndex1);
        auto coord3 = makeRefVec<3>(coord + beadCoordIndex2);
        auto r12 = coord2 - coord1;
        auto r23 = coord3 - coord2;
        auto l1 = magnitude(r12);
        auto l2 = magnitude(r23);
        auto cosTheta = dot(r12, r23) / (l1 * l2);

        return kbend * (1 - cosTheta);
    }


    void forces(
        const floatingpoint *coord, floatingpoint *force,
        Index bubbleCoordIndex, Index beadCoordIndex1, Index beadCoordIndex2,
        floatingpoint kbend
    ) const {
        auto coord1 = makeRefVec<3>(coord + bubbleCoordIndex);
        auto coord2 = makeRefVec<3>(coord + beadCoordIndex1);
        auto coord3 = makeRefVec<3>(coord + beadCoordIndex2);
        auto force1 = makeRefVec<3>(force + bubbleCoordIndex);
        auto force2 = makeRefVec<3>(force + beadCoordIndex1);
        auto force3 = makeRefVec<3>(force + beadCoordIndex2);

        auto r12 = coord2 - coord1;
        auto r23 = coord3 - coord2;
        auto l1 = magnitude(r12);
        auto l2 = magnitude(r23);
        auto dotr = dot(r12, r23);
        auto invL1 = 1 / l1;
        auto invL2 = 1 / l2;

        auto A = invL1 * invL2;
        auto B = dotr * invL1 * A * A * l2;
        auto C = dotr * invL2 * A * A * l1;

        force1 -= kbend * (A * r23 - B * r12);
        force2 -= kbend * (A * (r12 - r23) + B * r12 - C * r23);
        force3 -= kbend * (C * r23 - A * r12);
    }
    
};

} // namespace medyan

#endif



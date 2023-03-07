
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

#ifndef MEDYAN_AFMAttachmentHarmonic_h
#define MEDYAN_AFMAttachmentHarmonic_h

#include "common.h"
#include "Util/Math/Vec.hpp"

namespace medyan {

/// A harmonic potential used by the AFMAttachment template.
struct AFMAttachmentHarmonic {
    floatingpoint energy(
        const floatingpoint *coord,
        Index bubbleCoordIndex, Index beadCoordIndex,
        floatingpoint kstr, floatingpoint radius
    ) const {
        auto coord1 = makeRefVec<3>(coord + bubbleCoordIndex);
        auto coord2 = makeRefVec<3>(coord + beadCoordIndex);
        auto diff = distance(coord1, coord2) - radius;
        
        auto energy = (kstr / 2) * diff * diff;
        
        return energy;
    }

    void forces(
        const floatingpoint *coord, floatingpoint *force,
        Index bubbleCoordIndex, Index beadCoordIndex,
        floatingpoint kstr, floatingpoint radius
    ) const {
        auto coord1 = makeRefVec<3>(coord + bubbleCoordIndex);
        auto coord2 = makeRefVec<3>(coord + beadCoordIndex);
        auto force1 = makeRefVec<3>(force + bubbleCoordIndex);
        auto force2 = makeRefVec<3>(force + beadCoordIndex);

        auto dist = distance(coord1, coord2);
        auto f0 = kstr * (dist - radius) / dist;

        force1 += f0 * (coord2 - coord1);
        force2 += f0 * (coord1 - coord2);
    }
};

} // namespace medyan

#endif


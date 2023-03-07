
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

#ifndef MEDYAN_BubbleCylinderRepulsionExp_h
#define MEDYAN_BubbleCylinderRepulsionExp_h

#include <vector>

#include "common.h"
#include "MathFunctions.h"

namespace medyan {

/// A exponential repulsive potential used by the BubbleCylinderRepulsion template.
struct BubbleBeadRepulsionExp {
    floatingpoint energy(
        const floatingpoint* coord,
        Index bubbleCoordIndex, Index beadCoordIndex,
        floatingpoint krep, floatingpoint slen, floatingpoint radius
    ) const {
        auto bubbleCoord = makeRefVec<3>(coord + bubbleCoordIndex);
        auto beadCoord = makeRefVec<3>(coord + beadCoordIndex);

        auto dist = distance(bubbleCoord, beadCoord);
        auto effd = dist - radius;

        auto exponent = -effd / slen;
        return krep * std::exp(exponent);
    }

    void forces(
        const floatingpoint* coord, floatingpoint* force,
        Index bubbleCoordIndex, Index beadCoordIndex,
        floatingpoint krep, floatingpoint slen, floatingpoint radius
    ) const {
        auto bubbleCoord = makeRefVec<3>(coord + bubbleCoordIndex);
        auto beadCoord = makeRefVec<3>(coord + beadCoordIndex);
        auto bubbleForce = makeRefVec<3>(force + bubbleCoordIndex);
        auto beadForce = makeRefVec<3>(force + beadCoordIndex);

        auto dist = distance(bubbleCoord, beadCoord);
        auto effd = dist - radius;

        auto exponent = -effd / slen;
        auto forceFactor = krep * std::exp(exponent) / (dist * slen);

        bubbleForce += forceFactor * (bubbleCoord - beadCoord);
        beadForce += forceFactor * (beadCoord - bubbleCoord);
    }

    template< typename VT1, typename VT2 >
	floatingpoint loadForces(
        const VT1&    coord1,
        const VT2&    coord2,
        floatingpoint radius,
        floatingpoint kRep,
        floatingpoint screenLength
    ) const {

        const auto dist = medyan::distance(coord1, coord2);

        floatingpoint effd = dist - radius;

        floatingpoint R = -effd / screenLength;
        return kRep * exp(R) / screenLength;

    }
};

} // namespace medyan

#endif

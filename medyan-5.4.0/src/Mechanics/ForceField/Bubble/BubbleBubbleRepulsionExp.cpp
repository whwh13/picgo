
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

#include "BubbleBubbleRepulsionExp.h"

#include "MathFunctions.h"

namespace medyan {
using namespace mathfunc;

floatingpoint BubbleBubbleRepulsionExp::energy(
    const floatingpoint* coord,
    Index i1, Index i2, floatingpoint r1, floatingpoint r2,
    floatingpoint kRep, floatingpoint screenLength
) const {
    
    floatingpoint dist = distance(
        makeRefVec< 3 >(coord + i1),
        makeRefVec< 3 >(coord + i2)
    );
    
    floatingpoint effd = dist - r1 - r2;
    
    floatingpoint R = -effd / screenLength;
    return kRep * exp(R);
}

void BubbleBubbleRepulsionExp::forces(
    const floatingpoint* coord, floatingpoint* force,
    Index i1, Index i2,
    floatingpoint r1, floatingpoint r2,
    floatingpoint kRep, floatingpoint screenLength
) const {

    const auto coord1 = makeRefVec< 3 >(coord + i1);
    const auto coord2 = makeRefVec< 3 >(coord + i2);
    auto force1 = makeRefVec< 3 >(force + i1);
    auto force2 = makeRefVec< 3 >(force + i2);

    //get dist
    floatingpoint dist = distance(coord1, coord2);

    floatingpoint effd = dist - r1 - r2;
    
    floatingpoint R = -effd / screenLength;
    floatingpoint f0 = kRep * exp(R) / screenLength;

    //get norm
    auto norm = normalizedVector(coord2 - coord1);

    force1 -= f0 * norm;
    force2 += f0 * norm;
}

} // namespace medyan

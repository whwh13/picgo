
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

#ifndef MEDYAN_BubbleBubbleRepulsionExp_h
#define MEDYAN_BubbleBubbleRepulsionExp_h

#include <vector>
#include <cmath>

#include "common.h"

namespace medyan {

/// A exponential repulsive potential used by the BubbleBubbleRepulsion template.
struct BubbleBubbleRepulsionExp {
    
    floatingpoint energy(
        const floatingpoint* coord,
        Index i1, Index i2,
        floatingpoint, floatingpoint, floatingpoint, floatingpoint
    ) const;
    
    void forces(
        const floatingpoint* coord, floatingpoint* force,
        Index i1, Index i2,
        floatingpoint, floatingpoint, floatingpoint, floatingpoint
    ) const;
};

} // namespace medyan

#endif

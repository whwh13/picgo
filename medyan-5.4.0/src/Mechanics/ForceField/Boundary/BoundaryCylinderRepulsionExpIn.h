
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

#ifndef MEDYAN_BoundaryCylinderRepulsionExpIn_h
#define MEDYAN_BoundaryCylinderRepulsionExpIn_h

#include <vector>
#include <cmath>

#include "common.h"

namespace medyan {
//FORWARD DECLARATIONS
class Bead;

/// A exponential repulsive potential used by the BoundaryCylinderRepulsion template.
class BoundaryCylinderRepulsionExpIn {
    
public:
    floatingpoint energy(floatingpoint *coord, int *beadSet,
                         floatingpoint *krep, floatingpoint *slen, int *nneighbors);
    
    void forces(floatingpoint *coord, floatingpoint *f, int *beadSet,
                floatingpoint *krep, floatingpoint *slen, int *nneighbors);
    
    floatingpoint loadForces(floatingpoint r, floatingpoint krep , floatingpoint slen) const;
    
};

} // namespace medyan

#endif

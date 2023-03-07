
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

#ifndef MEDYAN_BoundaryCylinderAttachmentHarmonic_h
#define MEDYAN_BoundaryCylinderAttachmentHarmonic_h

#include <vector>
#include <cmath>

#include "common.h"
#include "Util/Math/Vec.hpp"

namespace medyan {
//FORWARD DECLARATIONS
class Bead;

/// A harmonic attractive potential used by the BoundaryCylinderAttachment template.
class BoundaryCylinderAttachmentHarmonic {
    
public:
    floatingpoint energy(
        floatingpoint *coord, int *beadSet,
        floatingpoint *kattr, const std::vector< medyan::Vec< 3, floatingpoint > >& pins
    ) const;
    
    void forces(
        floatingpoint *coord, floatingpoint *f, int *beadSet,
        floatingpoint *kattr, const std::vector< medyan::Vec< 3, floatingpoint > >& pins
    ) const;
};

} // namespace medyan

#endif

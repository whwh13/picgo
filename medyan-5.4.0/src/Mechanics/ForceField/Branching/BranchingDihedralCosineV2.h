
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

#ifndef MEDYAN_BranchingDihedralCosineV2_h
#define MEDYAN_BranchingDihedralCosineV2_h

#include "common.h"

namespace medyan {
//FORWARD DECLARATIONS
class Bead;

/// A cosine potential used by the BranchingDihedralTemplate.
class BranchingDihedralCosineV2 {

public:
	floatingpoint energy(floatingpoint *coord, size_t nint,
	                     unsigned int *beadSet, floatingpoint *kdih, floatingpoint *pos);

	void forces(floatingpoint *coord, floatingpoint *f, size_t nint,
	            unsigned int *beadSet, floatingpoint *kdih, floatingpoint *pos,
	            floatingpoint *stretchforce);
};

} // namespace medyan

#endif

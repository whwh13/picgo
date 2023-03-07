
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

#ifndef MEDYAN_BranchingBendingCosine_h
#define MEDYAN_BranchingBendingCosine_h

#include "common.h"

namespace medyan {
//FORWARD DECLARATIONS
class Bead;

/// A cosine potential used by the BranchingBending template.
class BranchingBendingCosine {
    
public:
    floatingpoint energy(floatingpoint *coord, int *beadSet,
                  floatingpoint *kbend, floatingpoint *eqt);
    
    void forces(floatingpoint *coord, floatingpoint *f, int *beadSet,
                floatingpoint *kbend, floatingpoint *eqt, floatingpoint *stretchforce);
#ifdef CUDAACCL
    void optimalblocksnthreads(int nint);

    floatingpoint* energy(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *kbend, floatingpoint *eqt, int *params);

    floatingpoint* energy(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *kbend, floatingpoint *eqt, floatingpoint *z, int *params);

    void forces(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *kbend, floatingpoint *eqt, int *params);
    void deallocate();
    static void checkforculprit();
    floatingpoint *gU_i;
    floatingpoint *gU_sum;
    char *gFF, *ginteraction;
    vector<int> blocksnthreadse;
    vector<int> blocksnthreadsez;
    vector<int> blocksnthreadsf;
    vector<int> bntaddvec2;
    cudaStream_t stream = NULL;
#endif
};

} // namespace medyan

#endif


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

#ifndef MEDYAN_BranchingPositionCosine_h
#define MEDYAN_BranchingPositionCosine_h

#include "common.h"

namespace medyan {
//FORWARD DECLARATIONS
class Bead;

/// A cosine potential used by the BranchingPosition template.
class BranchingPositionCosine {
    
public:
    floatingpoint energy(const floatingpoint *coord, unsigned int *beadSet,
                  floatingpoint *kpos, floatingpoint *pos) const;
    
    [[deprecated]] floatingpoint energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                  floatingpoint *kpos, floatingpoint *pos, floatingpoint d);
    
    void forces(const floatingpoint *coord, floatingpoint *f, unsigned int *beadSet,
                floatingpoint *kpos, floatingpoint *pos, floatingpoint *stretchforce) const;
#ifdef CUDAACCL
    void optimalblocksnthreads(int nint);
    floatingpoint* energy(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *kpos, floatingpoint *pos, int *params);

    floatingpoint* energy(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *kpos, floatingpoint *pos, floatingpoint *z, int *params);

    void forces(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *kpos, floatingpoint *pos, int *params);
    void deallocate();
    vector<int> blocksnthreadse;
    vector<int> blocksnthreadsez;
    vector<int> blocksnthreadsf;
    vector<int> bntaddvec2;
    static void checkforculprit();
    floatingpoint *gU_i;
    floatingpoint *gU_sum;
    char *gFF, *ginteraction;
    cudaStream_t stream = NULL;
#endif
};

} // namespace medyan

#endif

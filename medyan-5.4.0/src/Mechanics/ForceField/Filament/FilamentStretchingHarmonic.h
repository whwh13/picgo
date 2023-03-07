
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

#ifndef MEDYAN_FilamentStretchingHarmonic_h
#define MEDYAN_FilamentStretchingHarmonic_h

#include "common.h"

namespace medyan {
//FORWARD DECLARATIONS
class Bead;

/// A harmonic potential used by the FilamentStretching and MTOCAttachment template.
class FilamentStretchingHarmonic {
    
public:
    floatingpoint energy(floatingpoint *coord, int *beadSet,
                  floatingpoint *kstr, floatingpoint *eql);

    void forces(floatingpoint *coord, floatingpoint *f, int *beadSet,
                floatingpoint *kstr, floatingpoint *eql);
#ifdef CUDAACCL
    void optimalblocksnthreads(int nint, cudaStream_t stream);

    floatingpoint* energy(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *kstr, floatingpoint *eql, int *params);

    floatingpoint* energy(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *kstr, floatingpoint *eql, floatingpoint *z, int *params);

    void forces(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *kstr, floatingpoint *eql, int *params);
    void deallocate();
    static void checkforculprit();
    vector<int> blocksnthreadse;
    vector<int> blocksnthreadsez;
    vector<int> blocksnthreadsf;
    vector<int> bntaddvec2;
    floatingpoint *gU_i;
    floatingpoint *gU_sum;
    char *gFF, *ginteraction;
    cudaStream_t stream = NULL;
#endif
#ifdef CROSSCHECK
    floatingpoint energy(Bead*, Bead*, floatingpoint, floatingpoint);
    floatingpoint energy(Bead*, Bead*, floatingpoint, floatingpoint, floatingpoint);
    
    void forces(Bead*, Bead*, floatingpoint, floatingpoint);
    void forcesAux(Bead*, Bead*, floatingpoint, floatingpoint);
#endif
};

} // namespace medyan

#endif

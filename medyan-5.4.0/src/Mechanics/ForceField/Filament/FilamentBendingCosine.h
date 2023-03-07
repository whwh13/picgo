
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

#ifndef MEDYAN_FilamentBendingCosine_h
#define MEDYAN_FilamentBendingCosine_h

#include "common.h"

namespace medyan {
//FORWARD DECLARATIONS
class Bead;

/// A cosine potential used by the [FilamentBending](@ref FilamentBending) template.
class FilamentBendingCosine {

public:

    floatingpoint energy(floatingpoint *coord, std::size_t nint, int *beadSet,
                  floatingpoint *kbend, floatingpoint *eqt);
    
    void forces(floatingpoint *coord, floatingpoint *f, std::size_t nint, int *beadSet,
                floatingpoint *kbend, floatingpoint *eqt);
#ifdef CUDAACCL
    void optimalblocksnthreads(int nint, cudaStream_t stream);

    floatingpoint* energy(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *kbend, floatingpoint *eqt, int *params);

    floatingpoint* energy(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *kbend, floatingpoint *eqt, floatingpoint *z, int *params);

    void forces(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *kbend, floatingpoint *eqt, int *params);
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
#ifdef CROSSCHECK
    floatingpoint energy(Bead*, Bead*, Bead*, floatingpoint, floatingpoint);
    floatingpoint energy(Bead*, Bead*, Bead*, floatingpoint, floatingpoint, floatingpoint);
    
    void forces(Bead*, Bead*, Bead*, floatingpoint, floatingpoint);
    void forcesAux(Bead*, Bead*, Bead*, floatingpoint, floatingpoint);

#endif
};

} // namespace medyan

#endif

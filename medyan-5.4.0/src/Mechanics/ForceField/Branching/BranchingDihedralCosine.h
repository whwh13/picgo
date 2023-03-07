
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

#ifndef MEDYAN_BranchingDihedralCosine_h
#define MEDYAN_BranchingDihedralCosine_h

#include "common.h"

namespace medyan {
//FORWARD DECLARATIONS
class Bead;

/// A cosine potential used by the BranchingDihedralTemplate.
class BranchingDihedralCosine {

public:
    floatingpoint energy(
        const floatingpoint *coord, size_t nint,
        unsigned int *beadSet, floatingpoint *kdih, floatingpoint *pos
    ) const;

    void forces(
        const floatingpoint *coord, floatingpoint *f, size_t nint,
        unsigned int *beadSet, floatingpoint *kdih, floatingpoint *pos,
        floatingpoint *stretchforce
    ) const;

	void forcesNumericalinteraction(floatingpoint *coord, floatingpoint *f, size_t nint,
			unsigned int *beadSet, floatingpoint *kdih, floatingpoint *pos, int
			interactionID, double *Nforce);

	template <class dataType = double>
	dataType energyininteractionperturbed(floatingpoint *coord, size_t nint,
	                                           unsigned int *beadSet, floatingpoint *kdih,
	                                           floatingpoint *pos, int interactionID,
	                                           const int perturbcoord, const int
	                                           perturbaxis, double delta);

	void testdihedral();
#ifdef CUDAACCL
    void optimalblocksnthreads(int nint);

    floatingpoint* energy(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *kdih,
                   floatingpoint *pos, int *params);

    floatingpoint* energy(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *kdih, floatingpoint *pos,
                   floatingpoint *z, int *params);

    void forces(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *kdih, floatingpoint *pos, int *params);
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
/*private:
  int counterE = 0;
  int counterF = 0;
  double EmagerrorMvsN = 0.0;
  double FmagerrorMvsN[4] = {0.0, 0.0, 0.0, 0.0};
  double FdirerrorMvsN[4] = {0.0, 0.0, 0.0, 0.0};
  //
  double EmagerrorLvsN  = 0.0;
  double FmagerrorLvsN[4] = {0.0, 0.0, 0.0, 0.0};
  double FdirerrorLvsN[4] = {0.0, 0.0, 0.0, 0.0};*/
};

} // namespace medyan

#endif


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

#ifndef MEDYAN_FilamentStretchingHarmonicandBendingCosine_h
#define MEDYAN_FilamentStretchingHarmonicandBendingCosine_h

#include "common.h"

namespace medyan {
//FORWARD DECLARATIONS
class Bead;

/// A cosine potential used by the [FilamentBending](@ref FilamentBending) template.
class FilamentStretchingHarmonicandBendingCosine {

public:
	void energy(floatingpoint *coord, std::size_t nint, int *beadSet,
	            floatingpoint *kstr, floatingpoint *kbend, floatingpoint *eql,
	            floatingpoint *eqt, floatingpoint* totalenergy, const int startID, const
	            int endID, int threadID);
	void energy(floatingpoint *coord, int *beadSetsansbending, floatingpoint *kstrsansbending,
	            floatingpoint *eqlsansbending, floatingpoint* totalenergy, const int startID,
	            const int endID, int threadID);

	void forces(floatingpoint *coord,
	            floatingpoint *f, size_t nint, int *beadSet,
	            floatingpoint *kstr, floatingpoint *kbend,
	            floatingpoint *eql, floatingpoint *eqt);

	void forces(floatingpoint *coord,  floatingpoint *f,
				int *beadSetsansbending, floatingpoint *kstrsansbending,
				floatingpoint *eqlsansbending,
				const int startID, const int endID, int threadID);

    void energy(floatingpoint *coord, std::size_t nint, int * cylSet,
                floatingpoint *cyllengthset, floatingpoint *cylbenddotproduct,
                floatingpoint *kstr, floatingpoint *kbend, floatingpoint *eql,
                floatingpoint *eqt, floatingpoint* totalenergy,
                const int startID, const int endID, int threadID);

    void energy(floatingpoint *coord, int * cylSetcylsansbending,
                floatingpoint *cyllengthset, floatingpoint *kstrsansbending,
                floatingpoint *eqlsansbending, floatingpoint* totalenergy, const int startID,
                const int endID, int threadID);

    void forces(floatingpoint *coord, floatingpoint *f, size_t nint, int *beadSet,
                int *cylSet, floatingpoint *cyllengthset, floatingpoint *cylbenddotproduct,
                 floatingpoint *kstr, floatingpoint *kbend,
                floatingpoint *eql, floatingpoint *eqt);

    void forces(floatingpoint *coord,  floatingpoint *f,  int *beadSet,
                int * cylSetcylsansbending, floatingpoint *cyllengthset,
                int *beadSetsansbending, floatingpoint *kstrsansbending,
                floatingpoint *eqlsansbending,
                const int startID, const int endID, int threadID);
};

} // namespace medyan

#endif

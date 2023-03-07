
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

#ifndef MEDYAN_CylinderExclVolRepulsion_h
#define MEDYAN_CylinderExclVolRepulsion_h

#include "common.h"

namespace medyan {
//FORWARD DECLARATIONS
class Bead;

/// Represents a repulsive excluded volume potential used by the
/// CylinderExclVolume template.
class CylinderExclVolRepulsion {
    
public:
    

    floatingpoint energy(floatingpoint *coord, const int *beadSet, const floatingpoint *krep, const floatingpoint* eqLengths, int nint);

    void forces(floatingpoint *coord, floatingpoint *f, const int *beadSet, const floatingpoint *krep, const floatingpoint* eqLengths, int nint);
    
    //returns vec of cyl cyl interaction energies if they are above  cylthresh
    vector<tuple<floatingpoint, int, vector<tuple<floatingpoint*,floatingpoint*,floatingpoint*,floatingpoint*, floatingpoint>>>> getCylEnergies(){ return cylEnergies;};
    
    void clearCylEnergies(){
        cylEnergies.clear();
    }
    // c1 c2 c3 c4 energy
    vector<tuple<floatingpoint, int, vector<tuple<floatingpoint*,floatingpoint*,floatingpoint*,floatingpoint*, floatingpoint>>>> cylEnergies;
    
    vector<floatingpoint> uniqueTimes;

private:

	floatingpoint energyN(floatingpoint *coord, const int *beadSet,
	                      const floatingpoint *krep, const floatingpoint* eqLengths, int intID, bool movebeads = false);

	void forceN(floatingpoint *coord, floatingpoint *f, const int *beadSet,
	                      const floatingpoint *krep, const floatingpoint* eqLengths, int intID, bool movebeads = false);

	doubleprecision getenergyintegrand(doubleprecision& a, doubleprecision& b,
			doubleprecision& c, doubleprecision& d, doubleprecision& e, doubleprecision& F,
			doubleprecision s, doubleprecision t){
		doubleprecision r_2 = c + 2*s*e + s*s*a - 2*t*F - 2*s*t*d + t*t*b;
		return (1.0/(r_2 * r_2));
	}

	void getforceintegrand(doubleprecision& a, doubleprecision& b,
	                                   doubleprecision& c, doubleprecision& d, doubleprecision& e, doubleprecision& F,
	                                   doubleprecision s, doubleprecision t,
	                                   doubleprecision* integrand){
		doubleprecision r_2 = c + 2*s*e + s*s*a - 2*t*F - 2*s*t*d + t*t*b;
		doubleprecision x = (1.0/(r_2 * r_2 * r_2));
		integrand[0] = x;
		integrand[1] = s*x;
		integrand[2] = s*s*x;
		integrand[3] = s*t*x;
		integrand[4] = t*x;
		integrand[5] = t*t*x;
	}

#ifdef CUDAACCL
    void optimalblocksnthreads(int nint, cudaStream_t stream);
    void deallocate();
    floatingpoint* energy(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *krep, int *params);

    floatingpoint* energy(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *krep, floatingpoint *z, int *params);

    void forces(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *krep, int *params);
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
#ifdef CROSSCHECK
    floatingpoint energy(Bead*, Bead*, Bead*, Bead*, floatingpoint Krepuls);
    floatingpoint energy(Bead*, Bead*, Bead*, Bead*, floatingpoint Krepuls, floatingpoint d);
    
    void forces(Bead*, Bead*, Bead*, Bead*, floatingpoint Krepuls);
    void forcesAux(Bead*, Bead*, Bead*, Bead*, floatingpoint Krepuls);
#endif
};

} // namespace medyan

#endif


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

#include "BoundaryBubbleRepulsionExp.h"
#include "BoundaryBubbleRepulsion.h"
#include "BoundaryElement.h"

#include "Bead.h"

namespace medyan {
floatingpoint BoundaryBubbleRepulsionExp::energy(floatingpoint *coord, int *beadSet, floatingpoint *krep, floatingpoint *slen, int *nneighbors) {
    int nb, nc;
	floatingpoint *coord1, R, r, U_i;
	floatingpoint U = 0.0;
    int Cumnc=0;
    auto beList = BoundaryElement::getBoundaryElements();
    nb = beList.size();

    for (int ib = 0; ib < nb; ib++) {

        auto be = beList[ib];
        nc = nneighbors[ib];

        for (int ic = 0; ic < nc; ic++) {

            coord1 = &coord[beadSet[Cumnc + ic]];
            r = be->distance(coord1);

            //Boundary-Bubble repulsion only considers the distance between
            //boundary element and the bubble center.
            R = -r / slen[Cumnc + ic];
            U_i = krep[Cumnc + ic] * exp(R);

            U += U_i;
        }
        Cumnc += nc;
    }
    return U;

}

void BoundaryBubbleRepulsionExp::forces(floatingpoint *coord, floatingpoint *f,
		int *beadSet, floatingpoint *krep, floatingpoint *slen, int *nneighbors) {
    
    int nb, nc;
	floatingpoint *coord1, *force1, R, r, f0;
    auto beList = BoundaryElement::getBoundaryElements();
    nb = beList.size();
    int Cumnc=0;
    
    for (int ib = 0; ib < nb; ib++) {

        auto be = beList[ib];
        nc = nneighbors[ib];
        for(int ic = 0; ic < nc; ic++) {
            coord1 = &coord[beadSet[ Cumnc + ic]];
            force1 = &f[beadSet[ Cumnc + ic]];
            r = be->distance(coord1);
            auto norm = be->normal(coord1);

            R = -r / slen[Cumnc + ic];
            f0 = krep[Cumnc + ic] * exp(R)/ slen[Cumnc + ic];
            force1[0] += f0 *norm[0];
            force1[1] += f0 *norm[1];
            force1[2] += f0 *norm[2];

        }
        Cumnc+=nc;
    }
    
}

} // namespace medyan

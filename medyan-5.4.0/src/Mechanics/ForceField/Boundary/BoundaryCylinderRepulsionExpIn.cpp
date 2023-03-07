
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

#include "BoundaryCylinderRepulsionExpIn.h"
#include "BoundaryCylinderRepulsionIn.h"

#include "BoundaryElement.h"
#include "Bead.h"

#include "cross_check.h"
#include "Cylinder.h"
#include "MathFunctions.h"

namespace medyan {
using namespace mathfunc;

floatingpoint BoundaryCylinderRepulsionExpIn::energy(floatingpoint *coord, int *beadSet,
                                                   floatingpoint *krep, floatingpoint *slen, int *nneighbors) {
    
    int nb, nc;
    floatingpoint *coord1, R, r, U_i;
    floatingpoint U = 0.0;
    int Cumnc=0;
    auto beList = BoundaryElement::getBoundaryElements();
    nb = beList.size();
    
    //    for (int ib = 0; ib < nb; ib++) {
    //        auto be = beList[ib];
    //        be->printSelf();
    //    }
    
    for (int ib = 0; ib < nb; ib++) {
        
        auto be = beList[ib];
        nc = nneighbors[ib];
        
        for (int ic = 0; ic < nc; ic++) {
            
            coord1 = &coord[beadSet[Cumnc + ic]];
            r = be->distance(coord1);
            
            R = -r / slen[Cumnc + ic] + 100.0 / slen[Cumnc + ic];
            U_i = krep[Cumnc + ic] * exp(R);

            U += U_i;
        }
        Cumnc += nc;
    }
    return U;
}



void BoundaryCylinderRepulsionExpIn::forces(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                          floatingpoint *krep, floatingpoint *slen, int *nneighbors) {
    int nb, nc;
    floatingpoint *coord1, R, r, f0;
    floatingpoint *force1;
    //    floatingpoint *forcecopy;
    //    forcecopy = new floatingpoint[CGMethod::N];
    //    for(auto iter=0;iter<CGMethod::N;iter++)
    //        forcecopy[iter]=0.0;
    
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
            
            R = -r / slen[Cumnc + ic] + 100.0 / slen[Cumnc + ic];
            f0 = krep[Cumnc + ic] * exp(R)/ slen[Cumnc + ic];
            force1[0] += f0 *norm[0];
            force1[1] += f0 *norm[1];
            force1[2] += f0 *norm[2];
#ifdef CHECKFORCES_INF_NAN
            if(checkNaN_INF<floatingpoint>(force1, 0, 2)){
                cout<<"Boundary Cylinder Force becomes infinite. Printing data "<<endl;
                
                cout<<"Printing coords"<<endl;
                cout<<coord1[0]<<" "<<coord1[1]<<" "<<coord1[2]<<endl;
                cout<<"Printing force"<<endl;
                cout<<force1[0]<<" "<<force1[1]<<" "<<force1[2]<<endl;
                cout<<"Printing binary Coords"<<endl;
                printvariablebinary(coord1,0,2);
                cout<<"Printing binary Force"<<endl;
                printvariablebinary(force1,0,2);
                exit(EXIT_FAILURE);
            }
#endif
            
        }
        Cumnc+=nc;
    }
}

floatingpoint BoundaryCylinderRepulsionExpIn::loadForces(floatingpoint r, floatingpoint kRep, floatingpoint screenLength) const {
    
    floatingpoint R = -r/screenLength + 100.0 / screenLength;
    return kRep * exp(R)/screenLength;
    
}

} // namespace medyan

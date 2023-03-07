
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "FilamentBendingCosine.h"

#include "MathFunctions.h"
#include "Bead.h"

#ifdef CROSSCHECK
using namespace mathfunc;
floatingpoint FilamentBendingCosine::energy(Bead* b1, Bead* b2, Bead* b3,
                                     floatingpoint kBend, floatingpoint eqTheta){
    
    floatingpoint L1 = sqrt(scalarProduct(b1->coordinate, b2->coordinate,
                                   b1->coordinate, b2->coordinate));
    floatingpoint L2 = sqrt(scalarProduct(b2->coordinate, b3->coordinate,
                                   b2->coordinate, b3->coordinate));
    
    floatingpoint L1L2 = L1*L2;
    floatingpoint l1l2 = scalarProduct(b1->coordinate, b2->coordinate,
                                b2->coordinate, b3->coordinate);
    
    floatingpoint phi = safeacos(l1l2 / L1L2);
    floatingpoint dPhi = phi-eqTheta;

    return kBend * ( 1 - cos(dPhi) );
}

floatingpoint FilamentBendingCosine::energy(Bead* b1, Bead* b2, Bead* b3,
                                     floatingpoint kBend, floatingpoint eqTheta, floatingpoint d ){
    
    floatingpoint L1 = sqrt(scalarProductStretched(b1->coordinate, b1->force,
                                            b2->coordinate, b2->force,
                                            b1->coordinate, b1->force,
                                            b2->coordinate, b2->force, d));
    floatingpoint L2 = sqrt(scalarProductStretched(b2->coordinate, b2->force,
                                            b3->coordinate, b3->force,
                                            b2->coordinate, b2->force,
                                            b3->coordinate, b3->force, d));
    
    floatingpoint L1L2 = L1*L2;
    floatingpoint l1l2 = scalarProductStretched(b1->coordinate, b1->force,
                                         b2->coordinate, b2->force,
                                         b2->coordinate, b2->force,
                                         b3->coordinate, b3->force, d);
    
    floatingpoint phi = safeacos(l1l2 / L1L2);
    floatingpoint dPhi = phi-eqTheta;

    return kBend * ( 1 - cos(dPhi) );
}

void FilamentBendingCosine::forces(Bead* b1, Bead* b2, Bead* b3,
                                   floatingpoint kBend, floatingpoint eqTheta ){
    
    floatingpoint k = 0;
    floatingpoint L1 = sqrt(scalarProduct(b1->coordinate, b2->coordinate,
                                   b1->coordinate, b2->coordinate));
    floatingpoint L2 = sqrt(scalarProduct(b2->coordinate, b3->coordinate,
                                   b2->coordinate, b3->coordinate));
    floatingpoint l1l2 = scalarProduct(b1->coordinate, b2->coordinate,
                                b2->coordinate, b3->coordinate);
    
    //invL = 1/L;
    floatingpoint invL1 = 1/L1;
    floatingpoint invL2 = 1/L2;
    floatingpoint A = invL1*invL2;
    floatingpoint B = l1l2*invL1*A*A*L2;
    floatingpoint C = l1l2*invL2*A*A*L1;
    
    if (areEqual(eqTheta, 0.0)) k = kBend;
    
    else{
        floatingpoint phi = safeacos(l1l2 *A);
        floatingpoint dPhi = phi-eqTheta;
        
        k =  kBend* sin(dPhi)/sin(phi);
    }
    
    //force on i-1, f = k*(-A*l2 + B*l1):
    b1->force[0] +=  k * ((-b3->coordinate[0] + b2->coordinate[0])*A +
                          (b2->coordinate[0] - b1->coordinate[0])*B );
    b1->force[1] +=  k * ((-b3->coordinate[1] + b2->coordinate[1])*A +
                          (b2->coordinate[1] - b1->coordinate[1])*B );
    b1->force[2] +=  k * ((-b3->coordinate[2] + b2->coordinate[2])*A +
                          (b2->coordinate[2] - b1->coordinate[2])*B );
    
    
    //force on i, f = k*(A*(l1-l2) - B*l1 + C*l2):
    b2->force[0] +=  k *( (b3->coordinate[0] - 2*b2->coordinate[0] + b1->coordinate[0])*A -
                         (b2->coordinate[0] - b1->coordinate[0])*B +
                         (b3->coordinate[0] - b2->coordinate[0])*C );
    
    b2->force[1] +=  k *( (b3->coordinate[1] - 2*b2->coordinate[1] + b1->coordinate[1])*A -
                         (b2->coordinate[1] - b1->coordinate[1])*B +
                         (b3->coordinate[1] - b2->coordinate[1])*C );
    
    b2->force[2] +=  k *( (b3->coordinate[2] - 2*b2->coordinate[2] + b1->coordinate[2])*A -
                         (b2->coordinate[2] - b1->coordinate[2])*B +
                         (b3->coordinate[2] - b2->coordinate[2])*C );
    
    //force on i-1, f = k*(A*l - B*l2):
    b3->force[0] +=  k *( (b2->coordinate[0] - b1->coordinate[0])*A -
                         (b3->coordinate[0] - b2->coordinate[0])*C );
    
    b3->force[1] +=  k *( (b2->coordinate[1] - b1->coordinate[1])*A -
                         (b3->coordinate[1] - b2->coordinate[1])*C );
    
    b3->force[2] +=  k *( (b2->coordinate[2] - b1->coordinate[2])*A -
                         (b3->coordinate[2] - b2->coordinate[2])*C );
}

void FilamentBendingCosine::forcesAux(Bead* b1, Bead* b2, Bead* b3,
                                      floatingpoint kBend, floatingpoint eqTheta){
    
    floatingpoint k = 0;
    floatingpoint L1 = sqrt(scalarProduct(b1->coordinate, b2->coordinate,
                                   b1->coordinate, b2->coordinate));
    floatingpoint L2 = sqrt(scalarProduct(b2->coordinate, b3->coordinate,
                                   b2->coordinate, b3->coordinate));
    floatingpoint l1l2 = scalarProduct(b1->coordinate, b2->coordinate,
                                b2->coordinate, b3->coordinate);
    
    //invL = 1/L;
    floatingpoint invL1 = 1/L1;
    floatingpoint invL2 = 1/L2;
    floatingpoint A = invL1*invL2;
    floatingpoint B = l1l2*invL1*A*A*L2;
    floatingpoint C = l1l2*invL2*A*A*L1;
    
    if (areEqual(eqTheta, 0.0)) k = kBend;
    
    else{
        floatingpoint phi = safeacos(l1l2 *A);
        floatingpoint dPhi = phi-eqTheta;
        
        k =  kBend* sin(dPhi)/sin(phi);
    }
    
    //force on i-1, f = k*(-A*l2 + B*l1):
    b1->forceAux[0] +=  k * ( (-b3->coordinate[0] + b2->coordinate[0])*A +
                             (b2->coordinate[0] - b1->coordinate[0])*B );
    
    b1->forceAux[1] +=  k * ( (-b3->coordinate[1] + b2->coordinate[1])*A +
                             (b2->coordinate[1] - b1->coordinate[1])*B );
    
    b1->forceAux[2] +=  k* ( (-b3->coordinate[2] + b2->coordinate[2])*A +
                            (b2->coordinate[2] - b1->coordinate[2])*B );
    
    //force on i, f = k*(A*(l1-l2) - B*l1 + C*l2):
    b2->forceAux[0] +=  k *( (b3->coordinate[0] - 2*b2->coordinate[0] + b1->coordinate[0])*A -
                            (b2->coordinate[0] - b1->coordinate[0])*B +
                            (b3->coordinate[0] - b2->coordinate[0])*C );
    
    b2->forceAux[1] +=  k *( (b3->coordinate[1] - 2*b2->coordinate[1] + b1->coordinate[1])*A -
                            (b2->coordinate[1] - b1->coordinate[1])*B +
                            (b3->coordinate[1] - b2->coordinate[1])*C );
    
    b2->forceAux[2] +=  k *( (b3->coordinate[2] - 2*b2->coordinate[2] + b1->coordinate[2])*A -
                            (b2->coordinate[2] - b1->coordinate[2])*B +
                            (b3->coordinate[2] - b2->coordinate[2])*C );
    
    //force on i-1, f = k*(A*l - B*l2):
    b3->forceAux[0] +=  k *( (b2->coordinate[0] - b1->coordinate[0])*A -
                            (b3->coordinate[0] - b2->coordinate[0])*C );
    
    b3->forceAux[1] +=  k *( (b2->coordinate[1] - b1->coordinate[1])*A -
                            (b3->coordinate[1] - b2->coordinate[1])*C );
    
    b3->forceAux[2] +=  k *( (b2->coordinate[2] - b1->coordinate[2])*A -
                            (b3->coordinate[2] - b2->coordinate[2])*C );
    
}
#endif

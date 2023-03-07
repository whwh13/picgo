
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

#include "FilamentStretchingHarmonic.h"

#include "Bead.h"

#include "MathFunctions.h"

#ifdef CROSSCHECK
using namespace mathfunc;
floatingpoint FilamentStretchingHarmonic::energy(Bead* b1, Bead* b2,
                                          floatingpoint kStretch, floatingpoint eqLength){
    
    floatingpoint dist = twoPointDistance( b1->coordinate, b2->coordinate) - eqLength;
    
    return 0.5 * kStretch* dist * dist;

}

floatingpoint FilamentStretchingHarmonic::energy(Bead* b1, Bead* b2,
                                          floatingpoint kStretch, floatingpoint eqLength, floatingpoint d){
    
    floatingpoint distStretched = twoPointDistanceStretched(b1->coordinate,
                                                     b1->force, b2->coordinate,
                                                     b2->force, d) - eqLength;

    return 0.5 * kStretch * distStretched * distStretched;
}

void FilamentStretchingHarmonic::forces(Bead* b1, Bead* b2,
                                        floatingpoint kStretch, floatingpoint eqLength ){
    floatingpoint dist = twoPointDistance( b1->coordinate, b2->coordinate);
    floatingpoint invL = 1 / dist;
    
    floatingpoint f0 = kStretch * ( dist - eqLength ) * invL;
    
    //force on i
    b2->force[0] +=  f0 * ( b1->coordinate[0] - b2->coordinate[0] );
    b2->force[1] +=  f0 * ( b1->coordinate[1] - b2->coordinate[1] );
    b2->force[2] +=  f0 * ( b1->coordinate[2] - b2->coordinate[2] );
    
    // force i-1
    b1->force[0] +=  f0 * ( b2->coordinate[0] - b1->coordinate[0] );
    b1->force[1] +=  f0 * ( b2->coordinate[1] - b1->coordinate[1] );
    b1->force[2] +=  f0 * ( b2->coordinate[2] - b1->coordinate[2] );
}

void FilamentStretchingHarmonic::forcesAux(Bead* b1, Bead* b2,
                                           floatingpoint kStretch, floatingpoint eqLength ){
    
    floatingpoint dist = twoPointDistance( b1->coordinate, b2->coordinate);
    floatingpoint invL = 1 / dist;
    floatingpoint f0 = kStretch * ( dist - eqLength ) * invL;
    
    //force on i
    b2->forceAux[0] +=  f0 * ( b1->coordinate[0] - b2->coordinate[0] );
    b2->forceAux[1] +=  f0 * ( b1->coordinate[1] - b2->coordinate[1] );
    b2->forceAux[2] +=  f0 * ( b1->coordinate[2] - b2->coordinate[2] );
    
    // force i-1
    b1->forceAux[0] +=  f0 * ( b2->coordinate[0] - b1->coordinate[0] );
    b1->forceAux[1] +=  f0 * ( b2->coordinate[1] - b1->coordinate[1] );
    b1->forceAux[2] +=  f0 * ( b2->coordinate[2] - b1->coordinate[2] );
}
#endif

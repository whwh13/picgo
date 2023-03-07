
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

#include <cmath>
#include <algorithm>

#include "RateChangerImpl.h"

#include "SysParams.h"

namespace medyan {
float BrownianRatchet::changeRate(float bareRate, floatingpoint force) {

    force = min<floatingpoint>(force, (floatingpoint)100.0); //ceiling

    floatingpoint newRate = bareRate * exp( - force * _x / kT);
    
    return newRate;
}

float BrownianRatchet::getRateChangeFactor(floatingpoint force) {

    force = min<floatingpoint>(force, (floatingpoint)100.0); //ceiling

    return exp( - force * _x / kT);
}


float LinkerCatchSlip::changeRate(float bareRate, floatingpoint force) {
    
    return bareRate * (_a1 * exp(- force * _x1 / kT) +
                       _a2 * exp(  force * _x2 / kT));
}

float LinkerCatchSlip::getRateChangeFactor(floatingpoint force) {

    return (_a1 * exp(- force * _x1 / kT) +
                       _a2 * exp(  force * _x2 / kT));
}

float LinkerSlip::changeRate(float bareRate, floatingpoint force) {

    floatingpoint newRate = bareRate * exp( force * _x / kT);
    
    return newRate;
}

float LinkerSlip::getRateChangeFactor(floatingpoint force) {

    return exp( force * _x / kT);
}

float BranchSlip::changeRate(float bareRate, floatingpoint force) {

    floatingpoint newRate = bareRate * exp( force * _x / kT);
    
    return newRate;
}

float BranchSlip::getRateChangeFactor( floatingpoint force) {

    return exp( force * _x / kT);
}

float BranchSlipF::changeRate(float bareRate, floatingpoint force) {

	floatingpoint newRate = bareRate * exp( force/_F0);

	return newRate;
}

float BranchSlipF::getRateChangeFactor( floatingpoint force) {

	return exp( force / _F0);
}

float MotorCatch::numBoundHeads(float onRate, float offRate,
                                floatingpoint force, int numHeads) {

    _dutyRatio = onRate/(onRate+offRate);
#ifdef PLOSFEEDBACK
    return min<floatingpoint >((floatingpoint)numHeads,
    		numHeads * _dutyRatio + _gamma * force);
#else
    return min<floatingpoint>((floatingpoint)numHeads, numHeads * _dutyRatio + _beta * force / numHeads);
#endif
}

float MotorCatch::changeRate(float onRate, float offRate,
                             floatingpoint numHeads, floatingpoint force) {
//    cout<<"offRate "<<offRate<<endl;
//    return offRate;
    
    //calculate new rate
#ifdef PLOSFEEDBACK
    floatingpoint k_0 = _beta * offRate /numBoundHeads(onRate, offRate, force, numHeads);

    floatingpoint factor = exp(-force / (numBoundHeads(onRate, offRate, force, numHeads) * _F0));
#else
    floatingpoint k_0 = onRate * (numHeads) / (exp(std::log((onRate + offRate) / offRate) * numHeads)
                                        - (floatingpoint)1.0);

    floatingpoint factor = max<floatingpoint>((floatingpoint)0.1, exp(-force / (numBoundHeads(onRate,
    		offRate, force,
    		numHeads) * _F0)));
#endif
    
    floatingpoint newRate = k_0 * factor;

    return newRate;
}

float MotorStall::changeRate(float onRate, float offRate,
                             floatingpoint numHeads, floatingpoint force) {
//    cout<<"onRate "<<onRate<<endl;
//    return onRate*10;
    _dutyRatio = onRate/(onRate+offRate);
    //determine k_0
    float k_0 = ((1 - _dutyRatio) / _dutyRatio) * onRate * _stepFrac;

#if defined(PLOSFEEDBACK) || defined(PLOSSTALLFEEDBACK)
    floatingpoint newRate =  max((floatingpoint(0.0)), k_0 * (_F0 - force/numHeads)
                          / (_F0 + (force / (numHeads * _alpha))));
#else
    floatingpoint newRate =  max<floatingpoint>((floatingpoint)0.0, k_0 * (_F0 - force)
                               / (_F0 + (force / (_alpha))));
#endif

    return newRate;
}

} // namespace medyan

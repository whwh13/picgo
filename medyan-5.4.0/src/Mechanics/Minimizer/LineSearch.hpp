#ifndef MEDYAN_Mechanics_Minimizer_LineSearch_hpp
#define MEDYAN_Mechanics_Minimizer_LineSearch_hpp

#include <optional>
#include <vector>

#include "Mechanics/CUDAcommon.h"
#include "Rand.h"

namespace medyan {

struct LineSearchResult {
    // If lambda is zero, then the line search failed.
    floatingpoint lambda = 0;

    int numEnergyCall = 0;
    int numForceCall = 0;
    floatingpoint newEnergy = 0;
};

struct LambdaRunningAverageManager {
    // Track the past 100 lambdas.
    //@{
    floatingpoint lambdaRunningAverageProbability = 0;
    unsigned maxprevlambdacount = 10;
    std::vector<floatingpoint> previouslambdavec = std::vector<floatingpoint>(maxprevlambdacount,0.0);
    short headpos = 0; //position where the next lambda can be inserted.
    float sum = 0;//sum of the lambdas found in previouslambdavcec.
    bool runningaveragestatus = false;
    //@}

    void record(floatingpoint lambda) {
        //@{
	    sum = sum - previouslambdavec[headpos] + lambda;
        previouslambdavec[headpos] = lambda;

        if(headpos==maxprevlambdacount-1) {
	        headpos = 0;
            runningaveragestatus = Rand::randfloatingpoint(0,1) < lambdaRunningAverageProbability;
        }
        else
        	headpos++;
    }

    floatingpoint suggestLambdaMax() const {
        return sum / maxprevlambdacount;
    }
};

} // namespace medyan

#endif

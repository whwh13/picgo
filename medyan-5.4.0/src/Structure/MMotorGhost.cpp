
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

#include "Structure/MMotorGhost.h"

#include "SysParams.h"

namespace medyan {
void MMotorGhost::initializerestart(int motorType, floatingpoint eqLength,
		floatingpoint numBoundHeads){

    if(SysParams::RUNSTATE){
        LOG(ERROR) << "initializerestart Function from MMotorGhost class can only be "
                      "called "
                      "during restart phase. Exiting.";
        throw std::logic_error("Illegal function call pattern");
    }

    if(numBoundHeads > 0) {
        if(!SysParams::Mechanics().MStretchingK.empty())
            kStretch = SysParams::Mechanics().MStretchingK[motorType] * numBoundHeads;
    }
    this->eqLength = eqLength;
}

void MMotorGhost::setStretchingConstant(int motorType, floatingpoint numBoundHeads) {
    
    if(!SysParams::Mechanics().MStretchingK.empty())
        kStretch = SysParams::Mechanics().MStretchingK[motorType] * numBoundHeads;
}

} // namespace medyan

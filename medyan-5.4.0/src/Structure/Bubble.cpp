
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

#include "Bubble.h"

#include "SubSystem.h"
#include "Bead.h"

#include "SysParams.h"
#include "CUDAcommon.h"

namespace medyan {

void Bubble::printSelf()const {
    
    cout << endl;
    
    cout << "Bubble: ptr = " << this << endl;
    cout << "Bubble ID = " << getId() << endl;
    cout << "Bubble type = " << type_ << endl;
    cout << "Bubble coord = " << coord << endl;
    cout << "Bubble force = " << force << endl;
    cout << "Bubble radius = " << _radius << endl;
    
    cout << endl;
}

} // namespace medyan

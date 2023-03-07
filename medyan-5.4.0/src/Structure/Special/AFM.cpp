
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

#include "AFM.h"

#include "Bubble.h"
#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"

namespace medyan {

void AFM::printSelf() const {
    
    cout << endl;
    cout << "AFM Bubble index " << bubbleSysIndex_.value << endl;
    
    cout << endl;
}

} // namespace medyan

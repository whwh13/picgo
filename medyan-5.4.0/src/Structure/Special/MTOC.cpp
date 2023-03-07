
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

#include "MTOC.h"

#include "Bubble.h"
#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"

namespace medyan {

void MTOC::printSelf()const {
    
    cout << endl;
    cout << "MTOC Bubble index " << bubbleSysIndex_.value << endl;
    
    cout << endl;
}

} // namespace medyan

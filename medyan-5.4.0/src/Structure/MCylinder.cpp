
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

#include "MCylinder.h"

#include "MathFunctions.h"
#include "SysParams.h"

namespace medyan {
using namespace mathfunc;

MCylinder::MCylinder(short filamentType, floatingpoint eqLength){
    
    //Set equilibrium length and constants relative to full cylinder length
    setEqLength(filamentType, eqLength);
    //set excluded volume const
    if(!SysParams::Mechanics().VolumeK.empty())
        _kExVol = SysParams::Mechanics().VolumeK[filamentType];
    
    //set angle
    if(!SysParams::Mechanics().FBendingTheta.empty())
        _eqTheta = SysParams::Mechanics().FBendingTheta[filamentType];
}

void MCylinder::setEqLength(short filamentType, floatingpoint l) {
    _eqLength = l;
    floatingpoint fracCylinderSize = SysParams::Geometry().cylinderSize[filamentType] / l;
    // recalculate other constants
    if(!SysParams::Mechanics().FStretchingK.empty())
        _kStretch = SysParams::Mechanics().FStretchingK[filamentType] * fracCylinderSize;
    if(!SysParams::Mechanics().FBendingK.empty())
        _kBend = SysParams::Mechanics().FBendingK[filamentType] * fracCylinderSize;
}

} // namespace medyan

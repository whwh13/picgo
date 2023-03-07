
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

#include "MBranchingPoint.h"

#include "SysParams.h"
#include "Util/Io/Log.hpp"

namespace medyan {
MBranchingPoint::MBranchingPoint(int branchType) {
    
    //set parameters
    if(!SysParams::Mechanics().BrStretchingK.empty()) {
        _kStretch = SysParams::Mechanics().BrStretchingK[branchType];
        _eqLength = SysParams::Mechanics().BrStretchingL[branchType];
    }
    
    if(!SysParams::Mechanics().BrBendingK.empty()) {
        _kBend = SysParams::Mechanics().BrBendingK[branchType];
        _eqTheta = SysParams::Mechanics().BrBendingTheta[branchType];
    }
    
    if(!SysParams::Mechanics().BrDihedralK.empty())
        _kDihedr = SysParams::Mechanics().BrDihedralK[branchType];
 
    if(!SysParams::Mechanics().BrPositionK.empty())
        _kPosition = SysParams::Mechanics().BrPositionK[branchType];
}

void MBranchingPoint::initializerestart(floatingpoint eqLength){
    if(SysParams::RUNSTATE){
        LOG(ERROR) << "initializerestart Function from MBranchingPoint class can only be "
                      "called "
                      "during restart phase. Exiting.";
        throw std::logic_error("Illegal function call pattern");
    }
    _eqLength = eqLength;}

} // namespace medyan

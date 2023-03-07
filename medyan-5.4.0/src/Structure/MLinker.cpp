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

#include <stdexcept>

#include "Structure/MLinker.h"

#include "SysParams.h"
#include "Util/Io/Log.hpp"

namespace medyan {
void MLinker::initializerestart(floatingpoint eqLength){
    if(SysParams::RUNSTATE){
        LOG(ERROR) << "initializerestart Function from MLinker class can only be called "
                      "during restart phase. Exiting.";
        throw std::logic_error("Illegal function call pattern");
    }
    this->eqLength = eqLength;
}

} // namespace medyan

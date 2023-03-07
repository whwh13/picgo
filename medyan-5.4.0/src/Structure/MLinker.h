
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

#ifndef MEDYAN_MLinker_h
#define MEDYAN_MLinker_h

#include "common.h"

namespace medyan {
/// Represents the mechanical component of a Linker.

/*! The class describes interaction between 4 [Beads](@ref Bead) connected by a Linker,
 *  and its associated equilibrium constants. Initial length of a Linker is determined
 *  by the condition of zero initial stress, i.e., it calculated within the constructor
 *  at initiation. A Linker heads positions on a segment (between two consecutive beads 
 *  on a Filament) determined by two numbers (floatingpoint from 0 to 1) position1 and 
 *  position2 held in the parent Linker.
 */
struct MLinker {
    
    floatingpoint stretchForce = 0.0; ///< Stretching force of linker at current state
    
    floatingpoint eqLength = 0;  ///< Equilibrium length
    floatingpoint kStretch = 0;  ///< Stretching constant
    
    //@{
    /// Getter for constants
    floatingpoint getStretchingConstant() const {return kStretch;}
    floatingpoint getEqLength() const {return eqLength;}
    //@}

    void initializerestart(floatingpoint eqLength);
};

} // namespace medyan

#endif 

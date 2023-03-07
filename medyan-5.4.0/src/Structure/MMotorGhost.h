
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


#ifndef MEDYAN_MMotorGhost_h
#define MEDYAN_MMotorGhost_h

#include "SysParams.h"
#include "Util/Io/Log.hpp"
#include "common.h"

namespace medyan {
/// Represents a cross-link between [Filaments](@ref Filament) that can move by way of
/// chemical reactions.

/*!
 *  The class describes interaction between 4 [Beads](@ref Bead) connected by a 
 *  MotorGhost, and equilibrium constants. Initial length is determinated by the condition 
 *  of zero initial stress, i.e., it is calculated within the constructor at initiation. 
 *  A ghost motor heads positions on a segment (between two consecutive beads on a filament) 
 *  determined by two numbers (0 to 1) position1 and position2 (finite number of steps 
 *  before move to the next segment o--x-o- -> o---xo- -> o---ox) held in the parent
 *  MotorGhost. They can be changed as a result of chemical reaction, and we consider that 
 *  as the motor making a step.
 */
class MMotorGhost {
    
public:
    floatingpoint stretchForce = 0.0; ///< Stretching force of motor at current state

    floatingpoint eqLength = 0;  ///< Equilibrium length
    floatingpoint kStretch = 0;  ///< Stretching parameter

    void initializerestart(int motorType, floatingpoint eqLength,
                           floatingpoint numBoundHeads);
    
    //@{
    /// Getter for constants
    floatingpoint getStretchingConstant() const {return kStretch;}
    floatingpoint getEqLength() const {return eqLength;}
    //@}
    
    /// Reset the spring constant of the motor based on number of bound heads
    void setStretchingConstant(int motorType, floatingpoint numBoundHeads);
};

} // namespace medyan

#endif

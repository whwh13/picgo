
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

#ifndef MEDYAN_RateChanger_h
#define MEDYAN_RateChanger_h

#include "common.h"

namespace medyan {
/// Used to change Filament reaction rates based on forces in the network
/*!
 *  The FilamentRateChanger class is an abstract class which allows
 *  for Filament rate changing based on a given force. Different
 *  implementations of this class will have different rate changing models, 
 *  and will all implement the changeRate() function.
 */

class FilamentRateChanger {
    
public:
    /// Change the reaction rate based on a bare rate and given force.
    virtual float changeRate(float bareRate, floatingpoint force) = 0;

    virtual float getRateChangeFactor(floatingpoint force) = 0;
};

/// Used to change Linker reaction rates based on forces in the network
/*!
 *  The LinkerRateChanger class is an abstract class which allows
 *  for Linker rate changing based on a given force. Different 
 *  implementations of this class will have different rate changing models, 
 *  and will all implement the changeRate() function.
 */
class LinkerRateChanger {
    
protected:
    short _linkerType; ///< This linker type
    
public:
    LinkerRateChanger(short linkerType) : _linkerType(linkerType) {}
    
    /// Change the reaction rate based on a bare rate and given force.
    virtual float changeRate(float bareRate, floatingpoint force) = 0;

    virtual float getRateChangeFactor(floatingpoint force) = 0;
};

/// Used to change Linker reaction rates based on forces in the network
/*!
 *  The LinkerRateChanger class is an abstract class which allows
 *  for Linker rate changing based on a given force. Different
 *  implementations of this class will have different rate changing models,
 *  and will all implement the changeRate() function.
 */
class BranchRateChanger {
    
protected:
    short _branchType; ///< This linker type
    
public:
    BranchRateChanger(short branchType) : _branchType(branchType) {}
    
    /// Change the reaction rate based on a bare rate and given force.
    virtual float changeRate(float bareRate, floatingpoint force) = 0;

    virtual float getRateChangeFactor(floatingpoint force) = 0;
};

/// Used to change MotorGhost reaction rates based on forces in the network
/*!
 *  The MotorRateChanger class is an abstract class which allows
 *  for MotorGhost rate changing based on a given force. Different 
 *  implementations of this class will have different rate changing models, 
 *  and will all implement the changeRate() function.
 */
class MotorRateChanger {
    
protected:
    short _motorType; ///< This motor type
    
public:
    MotorRateChanger(short motorType) : _motorType(motorType) {}
    
    /// Calculate the number of bound heads in the ensemble, dependent on force
    virtual float numBoundHeads(float onRate, float offRate,
                                floatingpoint force, int numHeads) = 0;
    
    /// Change the reaction rate based on an on rate, off rate,
    /// number of heads, and given force.
    virtual float changeRate(float onRate, float offRate,
                             floatingpoint numBoundHeads, floatingpoint force) = 0;
};

} // namespace medyan

#endif

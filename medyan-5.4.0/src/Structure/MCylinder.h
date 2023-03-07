
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

#ifndef MEDYAN_MCylinder_h
#define MEDYAN_MCylinder_h

#include <vector>

#include "common.h"

namespace medyan {
//FORWARD DECLARATIONS
class Cylinder;

/// Used to hold mechanical properties of a Cylinder.

/*! The class describes interaction between 2 [Beads](@ref Bead) connected by a
 *  Cylinder in the same Filament, and its associated equilibrium constants. 
 */
class MCylinder {

private:
    Cylinder* _pCylinder;  ///< parent cylinder

    floatingpoint _eqLength; ///< Length of unstretched cylinder
    floatingpoint _eqTheta;  ///< Equilibrium value for angle in bending potential.
                      ///< For interaction between this cylinder and PREVIOUS
    floatingpoint _eqPhi;    ///< Equilibrium value of twisiting potential
    floatingpoint _kStretch; ///< Local stretching constant, describes axial stretching
                      ///< of a single cylinder
    floatingpoint _kBend;    ///< Local bending constant, which describes bending
                      ///< interaction between current and PREVIOUS cylinders
    floatingpoint _kTwist;   ///< Local twisting constant, which describes stretching
                      ///< interaction between current and PREVIOUS cylinders
    floatingpoint _kExVol;   ///< Local excluded volume constant, which describes
                      ///< excluded volume interactions between cylinders
    
    floatingpoint _currentLength; ///< The current length of the cylinder
    
public:
    /// Constructor sets equlilibrium length, and also adjusts other
    /// parameters according to this length
    MCylinder(short filamentType, floatingpoint eqLength);
    ~MCylinder() {};

    /// Set parent 
    void setCylinder(Cylinder* c) {_pCylinder = c;}
    Cylinder* getCylinder() {return _pCylinder;}
    
    /// Set the equlilibrium length, which changes mechanical constants accordingly
    void setEqLength(short filamentType, floatingpoint l);
    /// Get the current equlibrium length of this MCylinder
    floatingpoint getEqLength() {return _eqLength;}
    
    //@{
    /// Mechanical parameter management function
    void setEqTheta(floatingpoint theta) {_eqTheta = theta;}
    floatingpoint getEqTheta() {return _eqTheta;}
    
    void setEqPhi(floatingpoint phi) {_eqPhi = phi;}
    floatingpoint getEqPhi() {return _eqPhi;}
    
    void setStretchingConst(floatingpoint k) {_kStretch = k;}
    floatingpoint getStretchingConst() {return _kStretch;}
    
    void setBendingConst(floatingpoint k) {_kBend = k;}
    floatingpoint getBendingConst() {return _kBend;}
    
    void setTwistingConst(floatingpoint k) {_kTwist = k;}
    floatingpoint getTwistingConst() {return _kTwist;}
    
    void setExVolConst(floatingpoint k) {_kExVol = k;}
    floatingpoint getExVolConst() const {return _kExVol;}
    
    void setLength(floatingpoint l){_currentLength = l;}
    floatingpoint getLength() {return _currentLength;}
    //@}
    
};

} // namespace medyan

#endif

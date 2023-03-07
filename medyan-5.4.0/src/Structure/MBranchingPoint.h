
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

#ifndef MEDYAN_MBranchingPoint_h
#define MEDYAN_MBranchingPoint_h

#include "common.h"
#include "Util/Math/Vec.hpp"

namespace medyan {
//FORWARD DECLARATIONS
class BranchingPoint;

/// Represents the mechanical component of a BranchingPoint.

/*! The class describes interaction between 4 [Beads](@ref Bead) connected by a 
 *  BranchingPoint, and its associated equilibrium constants. The branch is a connection 
 *  between two [Cylinders](@ref Cylinder), where one resides on the mother Filament and
 *  the other is the actual branch. This is positioned on a Cylinder with the field
 *  position (between 0 and 1), held by the parent BranchingPoint.
 */
class MBranchingPoint {
    
public:

    /// Main constructor, sets constants
    MBranchingPoint(int branchType);
    
    //@{
    /// Getter for constants
    floatingpoint getStretchingConstant(){return _kStretch;}
    floatingpoint getEqLength(){return _eqLength;}
    
    floatingpoint getBendingConstant(){return _kBend;}
    floatingpoint getEqTheta(){return _eqTheta;}
    
    floatingpoint getDihedralConstant(){return _kDihedr;}
    
    floatingpoint getPositionConstant(){return _kPosition;}
    //@}
    
    /// Set parent
    void setBranchingPoint(BranchingPoint* BranchingPoint) {
        _pBranchingPoint = BranchingPoint;
    }
    /// Get parent
    BranchingPoint* getBranchingPoint() {return _pBranchingPoint;}

    void initializerestart(floatingpoint eqLength);

    Vec<3, FP> branchForce = {0.0, 0.0, 0.0};

private:
    floatingpoint _eqLength;  ///< Equilibrium length
    floatingpoint _kStretch;  ///< Stretching constant
    
    floatingpoint _eqTheta; ///< Bending equilibrium angle
    floatingpoint _kBend; //< Bending constant
    
    floatingpoint _kDihedr; ///< Twisting constant
    
    floatingpoint _kPosition; ///< Position constant
    
    BranchingPoint* _pBranchingPoint; ///< Pointer to parent branch point
    
};

} // namespace medyan

#endif

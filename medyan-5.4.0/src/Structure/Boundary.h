
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

#ifndef MEDYAN_Boundary_h
#define MEDYAN_Boundary_h

#include "common.h"

#include "BoundarySurface.h"

namespace medyan {

//FORWARD DECLARATIONS
class Compartment;
class SubSystem;

/// BoundaryShape is a shape enumeration.
enum class BoundaryShape {Cube, Capsule, Sphere, Cylinder};

/// BoundaryMove is a enum describing the movement of a boundary.
enum class BoundaryMove {None, Top, Bottom, Left, Right, Front, Back, All};

/// To store all [BoundarySurfaces](@ref BoundarySurface) that are in the SubSystem.
/*!
 *  The boundary class stores all [BoundarySurfaces](@ref BoundarySurface) in the given 
 *  shape. Its constructors can create basic boundary shapes (for now). Eventually will 
 *  be extended to more complex surfaces.
 */
class Boundary {
    
protected:
    SubSystem* _subSystem; ///< SubSystem ptr
    
    /// Vector of boundarysurfaces (could be different implementations)
    vector<unique_ptr<BoundarySurface>> _boundarySurfaces;
    
    BoundaryShape _shape; ///< Shape of boundary
    vector<BoundaryMove> _move;   ///< Movement of boundary
    
    short _nDim; ///< Dimensionality
    
public:
    Boundary(SubSystem* s, int nDim, BoundaryShape shape, vector<BoundaryMove> move)
        : _subSystem(s), _shape(shape), _move(move), _nDim(nDim) {};
    
    ~Boundary() {};

    /// Get shape of this boundary
    BoundaryShape getShape() {return _shape;}
    
    /// Get boundary surfaces
    const vector<unique_ptr<BoundarySurface>>& getBoundarySurfaces() {
        return _boundarySurfaces;
    }
    
    /// Check if coordinates are within boundary
    virtual bool within(const vector<floatingpoint>& coordinates) = 0;
    
    /// Check if a compartment is within boundary
    /// @note - this checks if ANY part of the compartment volume
    ///         is within the boundary.
    virtual bool within(Compartment* C) = 0;
    
    /// Get the distance from the boundary. Returns the distance from
    /// closest boundary element in the boundary.
    /// Will return infinity if outside of the boundary.
    virtual floatingpoint distance(const vector<floatingpoint>& coordinates) = 0;

    // Returns the distance from the boundary element in the lower boundary
    virtual floatingpoint lowerdistance(const vector<floatingpoint>& coordinates) = 0;
    // Returns the distance from the boundary element in the side boundary
    virtual floatingpoint sidedistance(const vector<floatingpoint>& coordinates) = 0;
    virtual floatingpoint getboundaryelementcoord(int bidx) = 0;
    ///Move a given part of a boundary a given distance
    ///@note a negative distance denotes movement towards the center of the grid.
    virtual void move(vector<floatingpoint> dist) = 0;
    
    //Give a normal to the plane (pointing inward) at a given point
    virtual vector<floatingpoint> normal(vector<floatingpoint>& coordinates) = 0;

    //returns volume of the enclosed volume. Note. If there are moving boundaries, this
    // volume MAY not the same as volume specified in systeminputfile.
    virtual void volume() = 0;

    inline static floatingpoint systemvolume = 0;
    
};

} // namespace medyan

#endif


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

#ifndef MEDYAN_BoundaryImpl_h
#define MEDYAN_BoundaryImpl_h

#include <stdexcept> // logic_error
#include <vector>
#include <limits>

#include "common.h"

#include "Boundary.h"
#include "Util/Io/Log.hpp"

namespace medyan {
///FORWARD DECLARATIONS
class Compartment;

/// A cubic Boundary implementation.
class BoundaryCubic: public Boundary {
    
public:
    ///Default constructor, this will create a cube with given
    ///corners at edges of current CompartmentGrid
    BoundaryCubic(SubSystem* s, vector<BoundaryMove> move);
    
    virtual bool within(Compartment* C);
    virtual bool within(const vector<floatingpoint>& coordinates);

    virtual floatingpoint distance(const vector<floatingpoint>& coordinates);

    virtual floatingpoint lowerdistance(const vector<floatingpoint>& coordinates);
    virtual floatingpoint sidedistance(const vector<floatingpoint>& coordinates);

    virtual floatingpoint getboundaryelementcoord(int i);

    virtual void move(vector<floatingpoint> dist);
    
    ///Returns the normal inward at this coordinate
    //rule - takes the closest wall's normal inward.
    virtual vector<floatingpoint> normal(vector<floatingpoint>& coordinates);

    virtual void volume();
};

/// A spherical Boundary implementation.
class BoundarySpherical: public Boundary {
    
public:
    ///Default constructor, will create an sphere with given diameter
    ///@param diameter - diameter of sphere
    BoundarySpherical(SubSystem* s, floatingpoint diameter, vector<BoundaryMove> move);
    
    ///@note - not yet implemented correctly. Essentially checks
    ///        if the midpoint of the compartment is within the boundary.
    virtual bool within(Compartment* C);
    virtual bool within(const vector<floatingpoint>& coordinates);

    virtual floatingpoint distance(const vector<floatingpoint>& coordinates);
    
    //Qin
    virtual floatingpoint lowerdistance(const vector<floatingpoint>& coordinates);
    virtual floatingpoint sidedistance(const vector<floatingpoint>& coordinates);
    virtual floatingpoint getboundaryelementcoord(int i) {
        LOG(ERROR) << "Function is not implemented.";
        throw std::logic_error("Function not implemented");
    }
    ///@note - not yet implemented.
    virtual void move(vector<floatingpoint> dist) {}
    
    ///Returns the normal inward at this coordinate
    virtual vector<floatingpoint> normal(vector<floatingpoint>& coordinate);

    virtual void volume();
};

/// A capsule Boundary implementation.
class BoundaryCapsule: public Boundary {
    
public:
    /// Default constructor, will create a capsule with given diameter, and height equal
    /// to current grid.
    /// @param diameter - diameter of capsule (will set half sphere radii as well as
    /// cylinder radius)
    BoundaryCapsule(SubSystem* s, floatingpoint diameter, vector<BoundaryMove> move);
    
    ///@note - not yet implemented correctly. Essentially checks
    ///        if the midpoint of the compartment is within the boundary.
    virtual bool within(Compartment* C);
    virtual bool within(const vector<floatingpoint>& coordinates);

    virtual floatingpoint distance(const vector<floatingpoint>& coordinates);
    virtual floatingpoint getboundaryelementcoord(int i) {
        LOG(ERROR) << "Function is not implemented.";
        throw std::logic_error("Function not implemented");
    }

    virtual floatingpoint lowerdistance(const vector<floatingpoint>& coordinates);
    virtual floatingpoint sidedistance(const vector<floatingpoint>& coordinates);
    
    ///@note - Not yet implemented.
    virtual void move(vector<floatingpoint> dist) {}
    
    ///Returns the normal inward at this coordinate
    //@note - Not yet implemented.
    virtual vector<floatingpoint> normal(vector<floatingpoint>& coordinate) {return vector<floatingpoint>{0,0,0};}

    virtual void volume();
};

/// A cylinder Boundary implementation.
class BoundaryCylinder: public Boundary {
    
public:
    /// Default constructor, will create a capsule with given diameter, and height equal
    /// to current grid.
    /// @param diameter - diameter of capsule (will set half sphere radii as well as
    /// cylinder radius)
    BoundaryCylinder(SubSystem* s, floatingpoint diameter, vector<BoundaryMove> move);
    
    ///@note - not yet implemented correctly. Essentially checks
    ///        if the midpoint of the compartment is within the boundary.
    virtual bool within(Compartment* C);
    virtual bool within(const vector<floatingpoint>& coordinates);

    virtual floatingpoint distance(const vector<floatingpoint>& coordinates);
    virtual floatingpoint getboundaryelementcoord(int i) {
        LOG(ERROR) << "getboundaryelementcoord Function is not implemented.";
        throw std::logic_error("getboundaryelementcoord Function not implemented");
    }
    //Qin
    virtual floatingpoint lowerdistance(const vector<floatingpoint>& coordinates);
    virtual floatingpoint sidedistance(const vector<floatingpoint>& coordinates);
    
    ///@note - Not yet implemented.
    virtual void move(vector<floatingpoint> dist) {}
    
    ///Returns the normal inward at this coordinate
    //@note - Not yet implemented.
    virtual vector<floatingpoint> normal(vector<floatingpoint>& coordinate) {
        LOG(ERROR) << "normal Function is not implemented.";
        throw std::logic_error("normal Function not implemented");
        return vector<floatingpoint>{0,0,0};}

    virtual void volume();
};

} // namespace medyan

#endif

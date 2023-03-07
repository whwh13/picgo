
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

#ifndef MEDYAN_Trackable_h
#define MEDYAN_Trackable_h

#include "common.h"

namespace medyan {

//FORWARD DECLARATIONS
class SubSystem;
class Controller;

/// An abstract base class for a trackable object in the SubSystem.

/*!
 *  Every class extending Trackable will provide a container to track, 
 *  add, and remove instances of its class from the SubSystem. In general,
 *  this is done using the Database class.
 *
 *  A class extending Trackable can either be Movable, Reactable, DynamicNeighbor, or Neighbor.
 *  These concrete subclasses will implement the boolean functions which define which 
 *  extension of this class is being used.
 *  Dynamic neighbor refers to CylinderCylinderNL while Neighbor refers to
 *  BoundaryCylinderNL
 *  Filament _dneighbor false _neighbor false
 *  Bead _dneighbor true _neighbor false
 *  Cylinder _dneighbor true _neighbor false
 *  BoundaryElement _dneighbor false _neighbor true
 */
class Trackable {
    
friend class Controller;
friend class SubSystem;
    
protected:
    /// Constructor sets boolean values
    Trackable(bool movable   = false,
              bool reactable = false,
              bool dneighbor = false,
              bool neighbor  = false)
    
        : _movable(movable), _reactable(reactable),
          _dneighbor(dneighbor), _neighbor(neighbor) { };
    
    inline static SubSystem* _subSystem = nullptr; ///< A subsystem pointer for every trackable
    
    //@{
    /// Object type
    bool _movable;
    bool _reactable;
    bool _dneighbor;
    bool _neighbor;
    //@}
    
public:
    //@{
    /// Add or remove this Trackable element from the SubSystem.
    /// @note - this update could be due to an updated force minimization,
    /// a set of chemical steps, or any other event in the SubSystem.
    virtual void addToSubSystem() = 0;
    virtual void removeFromSubSystem() = 0;
    //@}
    
    ///Destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~Trackable() noexcept {}
};

} // namespace medyan

#endif


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

#ifndef MEDYAN_BoundaryElement_h
#define MEDYAN_BoundaryElement_h

#include <iostream>

#include "common.h"

#include "Database.h"
#include "Trackable.h"
#include "Neighbor.h"
#include "Component.h"

namespace medyan {

/// Represents an element of a BoundarySurface.
/*!
 * The BoundaryElement class is a representation of a BoundarySurface element, which can
 * interact with other elements in the system, including other BoundaryElements as well 
 * as [Beads] (@ref Bead) in [Filaments](@ref Filament) and [Bubbles](@ref Bubble). 
 * Together, a collection of boundary elements make up a BoundarySurface. 
 *
 * Extending the Neighbor class, all instances can be kept in 
 * [NeighborLists](@ref NeighborList).
 */
class BoundaryElement : public Component, public Trackable, public Neighbor,
    public Database< BoundaryElement, false > {

friend class BoundaryCubic;
friend class BoundarySpherical;
friend class BoundaryCapsule;
friend class BoundaryCylinder;

protected:
    
    vector<floatingpoint> _coords; ///< coordinates
    
    floatingpoint _kRep; ///< Repulsion constant
    floatingpoint _r0; ///< Screening length
    
public:
    /// Default constructor
    BoundaryElement(vector<floatingpoint> coords, floatingpoint kRepuls, floatingpoint screenLength)
    
        : Trackable(false, false, false, true),
          _coords(coords), _kRep(kRepuls), _r0(screenLength) {}
    
    /// Destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~BoundaryElement() noexcept {}
    
    ///return coordinates of boundary element
    const vector<floatingpoint>& getCoords() {return _coords;}
    
    ///update the coordinates of the boundary element
    virtual void updateCoords(const vector<floatingpoint> newCoords) = 0;
    
    //@{
    /// Implement for all boundary elements
    /// Returns the distance from a given point to this boundary element
    /// @return - 1) positive number if point is within boundary element
    ///           2) Negative number if point is outside boundary element
    ///           3) Infinity if point is not in domain of this boundary element
    virtual floatingpoint distance(const vector<floatingpoint>& point) = 0;
    virtual floatingpoint distance(floatingpoint const *point) = 0;
    //@}

    //@{
    virtual floatingpoint lowerdistance(const vector<floatingpoint>& point) = 0;
    virtual floatingpoint sidedistance(const vector<floatingpoint>& point) = 0;

    /// Returns stretched distance, similar to distance above
    virtual floatingpoint stretchedDistance(const vector<floatingpoint>& point,
                                     const vector<floatingpoint>& force,
                                     floatingpoint d) = 0;
    virtual floatingpoint stretchedDistance(floatingpoint const *point,
                                     floatingpoint const *force, floatingpoint d)
    = 0;
    //@}


    //@{
    /// Returns normal vector of point to plane
    virtual const vector<floatingpoint> normal(const vector<floatingpoint> &point) = 0;
    virtual const vector<floatingpoint> normal(const floatingpoint *point) = 0;
    virtual const void elementeqn(floatingpoint* var) = 0;
    //@}

    //@{
    /// Getter for mechanical parameters
    virtual floatingpoint getRepulsionConst() {return _kRep;}
    virtual floatingpoint getScreeningLength() {return _r0;}
    //@}
    
    //@{
    /// SubSystem management, inherited from Trackable
    // Does nothing
    virtual void addToSubSystem() { }
    virtual void removeFromSubSystem() {}
    //@}
    
    /// Get all instances of this class from the SubSystem
    static const vector<BoundaryElement*>& getBoundaryElements() {
        return getElements();
    }
    /// Get the number of boundary elements in this system
    static int numBoundaryElements() {
        return getElements().size();
    }
    
    virtual void printSelf()const;
    
    //GetType implementation just returns zero (no boundary element types yet)
    virtual int getType() {return 0;}


};

} // namespace medyan

#endif

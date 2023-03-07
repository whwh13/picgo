
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

#ifndef MEDYAN_NeighborList_h
#define MEDYAN_NeighborList_h

#include "common.h"

namespace medyan {
//FORWARD DECLARATIONS
class Neighbor;
class DynamicNeighbor;

/// To hold an external neighbor list of general type.

/*!
 *  This class is used to hold any neighbor list. Contains a map of neighbors as well as
 *  min and max cutoffs for generation of the list. This class is abstract and must be
 *  implemented by writing functionality to add and update a neighbor.
 * 
 *  The neighbor list contains a function to reset, which uses the databases to clear
 *  and update the list.
 */

class NeighborList {

protected:
    float _rMax;  ///< max distance cutoff
    float _rMin;  ///< min distance cutoff

public:
    ///Constructor and destructor
    NeighborList(float rMax = 0.0, float rMin = 0.0) : _rMax(rMax), _rMin(rMin) {}

    ///Destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~NeighborList() noexcept {}

    /// Add neighbor
    virtual void addNeighbor(Neighbor* n) = 0;
    /// Remove a neighbor if possible
    virtual void removeNeighbor(Neighbor* n) = 0;

    /// Add a dynamic neighbor to the system.
    /// For BoundaryElementNeighborList list, this will be Bead.
    /// For CylinderNeighborList, all elements in neighbors list are dynamic.
    virtual void addDynamicNeighbor(DynamicNeighbor* n) = 0;

    /// Remove a dynamic neighbor from the system.
    /// For BoundaryElementNeighborList list, this will be Bead.
    /// For CylinderNeighborList, all elements in neighbors list are dynamic.
    virtual void removeDynamicNeighbor(DynamicNeighbor* n) = 0;

    /// Re-initialize the neighborlist
    virtual void reset() = 0;

};

} // namespace medyan

#endif

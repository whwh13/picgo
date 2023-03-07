
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_HybridNeighborList_h
#define MEDYAN_HybridNeighborList_h
#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
#include "common.h"

namespace medyan {
//FORWARD DECLARATIONS
class Neighbor;
class DynamicNeighbor;

/// To hold an external hybrid neighbor list of general type.

/*!
 *  This class is used to hold any hybrid neighbor list. Contains a map of neighbors as
 *  well as min and max cutoff vectors for generation of the list. This class is abstract
 *  and must be implemented by writing functionality to add and update a neighbor.
 *  The neighbor list contains a function to reset, which uses the databases to clear
 *  and update the list.
 */

class HybridNeighborList {

protected:
    vector<vector<float>> _rMaxsqvec; //squared maxdistance cutoff
    vector<vector<float>> _rMinsqvec;//squared mindistance cutoff
    float _smallestrMinsq = 10000000.0;
    float _largestrMaxsq = 0.0;
    float _maxcylindersize = 0.0;
    float _maxcylsize = 0.0;
    vector<vector<short>> _filamentIDvec;//filament ID pairs considered
    vector<vector<bool>> uniquestatusvec; //if the NL has to be unique, neighbor lists
    // distances
    // that are appeneded won't be compared against these.
    vector<vector<bool>>_fullstatusvec;
    vector<vector<short>> HNLIDvec; //Hybrid NL ID to track the total number of
    // neighborlists in the system
    int totaluniquefIDpairs = 0;
public:
    ///Constructor and destructor
    HybridNeighborList() {}

    ///Destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~HybridNeighborList() noexcept {}

    virtual void initializeHybridNeighborList() = 0;
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

    //create a virtual function to append min max values;
    virtual short setneighborsearchparameters(short ftype1 =0, short ftype2 =0,bool
                            uniquestatus = false, bool fullstatus = false, float rMax =0.0,
                                              float rMin = 0.0) = 0;
    static ofstream _crosscheckdumpFileNL;
};
#endif

} // namespace medyan

#endif

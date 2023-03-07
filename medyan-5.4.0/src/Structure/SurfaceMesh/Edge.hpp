#ifndef MEDYAN_Structure_SurfaceMesh_Edge_Hpp
#define MEDYAN_Structure_SurfaceMesh_Edge_Hpp

/*
 
An unordered edge contains 2 vertices.
 
*/

#include "common.h"
#include "Structure/CellList.hpp"
#include "Util/Math/Vec.hpp"
#include "Util/StableVector.hpp"

namespace medyan {

// Forward declarations.
class Compartment;
class Membrane;

class Edge {

private:

    StableVectorIndex<Membrane> parentSysIndex_ {};
    Index topoIndex_ = 0; // Index in the meshwork topology.

public:
    medyan::CellListElementUser< StableVector<Edge>::Index, Compartment* > cellElement;
    StableVectorIndex<Edge> sysIndex {};

    medyan::Vec< 3, floatingpoint > coordinate {}; // Coordinate of the mid point, updated with updateCoordiante()


    Edge() = default;
    Edge(const Edge&) = default;

    void setParentSysIndex(StableVectorIndex<Membrane> parentSysIndex) { parentSysIndex_ = parentSysIndex; }
    StableVectorIndex<Membrane> getParentSysIndex() const { return parentSysIndex_; }
    template< typename Context >
    Membrane& getParent(Context& sys) const { return sys.membranes[parentSysIndex_]; }

    void setTopoIndex(Index index) { topoIndex_ = index; }
    auto getTopoIndex() const { return topoIndex_; }

    Compartment* getCompartment() const { return cellElement.manager->getHead(cellElement); }

};

} // namespace medyan

#endif

#ifndef MEDYAN_Structure_SurfaceMesh_Triangle_Hpp
#define MEDYAN_Structure_SurfaceMesh_Triangle_Hpp

#include "common.h"
#include "Structure/CellList.hpp"
#include "Structure/DynamicNeighbor.h"
#include "Structure/SurfaceMesh/MTriangle.hpp"
#include "Util/Math/Vec.hpp"
#include "Util/StableVector.hpp"

namespace medyan {

// Forward declarations.
class Compartment;
class Membrane;

/******************************************************************************
Triangles are the only element in the meshwork that has area and act as patches
that construct the surface.

The triangle patches have geometric and mechanical properties.

The Triangle class has pointers to the vertices and edges.
******************************************************************************/
class Triangle : public DynamicNeighbor {

private:
    StableVectorIndex<Membrane> parentSysIndex_ {};
    Index _topoIndex = 0; // Index in the meshwork topology.


public:
    medyan::CellListElementUser< medyan::StableVector<Triangle>::Index, Compartment* > cellElement;
    medyan::StableVector<Triangle>::Index sysIndex {};

    // Stores triangle mechanical data
    MTriangle mTriangle;

    Vec< 3, floatingpoint > coordinate; // Coordinate of the center point, updated with updateCoordiante()

    Triangle() = default;
    Triangle(const Triangle&) = default;

    void setParentSysIndex(StableVectorIndex<Membrane> parentSysIndex) { parentSysIndex_ = parentSysIndex; }
    StableVectorIndex<Membrane> getParentSysIndex() const { return parentSysIndex_; }
    template< typename Context >
    Membrane& getParent(Context& sys) const { return sys.membranes[parentSysIndex_]; }

    void setTopoIndex(Index index) { _topoIndex = index; }
    auto getTopoIndex() const { return _topoIndex; }

    Compartment* getCompartment() const { return cellElement.manager->getHead(cellElement); }


};

} // namespace medyan

#endif

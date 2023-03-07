#ifndef MEDYAN_Membrane_hpp
#define MEDYAN_Membrane_hpp

#include <array>
#include <limits> // numeric_limits
#include <memory>
#include <stdexcept> // runtime_error
#include <vector>

#include "Structure/SurfaceMesh/MembraneMeshAttribute.hpp"
#include "Structure/SurfaceMesh/MMembrane.hpp"
#include "SysParams.h"

namespace medyan {


/******************************************************************************
Topologically, the membrane is represented by a 2d surface with 2 sides (which
means no Klein bottles are allowed!). The surface is constructed by
interconnected vertices, edges and triangles.

The Membrane class is a manager for constructing the mesh and computing the
geometry. It contains a meshwork instance that's responsible for adding and
removing vertices, edges (halfedges) and triangles to/from the SubSystem.
However, the ownership of all elements is in this Membrane class through
inheriting Composite.
******************************************************************************/
class Membrane {
public:
    using MeshAttributeType = MembraneMeshAttribute;
    using CoordinateType = typename MeshAttributeType::CoordinateType;
    using MeshType = HalfEdgeMesh< MeshAttributeType >;

private:

    MeshType  mesh_;

public:
    // A copy of the membrane setup. Should not be modified once initialized.
    MembraneSetup setup;

    MMembrane mMembrane; // Mechanical parameters


    Membrane() = default;

    /// Get vector of triangles/edges/vertices that this membrane contains.
    const auto& getMesh() const { return mesh_; }
    auto&       getMesh()       { return mesh_; }

    // Get membrane setup.
    const auto& getSetup() const { return setup; }


    // Print self information
    void printSelf() const {
        using namespace std;

        cout << endl;

        cout << "Membrane type = " << getSetup().type << endl;
        cout << "Number of vertices, edges, half edges, triangles, borders =\n  "
            << mesh_.numVertices() << ' ' << mesh_.numEdges() << ' ' << mesh_.numHalfEdges() << ' '
            << mesh_.numTriangles() << ' ' << mesh_.numBorders() << endl;

        cout << endl;
    }


    /**************************************************************************
    Topological
    **************************************************************************/
    bool isClosed() const { return mesh_.isClosed(); }


};

} // namespace medyan

#endif

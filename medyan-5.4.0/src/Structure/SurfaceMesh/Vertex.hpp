#ifndef MEDYAN_Structure_SurfaceMesh_Vertex_hpp
#define MEDYAN_Structure_SurfaceMesh_Vertex_hpp

#include <memory> // unique_ptr
#include <vector>

#include <Eigen/Dense>

#include "Chemistry/ReactionDy.hpp"
#include "Chemistry/SpeciesContainer.h"
#include "MathFunctions.h"
#include "Util/StableVector.hpp"

namespace medyan {
// Forward declarations.
class Compartment;
struct Membrane;


struct MVertex {
    // Equilibrium area.
    double eqArea = 1;
};

// CVertex represents a cell around a vertex related to chemistry, and owns
//   - all the diffusing species in the cell
//   - all reactions with only diffusing species in this cell
//
// Note:
//   - diffusion reactions between the cells will be stored in CHalfEdge
struct CVertex {
    using ReactionContainer = std::vector< std::unique_ptr< ReactionDy > >;

    SpeciesPtrContainerVector species;
    ReactionContainer         reactions;
    // Note: adsorption reactions might be optimized to be compartment-based (one reaction plus callback), which will no longer use this container.
    ReactionContainer         adsorptionReactions;
    ReactionContainer         desorptionReactions;

    // The accessible area for protein adsorption.
    // Updated by adsorption/desorption reaction callbacks and their rate update.
    // This value can be negtive occasionally, where it should be treated as zero.
    FP                        accessibleArea = 0;
};

struct Vertex {
    
    using CoordinateType = Vec< 3, floatingpoint >;

    StableVectorIndex<Vertex> sysIndex {};

    // The index when looped in membrane-vertex sequence.
    // Unlike loopIndex in other structures, this one is not directly related to the minimization DOF serialization.
    // Updated during DOF serialization.
    Index                     loopIndex = 0;

    // The parent membrane system index.
    StableVectorIndex<Membrane> parentSysIndex {};
    // The index corresponding to the vertex in the membrane mesh.
    Index                       topoIndex = 0;

    // The counter for other attachments referencing this vertex.
    // During vertex finalization, if this counter is not zero, an error should be thrown.
    Size                        attachmentRefCount = 0;

    // Pinned vertex does not move during mechanical equilibration.
    // Updated during DOF serialization.
    bool                      pinned = false;

    CoordinateType coord {};
    CoordinateType force {};

    MVertex mVertex; // vertex mechanical information
    CVertex cVertex; // vertex chemical information

    // Vertex as an element in a compartment.
    medyan::CellListElementUser< medyan::StableVector<Vertex>::Index, Compartment* > cellElement;


    // Accessors and mutators.
    void setParentSysIndex(StableVectorIndex<Membrane> parentSysIndex) { this->parentSysIndex = parentSysIndex; }
    auto getParentSysIndex() const { return parentSysIndex; }
    void setTopoIndex(Index index) { topoIndex = index; }
    auto getTopoIndex() const { return topoIndex; }

    auto getAttachmentRefCount() const { return attachmentRefCount; }
};

// Structure used in meshless membrane representation.
struct MeshlessSpinVertex {
    using CoordinateType = Eigen::Vector3d;

    // Managed indices.
    //----------------------------------

    // Stable index. This will not change during its lifetime.
    // Can be used anywhere.
    // Should never be updated.
    StableVectorIndex<MeshlessSpinVertex>           sysIndex {};
    // Looping index. This may change, but will be contiguous for all such vertices.
    // Used as the sequence in mechanical vectorization.
    // Updated during DOF serialization.
    Index                                           loopIndex = 0;
    // Index in vertex-vertex cell list for energy minimization purpose.
    // Used when updating cell list for this vertex is necessary.
    // Updated during cell list building.
    Index                                           meshlessSpinVertexCellListIndex = 0;

    CoordinateType coord {};
    // Unit vector representing orientation.
    // Note:
    // - θ ∈ [0, π]
    // - ϕ ∈ [0, 2π)
    // - (x, y, z) = (cos(ϕ)sin(θ), sin(ϕ)sin(θ), cos(θ))
    // - Values can be made out of range by minimization algorithms, which is okay.
    double         theta = 0;
    double         phi = 0;
};


constexpr auto initializeVertex = [](
    auto&&, auto&, Vertex& v, StableVectorIndex<Vertex> index,
    const Vertex::CoordinateType& coord,
    StableVectorIndex<Membrane>   parentSysIndex,
    Index                         topoIndex
) {
    v.coord = coord;
    v.setParentSysIndex(parentSysIndex);
    v.setTopoIndex(topoIndex);
    // Register the stable index inside vertex.
    v.sysIndex = index;
};
constexpr auto finalizeVertex = [](auto&&, auto&, Vertex& v, StableVectorIndex<Vertex> index) {
    if(v.attachmentRefCount != 0) {
        log::error("Vertex at index {} has non-zero attachment reference count {}.", index.value, v.attachmentRefCount);
        throw std::runtime_error("Vertex has non-zero attachment reference count.");
    }
};

constexpr auto initializeMeshlessSpinVertex = [](auto&&, auto&, MeshlessSpinVertex& v, auto index) {
    // Register the stable index inside vertex.
    v.sysIndex = index;
};
constexpr auto finalizeMeshlessSpinVertex = [](auto&&, auto&, MeshlessSpinVertex&, auto) {};

} // namespace medyan

#endif

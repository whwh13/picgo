#ifndef MEDYAN_Structure_SurfaceMesh_FixedVertexAttachment_hpp
#define MEDYAN_Structure_SurfaceMesh_FixedVertexAttachment_hpp

#include "common.h"
#include "Util/StableVector.hpp"
#include "Util/Math/Vec.hpp"

namespace medyan {

// Forward declarations.
struct Vertex;

struct FixedVertexAttachment {
    // The vertex used.
    StableVectorIndex<Vertex> vertexSysIndex {};

    Vec<3, FP> coord {};

    // Mechanical properties.
    FP kStretch = 0;
};

constexpr auto initializeFixedVertexAttachment = [](
    auto&&, auto& sys, FixedVertexAttachment& attachment, auto index,
    const FixedVertexAttachment& refAttachment
) {
    attachment = refAttachment;

    // Increase ref counter in the vertex.
    ++ sys.vertices[attachment.vertexSysIndex].attachmentRefCount;
};
constexpr auto finalizeFixedVertexAttachment = [](auto&&, auto& sys, FixedVertexAttachment& attachment, auto) {
    // Decrease ref counter in the vertex.
    -- sys.vertices[attachment.vertexSysIndex].attachmentRefCount;
};

} // namespace medyan

#endif

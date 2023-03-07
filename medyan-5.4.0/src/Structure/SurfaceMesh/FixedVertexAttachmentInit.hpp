#ifndef MEDYAN_Structure_SurfaceMesh_FixedVertexAttachmentInit_hpp
#define MEDYAN_Structure_SurfaceMesh_FixedVertexAttachmentInit_hpp

#include "Structure/SubSystem.h"
#include "Structure/SurfaceMesh/FixedVertexAttachment.hpp"
#include "SysParams.h"

namespace medyan {

std::optional<FixedVertexAttachment> createFixedVertexAttachment(SubSystem& sys, const FixedVertexAttachmentInitSearch& init) {
    std::optional<FixedVertexAttachment> res;

    StableVectorIndex<Vertex> minDistVertexSysIndex { -1 };
    FP minDist2 = inffp;

    // Find the closest vertex.
    for(auto it = sys.vertices.begin(); it != sys.vertices.end(); ++it) {
        const auto dist2 = distance2(it->coord, init.coord);
        if(dist2 < init.range * init.range && dist2 < minDist2) {
            minDistVertexSysIndex = sys.vertices.indexat(it);
            minDist2 = dist2;
        }
    }

    if(minDistVertexSysIndex.value >= 0) {
        // Found.
        res = FixedVertexAttachment {
            minDistVertexSysIndex,
            init.coord,
            init.kStretch,
        };
    }

    return res;
}


// Returns number of attachments added.
template< typename ContextFunc >
Size createFixedVertexAttachments(SubSystem& sys, const FixedVertexAttachmentInitSearchMultiple& inits) {

    Size count = 0;

    for(int i = 0; i < inits.attachments.size(); ++i) {
        const auto& init = inits.attachments[i];
        if(auto attachment = createFixedVertexAttachment(sys, init)) {
            ContextFunc{}.template emplaceTrackable<FixedVertexAttachment>(sys, *attachment);
            ++count;
        }
        else if(inits.throwOnNotFound) {
            log::error("Could not find a vertex to attach within the specified range.");
            throw std::runtime_error("Could not find a vertex within the specified range.");
        }
        else {
            log::warn("Could not find a vertex to attach within the specified range.");
        }
    }

    return count;
}

} // namespace medyan

#endif

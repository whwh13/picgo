#ifndef MEDYAN_Structure_SurfaceMesh_HalfEdge_hpp
#define MEDYAN_Structure_SurfaceMesh_HalfEdge_hpp

#include "Chemistry/ReactionDy.hpp"

namespace medyan {

struct CHalfEdge {
    struct EachReaction {
        std::unique_ptr< ReactionDy > reaction;

        // Sentinel: -1. The index of diffusing species in the membrane.
        int membraneDiffusingSpeciesIndex = -1;
    };

    using ReactionContainer = std::vector< EachReaction >;

    ReactionContainer diffusionReactions;
};

struct HalfEdge {
    CHalfEdge cHalfEdge;
};

} // namespace medyan

#endif

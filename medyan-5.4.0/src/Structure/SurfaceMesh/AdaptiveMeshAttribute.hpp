#ifndef MEDYAN_AdaptiveMeshAttribute_hpp
#define MEDYAN_AdaptiveMeshAttribute_hpp

#include "MathFunctions.h"

namespace medyan {
// Additional attributes needed for meshwork
struct AdaptiveMeshAttribute {
    struct VertexAttribute {
        double size;
        double sizeAux; // Used in diffusing
        medyan::Vec3 unitNormal;
    };
    struct HalfEdgeAttribute {
    };
    struct EdgeAttribute {
        double eqLength;
    };
    struct TriangleAttribute {
    };
};

} // namespace medyan

#endif

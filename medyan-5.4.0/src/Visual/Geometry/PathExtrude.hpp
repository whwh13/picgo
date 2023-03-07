#ifndef MEDYAN_Visual_Geometry_PathExtrude_Hpp
#define MEDYAN_Visual_Geometry_PathExtrude_Hpp

#include <array>
#include <cstdint> // uint_fast8_t
#include <tuple>
#include <type_traits>
#include <variant> // monostate
#include <vector>

#include <Eigen/Core>

#include "common.h"
#include "Util/Math/Vec.hpp"

namespace medyan::visual {

// This function transforms a path to a mesh of tubes.
// The return mesh is represented as GL_TRIANGLES compatible vertex coord list and index list.
//
// Parameters
//   - coords:  the container where the coordinates of beads can be found.
//   - indices: the bead indices on the path in the coords container.
//   - radius:  the radius of the tube.
//   - sides:   the number of sides of the tube.
template<
    // Manually specified types.
    typename Float,               // Resulting floating point type.
    typename AttributeType,       // Type of attribute in the resulting array. It can be void, in which case the attribute is not used.
    // Deduced.
    typename CoordContainer,      // Container of coordinates.
    typename FuncCoordData,       // coord_container, index -> ptrcoord.
    typename AttributeContainer,  // Container of attributes. May be the same as CoordContainer.
    typename FuncGetAttribute,    // attrib_container, index -> attribute.
    typename IndexContainer,      // Container of indices.
    std::enable_if_t< std::is_floating_point_v<Float> >* = nullptr
>
inline auto pathExtrudeGenerate(
    const CoordContainer& coords,
    FuncCoordData&&       funcCoordData,
    const AttributeContainer& attribs,
    FuncGetAttribute&&    funcGetAttribute,
    const IndexContainer& indices,
    Float                 radius,
    int                   sides
) {
    using CoordType = Vec<3, Float>;
    constexpr bool useAttrib = !std::is_same_v<AttributeType, void>;
    using VecAttrib = std::conditional_t<
        useAttrib,
        std::vector< AttributeType >,
        std::monostate // Dummy type.
    >;


    constexpr CoordType a0 { 1.0, 0.0, 0.0 };
    constexpr CoordType a1 { 0.0, 1.0, 0.0 }; // Unused unless first segment is parallel to a0

    const Size numVertices = indices.size();
    const Size numTubeVertices = sides * numVertices;
    const Size numTriangles = (numVertices - 1) * 2 * sides;

    // Results
    std::vector< CoordType > vertices;      vertices.reserve(numTubeVertices);
    std::vector< CoordType > vertexNormals; vertexNormals.reserve(numTubeVertices);
    VecAttrib attributes;                   if constexpr(useAttrib) { attributes.reserve(numTubeVertices); }
    std::vector< std::array< int, 3 > > triInd(numTriangles);

    if(indices.size() < 2) {
        if constexpr(useAttrib) {
            return std::make_tuple(std::move(vertices), std::move(vertexNormals), std::move(attributes), std::move(triInd));
        } else {
            return std::make_tuple(std::move(vertices), std::move(vertexNormals), std::move(triInd));
        }
    }

    const auto getCoords = [&](Index index) {
        return makeRefVec<3>(funcCoordData(coords, index));
    };

    // Locate first circle
    CoordType seg ( normalizedVector(getCoords(indices[1]) - getCoords(indices[0])) );
    auto n0 = cross(seg, a0);
    if(magnitude2(n0) == (Float)0) n0 = cross(seg, a1);
    normalize(n0);
    auto n1 = normalizedVector(cross(n0, seg));
    auto segn = seg;

    for(Index j = 0; j < sides; ++j) {
        const Float a = j * 2 * M_PI / sides;
        const Float cosa = std::cos(a);
        const Float sina = std::sin(a);
        const auto point = CoordType(getCoords(indices[0])) + n0 * (radius * cosa) + n1 * (radius * sina);
        const auto un    = n0 * cosa + n1 * sina;
        vertices.push_back(point);
        vertexNormals.push_back(un);
        if constexpr(useAttrib) { attributes.push_back(funcGetAttribute(attribs, indices[0])); }
    }

    // Propagate circles
    for(Index i = 1; i < numVertices; ++i) {
        segn = normalizedVector(i == numVertices - 1 ? getCoords(indices[i]) - getCoords(indices[i-1]) : getCoords(indices[i+1]) - getCoords(indices[i]));
        const auto t = normalizedVector(seg + segn);
        const auto dot_seg_t = dot(seg, t);

        for(Index j = 0; j < sides; ++j) {
            // Solve for p_new given:
            //   dot(p_new - coords[indices[i]], t) = 0
            //   p_new = x * seg + pp
            const auto pp = vertices[(i-1) * sides + j];
            const Float x = dot(getCoords(indices[i]) - pp, t) / dot_seg_t;
            vertices.push_back(x * seg + pp);
            vertexNormals.push_back(normalizedVector(vertices.back() - CoordType(getCoords(indices[i]))));
            if constexpr(useAttrib) { attributes.push_back(funcGetAttribute(attribs, indices[i])); }
        }

        seg = segn;
    }

    // Make indices (GL_TRIANGLES)
    for(Index i = 0; i < numVertices - 1; ++i) {
        for(Index j = 0; j < sides; ++j) {
            // First triangle
            triInd[i * (2 * sides) + 2 * j][0] = (i    ) * sides + j;
            triInd[i * (2 * sides) + 2 * j][1] = (i + 1) * sides + j;
            triInd[i * (2 * sides) + 2 * j][2] = (i    ) * sides + (j + 1) % sides;
            // Second triangle
            triInd[i * (2 * sides) + 2 * j + 1][0] = (i    ) * sides + (j + 1) % sides;
            triInd[i * (2 * sides) + 2 * j + 1][1] = (i + 1) * sides + j;
            triInd[i * (2 * sides) + 2 * j + 1][2] = (i + 1) * sides + (j + 1) % sides;
        }
    }

    if constexpr(useAttrib) {
        return std::make_tuple(std::move(vertices), std::move(vertexNormals), std::move(attributes), std::move(triInd));
    } else {
        return std::make_tuple(std::move(vertices), std::move(vertexNormals), std::move(triInd));
    }

}

// This special case
// - Uses vector of Vec's.
// - Does not use attributes.
template<
    typename Float,
    typename CoordContainer,
    typename IndexContainer,
    std::enable_if_t<
        std::is_floating_point_v<Float> &&
        !std::is_base_of_v< Eigen::MatrixBase<CoordContainer>, CoordContainer > &&
        IsVecLike<typename CoordContainer::value_type>::value
    >* = nullptr
>
inline auto pathExtrudeGenerate(
    const CoordContainer& coords,
    const IndexContainer& indices,
    Float                 radius,
    int                   sides
) {
    return pathExtrudeGenerate<Float, void>(
        coords,
        [](const CoordContainer& coords, Index index) { return coords[index].data(); },
        coords,
        [](const CoordContainer& coords, Index index) {},
        indices,
        radius,
        sides
    );
}

// This special case
// - Uses Eigen matrices.
// - The attribute container is the coordinate container.
// - The attribute type is Float.
template<
    typename Float,
    typename CoordContainer,
    typename FuncGetAttribute,
    typename IndexContainer,
    std::enable_if_t<
        std::is_floating_point_v<Float> &&
        std::is_base_of_v< Eigen::MatrixBase<CoordContainer>, CoordContainer >
    >* = nullptr
>
inline auto pathExtrudeGenerateWithAttrib(
    const CoordContainer& coords,
    FuncGetAttribute&&    funcGetAttribute,
    const IndexContainer& indices,
    Float                 radius,
    int                   sides
) {
    return pathExtrudeGenerate<Float, Float>(
        coords,
        [](const CoordContainer& coords, Index index) { return &coords(0, index); },
        coords,
        std::forward<FuncGetAttribute>(funcGetAttribute),
        indices,
        radius,
        sides
    );
}


// Provide an estimate on the number of triangles returned by "generate" function.
inline int pathExtrudeEstimateNumTriangles(int numVertices, int sides) {
    if(numVertices < 2) return 0;

    const int numSegments = numVertices - 1;
    return numSegments * sides * 2;
}

} // namespace medyan::visual

#endif

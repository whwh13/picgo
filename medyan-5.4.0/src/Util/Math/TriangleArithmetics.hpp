#ifndef MEDYAN_Util_Math_TriangleArithmetics_Hpp
#define MEDYAN_Util_Math_TriangleArithmetics_Hpp

#include <stdexcept> // logic_error

#include "Util/Math/Vec.hpp"

namespace medyan {

//-----------------------------------------------------------------------------
// Calculate the min distance between a point and any point on a triangle
//-----------------------------------------------------------------------------
template<
    typename VT1, typename VT2,
    std::enable_if_t< VT1::vec_size == 3 && VT2::vec_size == 3 >* = nullptr
> inline auto trianglePointDistance(
    const VT1& v0, const VT1& v1, const VT1& v2,
    const VT2& p
) {
    using Float = std::common_type_t< typename VT1::float_type, typename VT2::float_type >;

    Float d;

    //-------------------------------------------------------------------------
    // Calculate the barycentric coordinate of the projection point p'
    //
    // See Heidrich 2005, Computing the Barycentric Coordinates of a Projected
    // Point.
    //-------------------------------------------------------------------------
    const auto r01 = v1 - v0;
    const auto r02 = v2 - v0;
    const auto r0p = p - v0;
    const auto cp = cross(r01, r02);
    const auto oneOver4AreaSquared = (Float)1.0 / magnitude2(cp);

    const auto b1 = dot(cross(r0p, r02), cp) * oneOver4AreaSquared;
    const auto b2 = dot(cross(r01, r0p), cp) * oneOver4AreaSquared;
    const auto b0 = (Float)1.0 - b1 - b2;

    if(b0 >= 0 && b1 >= 0 && b2 >= 0) {
        // p' is inside the triangle
        d = std::abs(dot(normalizedVector(cp), r0p));

    } else {
        // p' is outside the triangle
        const Vec< 3, Float > r2 {
            distance2(v1, v2),
            distance2(v2, v0),
            distance2(v0, v1)
        };
        const auto r1p = p - v1;
        const auto r2p = p - v2;
        const auto r12 = v2 - v1;
        const auto dot_1p_12 = dot(r1p, r12);
        const auto dot_2p_20 = -dot(r2p, r02);
        const auto dot_0p_01 = dot(r0p, r01);

        if(b0 < 0 && dot_1p_12 >= 0 && dot_1p_12 <= r2[0]) {
            // On edge 12
            d = magnitude(cross(r1p, r12)) / std::sqrt(r2[0]);
        } else if(b1 < 0 && dot_2p_20 >= 0 && dot_2p_20 <= r2[1]) {
            // On edge 20
            d = magnitude(cross(r2p, r02)) / std::sqrt(r2[1]);
        } else if(b2 < 0 && dot_0p_01 >= 0 && dot_0p_01 <= r2[2]) {
            // On edge 01
            d = magnitude(cross(r0p, r01)) / std::sqrt(r2[2]);
        } else if(dot_0p_01 < 0 && dot_2p_20 > r2[1]) {
            // On vertex 0
            d = distance(v0, p);
        } else if(dot_1p_12 < 0 && dot_0p_01 > r2[2]) {
            // On vertex 1
            d = distance(v1, p);
        } else if(dot_2p_20 < 0 && dot_1p_12 > r2[0]) {
            // On vertex 2
            d = distance(v2, p);
        } else {
            // The program should never come here
            throw std::logic_error("Unknown case of point projection on the plane of triangle.");
        }
    }

    return d;
}

} // namespace medyan

#endif

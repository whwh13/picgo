#ifndef MEDYAN_Util_Math_RayTriangleIntersect_Hpp
#define MEDYAN_Util_Math_RayTriangleIntersect_Hpp

#include <type_traits> // common_type, enable_if_t

#include "MathFunctions.h"
#include "utility.h" // areEqual

namespace medyan::ray_tracing {

// Given point O and ray vector D, and a triangle defined by points
// (A, B, C), test whether the ray intersects the triangle.
// Let the intersection between line (O, D) and plane (A, B, C) be P,
// the P could be represented by P = A + u(B-A) + v(C-A), as well as
// P = O + tD.

/// Stores intersect result: [Ray intersects triangle, t, u, v]
template< typename Float = double >
struct RayTriangleIntersectResult {
    bool intersect;
    Float t, u, v;
};

template <typename Derived>
struct RayTriangleIntersectBase {

    template<
        typename VTO, typename VTD,
        typename VTA, typename VTB, typename VTC,
        std::enable_if_t<
            VTO::vec_size == 3 && VTD::vec_size == 3 &&
            VTA::vec_size == 3 && VTB::vec_size == 3 && VTC::vec_size == 3
        >* = nullptr
    > auto operator()(
        const VTO& o, const VTD& d,
        const VTA& a, const VTB& b, const VTC& c
    ) const {
        return static_cast<const Derived*>(this)->intersect(o, d, a, b, c);
    }
};

// A fast ray-triangle intersection algorithm by Moller and Trumbore
// @inproceedings{moller2005fast,
//   title={Fast, minimum storage ray/triangle intersection},
//   author={M{\"o}ller, Tomas and Trumbore, Ben},
//   booktitle={ACM SIGGRAPH 2005 Courses},
//   pages={7},
//   year={2005},
//   organization={ACM}
// }
template< bool stopWhenNotIntersect = true >
struct MollerTrumboreIntersect
    : RayTriangleIntersectBase< MollerTrumboreIntersect< stopWhenNotIntersect > >
{

    // If stopWhenNotIntersect is true and result.intersect is false,
    // the rest of the result is undefined.
    template<
        typename VTO, typename VTD,
        typename VTA, typename VTB, typename VTC,
        std::enable_if_t<
            VTO::vec_size == 3 && VTD::vec_size == 3 &&
            VTA::vec_size == 3 && VTB::vec_size == 3 && VTC::vec_size == 3
        >* = nullptr
    > auto intersect(
        const VTO& o, const VTD& d,
        const VTA& a, const VTB& b, const VTC& c
    ) const {

        using namespace mathfunc;
        using Float = std::common_type_t<
            typename VTO::float_type, typename VTD::float_type,
            typename VTA::float_type, typename VTB::float_type, typename VTC::float_type
        >;

        RayTriangleIntersectResult< Float > res;

        const auto rab = b - a;
        const auto rac = c - a;
        const auto d_x_rac = cross(d, rac);
        const auto det = dot(d_x_rac, rab);

        if(areEqual(det, 0.0)) {
            res.intersect = false;
            if(stopWhenNotIntersect) return res;
        }

        const auto invDet = static_cast<Float>(1.0) / det;

        const auto rao = o - a;
        res.u = dot(d_x_rac, rao) * invDet;
        if(res.u < 0 || res.u > 1) {
            res.intersect = false;
            if(stopWhenNotIntersect) return res;
        }

        const auto rao_x_rab = cross(rao, rab);
        res.v = dot(rao_x_rab, d) * invDet;
        if(res.v < 0 || res.u + res.v > 1) {
            res.intersect = false;
            if(stopWhenNotIntersect) return res;
        }

        res.intersect = true;
        res.t = dot(rao_x_rab, rac) * invDet;
        return res;

    }
};

} // namespace medyan::ray_tracing

#endif

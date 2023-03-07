#ifndef MEDYAN_Structure_SurfaceMesh_MeshTriangleQuality_hpp
#define MEDYAN_Structure_SurfaceMesh_MeshTriangleQuality_hpp

#include <limits>
#include <type_traits> // enable_if, is_floating_point

#include "MathFunctions.h"

namespace medyan {
enum class TriangleQualityCriteria {
    radiusRatio     // Circumradius / (2 * Inradius), range [1, inf)
};
template< TriangleQualityCriteria > struct TriangleQuality;
template<> struct TriangleQuality< TriangleQualityCriteria::radiusRatio > {
    static constexpr double best  = 1.0;
    static constexpr double worst = std::numeric_limits<double>::infinity();

    static constexpr bool better   (double q1, double q2) { return q1 < q2; }
    static constexpr auto betterOne(double q1, double q2) { return better(q1, q2) ? q1 : q2; }
    static constexpr bool worse    (double q1, double q2) { return q1 > q2; }
    static constexpr auto worseOne (double q1, double q2) { return worse(q1, q2) ? q1 : q2; }
    static constexpr auto improvement(double q0, double q1) { return q0 / q1; }

    template<
        typename VT0, typename VT1, typename VT2,
        std::enable_if_t<
            (VT0::vec_size > 0) &&
            VT0::vec_size == VT1::vec_size &&
            VT0::vec_size == VT2::vec_size
        >* = nullptr
    > auto operator()(const VT0& v0, const VT1& v1, const VT2& v2) const {
        using namespace mathfunc;
        const auto d0 = distance(v1, v2);
        const auto d1 = distance(v2, v0);
        const auto d2 = distance(v0, v1);
        return operator()(d0, d1, d2);
    }

    template<
        typename F0, typename F1, typename F2,
        std::enable_if_t<
            std::is_floating_point_v< F0 > &&
            std::is_floating_point_v< F1 > &&
            std::is_floating_point_v< F2 >
        >* = nullptr
    > auto operator()(F0 d0, F1 d1, F2 d2) const {
        const auto p = 0.5 * (d0 + d1 + d2);
        // Note that the abs is needed to avoid extreme cases where the result is negative.
        return std::abs(d0 * d1 * d2 / (8 * (p - d0) * (p - d1) * (p - d2)));
    }
};

} // namespace medyan

#endif

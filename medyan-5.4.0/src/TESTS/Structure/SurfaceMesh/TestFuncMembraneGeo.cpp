#include <catch2/catch.hpp>

#include "Structure/SurfaceMesh/FuncMembraneGeo.hpp"

namespace medyan {

TEST_CASE("Membrane geometry functions", "[Mesh]") {

    SECTION("Triangle area.") {
        const Vec3d v1 { 0, 0, 0 };
        const Vec3d v2 { 1, 0, 0 };
        const Vec3d v3 { 0, 1, 0 };

        {
            const auto area = medyan::area(v1, v2, v3);
            CHECK(area == Approx(0.5));
        }
        {
            const auto [area, da1, da2, da3] = areaAndDerivative(v1, v2, v3);
            CHECK(area == Approx(0.5));
            CHECK(da1[0] == Approx(-0.5));
            CHECK(da1[1] == Approx(-0.5));
            CHECK(da1[2] == 0.0);
            CHECK(da2[0] == Approx(0.5));
            CHECK(da2[1] == 0.0);
            CHECK(da2[2] == 0.0);
            CHECK(da3[0] == 0.0);
            CHECK(da3[1] == Approx(0.5));
            CHECK(da3[2] == 0.0);
        }
    }
}

} // namespace medyan

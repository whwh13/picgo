#include <catch2/catch.hpp>

#include "Structure/SurfaceMesh/SurfaceMeshGenerator.hpp"

namespace medyan {

TEST_CASE("Surface mesh generator - Marching Tetrahedra", "[Mesh][MeshGen]") {
    using namespace medyan::mesh_gen;

    SECTION("Full sphere") {
        MarchingTetrahedraGenerator<> gen(
            // Box size.
            1.0,
            // Box origin.
            { -10.0, -10.0, -10.0 },
            // Num boxes.
            { 20, 20, 20 }
        );
        auto result = gen([](const auto& p) {
            return p[0] * p[0] + p[1] * p[1] + p[2] * p[2] - 8 * 8;
        });

        // Check the number of elements.
        const auto nv = result.vertexCoordinateList.size();
        const auto nt = result.triangleList.size();
        REQUIRE(nt % 2 == 0);
        const auto neExpected = nt / 2 * 3;
        CHECK(nv + nt - neExpected == 2);
    }

    SECTION("Donut") {
        MarchingTetrahedraGenerator<> gen(
            1.0,
            { -10.0, -10.0, -10.0 },
            { 20, 20, 20 }
        );
        auto result = gen([](const auto& p) {
            constexpr double r = 5;
            constexpr double a2 = 2.5 * 2.5;
            const auto xy = std::sqrt(p[0] * p[0] + p[1] * p[1]);
            return (xy - r) * (xy - r) + p[2] * p[2] - a2;
        });

        // Check the number of elements.
        const auto nv = result.vertexCoordinateList.size();
        const auto nt = result.triangleList.size();
        REQUIRE(nt % 2 == 0);
        const auto neExpected = nt / 2 * 3;
        CHECK(nv + nt - neExpected == 0);
    }
}

} // namespace medyan

#include <catch2/catch.hpp>

#include "Visual/Geometry/PathExtrude.hpp"

namespace medyan {

TEST_CASE("Geometry path extrusion test", "[Visual]") {
    using namespace medyan::visual;

    // Triangle path.
    {
        constexpr std::array< Vec<3, float>, 3 > coords {{
            { 0.0f, 0.0f, 0.0f },
            { 1.0f, 0.0f, 0.0f },
            { 0.0f, 1.0f, 0.0f },
        }};
        constexpr std::array< int, 4 > indices { 0, 1, 2, 0 };
        constexpr double radius = 0.5;
        constexpr int sides = 3;

        auto [vertices, vertexNormals, triInd] = pathExtrudeGenerate<double>(
            coords,
            indices,
            radius,
            sides
        );

        CHECK(vertices.size() == sides * indices.size());
        CHECK(vertexNormals.size() == vertices.size());
        CHECK(triInd.size() == pathExtrudeEstimateNumTriangles(indices.size(), sides));
    }

    // Segment path.
    {
        // With more vertex attributes.
        Eigen::MatrixXd coords(4, 2);
        coords <<
            0.0, 0.0,
            1.0, 0.0,
            0.0, 1.0,
            1.0, 1.0;
        const std::vector< long > indices { 1, 0 };
        constexpr double radius = 1.0f;
        constexpr int sides = 6;

        auto [vertices, vertexNormals, attributes, triInd] = pathExtrudeGenerateWithAttrib<float>(
            coords,
            [](const auto& attributes, Index index) { return 0; },
            indices,
            radius,
            sides
        );

        CHECK(vertices.size() == sides * indices.size());
        CHECK(vertexNormals.size() == vertices.size());
        CHECK(attributes.size() == vertices.size());
        CHECK(triInd.size() == pathExtrudeEstimateNumTriangles(indices.size(), sides));
    }
}

} // namespace medyan

#include "catch2/catch.hpp"

#include "Util/Math/TriangleArithmetics.hpp"

TEST_CASE("Triangle point distance test", "[TrianglePointDistance]") {
    using namespace medyan;

    Vec3d ta0 { 0.0, 0.0, 0.0 };
    Vec3d tb0 { 2.0, 0.0, 0.0 };
    Vec3d tc0 { 0.0, 2.0, 0.0 };

    Vec3d p0 { 0.5, 0.5, 2.0 }; // inside, 2
    Vec3d p1 { -1.0, -2.0, 2.0 }; // a, 3
    Vec3d p2 { 3.0, -2.0, 2.0 }; // b, 3
    Vec3d p3 { -2.0, 3.0, 2.0 }; // c, 3
    Vec3d p4 { 3.0, 3.0, 1.0 }; // bc, 3
    Vec3d p5 { -3.0, 1.0, 4.0 }; // ca, 5
    Vec3d p6 { 1.0, -3.0, 4.0 }; // ab, 5

    SECTION("Nearest distance") {
        const auto getDistance = [&](const Vec3d p) { return trianglePointDistance(ta0, tb0, tc0, p); };

        const auto res0 = getDistance(p0); CHECK(res0 == Approx(2.0));
        const auto res1 = getDistance(p1); CHECK(res1 == Approx(3.0));
        const auto res2 = getDistance(p2); CHECK(res2 == Approx(3.0));
        const auto res3 = getDistance(p3); CHECK(res3 == Approx(3.0));
        const auto res4 = getDistance(p4); CHECK(res4 == Approx(3.0));
        const auto res5 = getDistance(p5); CHECK(res5 == Approx(5.0));
        const auto res6 = getDistance(p6); CHECK(res6 == Approx(5.0));
    }
}

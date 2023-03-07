#include "catch2/catch.hpp"

#include "Util/Math/RayTriangleIntersect.hpp"

using namespace medyan;

TEST_CASE("Ray triangle intersection test", "[Ray Triangle Intersect]") {
    using namespace ray_tracing;

    Vec3 ta0 { 0.0, 0.0, 0.0 };
    Vec3 tb0 { 2.0, 0.0, 0.0 };
    Vec3 tc0 { 0.0, 2.0, 0.0 };

    Vec3 o0  { 0.0, 0.0, 2.0 };
    Vec3 d0  { 1.0, 1.0, 0.0 }; // parallel
    Vec3 d1  { 1.0, -1.0, -2.0 }; // outside
    Vec3 d2  { -1.0, 1.0, -2.0 }; // outside
    Vec3 d3  { 2.0, 2.0, -2.0 }; // outside
    Vec3 d4  { 0.5, 0.5, -2.0 }; // inside
    Vec3 d5  { -0.5, -0.5, 2.0 }; // inside, negative

    SECTION("Moller Trumbore method") {
        MollerTrumboreIntersect<> f{};

        // parallel case
        const auto res0 = f(o0, d0, ta0, tb0, tc0);
        CHECK_FALSE(res0.intersect);

        // outside cases
        const auto res1 = f(o0, d1, ta0, tb0, tc0);
        const auto res2 = f(o0, d2, ta0, tb0, tc0);
        const auto res3 = f(o0, d3, ta0, tb0, tc0);
        CHECK_FALSE(res1.intersect);
        CHECK_FALSE(res2.intersect);
        CHECK_FALSE(res3.intersect);

        // intersect cases
        const auto res4 = f(o0, d4, ta0, tb0, tc0);
        CHECK(res4.intersect);
        CHECK(res4.t == Approx(1.0));
        CHECK(res4.u == Approx(0.25));
        CHECK(res4.v == Approx(0.25));

        const auto res5 = f(o0, d5, ta0, tb0, tc0);
        CHECK(res5.intersect);
        CHECK(res5.t == Approx(-1.0));
        CHECK(res5.u == Approx(0.25));
        CHECK(res5.v == Approx(0.25));
    }
}

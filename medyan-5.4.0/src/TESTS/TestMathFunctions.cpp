#include <catch2/catch.hpp>

#include "MathFunctions.h"

TEST_CASE("Math functions", "[Math]") {
    using namespace std;
    using namespace medyan;
    using namespace mathfunc;

    SECTION("FENE potential") {
        {
            INFO("FENE potential energy values");
            CHECK(fene(FeneDistSq{ 0.0 }, 5.0, 100.0f) == 0.0);
            CHECK(fene(FeneDistSq{ 2500.0 }, 5.0, 100.0f) > 0);
            CHECK(isfinite(fene(FeneDistSq{ 2500.0 }, 5.0, 100.0f)));
            CHECK(fene(FeneDistSq{ 10000.0 }, 5.0, 100.0f) == std::numeric_limits<double>::infinity());
            CHECK(fene(FeneDistSq{ 40000.0 }, 5.0, 100.0f) == std::numeric_limits<double>::infinity());
        }
        {
            INFO("FENE gradient values");
            CHECK(dFeneCoeff(FeneDistSq{ 0.0 }, 1.0, 200.0) == 1.0);
            CHECK(dFeneCoeff(FeneDistSq{ 10000.0 }, 1.0, 200.0f) > 0);
            CHECK(isfinite(fene(FeneDistSq{ 10000.0 }, 1.0, 200.0f)));
            CHECK(dFeneCoeff(FeneDistSq{ 40000.0 }, 1.0, 200.0) == std::numeric_limits<double>::infinity());
            CHECK(dFeneCoeff(FeneDistSq{ 160000.0 }, 1.0, 200.0) == std::numeric_limits<double>::infinity());
        }
        {
            INFO("FENE value/gradient consistency");
            {
                INFO("1D case.");
                const double d = 100.0;
                const double d0 = 50.0;
                const double dd = 1e-6;

                const double e1 = fene(FeneDistSq{ (d + dd - d0) * (d + dd - d0) }, 1.0, 100.0);
                const double e2 = fene(FeneDistSq{ (d - dd - d0) * (d - dd - d0) }, 1.0, 
                100.0);
                const double de = dFeneCoeff(FeneDistSq{ (d - d0) * (d - d0) }, 1.0, 100.0) * (d - d0);

                REQUIRE(de > 0);
                REQUIRE(isfinite(de));

                CHECK(abs(1 - (e1 - e2) / (2 * de * dd)) < 1e-5);
            }
            {
                INFO("3D case.");
                const Vec3d d { -30.0, 0.0, 30.0 };
                const Vec3d d0 { 40.0, 20.0, 0.0 };
                const Vec3d dd { 1e-7, -1e-6, 5e-7 };

                const double e1 = fene(FeneDistSq{ magnitude2(d + dd - d0) }, 1.0, 300.0);
                const double e2 = fene(FeneDistSq{ magnitude2(d - dd - d0) }, 1.0, 300.0);
                const double deCoeff = dFeneCoeff(FeneDistSq{ magnitude2(d - d0) }, 1.0, 300.0);
                REQUIRE(deCoeff > 0);
                REQUIRE(isfinite(deCoeff));
                const auto de = deCoeff * (d - d0);

                CHECK(abs(1 - (e1 - e2) / (2 * dot(de, dd))) < 1e-5);
            }
        }
    }

    SECTION("Simple auxiliary functions") {
        {
            INFO("Integer ceiling division")
            CHECK(ceildiv(-1, 1) == -1);
            CHECK(ceildiv(0, 1) == 0);
            CHECK(ceildiv(1, 1) == 1);

            CHECK(ceildiv(-3, 3) == -1);
            CHECK(ceildiv(-1, 3) == 0);
            CHECK(ceildiv(0, 3) == 0);
            CHECK(ceildiv(1, 3) == 1);
            CHECK(ceildiv(3, 3) == 1);

            CHECK(ceildiv(-3, -3) == 1);
            CHECK(ceildiv(-1, -3) == 1);
            CHECK(ceildiv(0, -3) == 0);
            CHECK(ceildiv(1, -3) == 0);
            CHECK(ceildiv(3, -3) == -1);
        }
    }
}

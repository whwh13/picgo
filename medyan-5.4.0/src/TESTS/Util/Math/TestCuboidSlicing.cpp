#include "catch2/catch.hpp"

#include "MathFunctions.h"
#include "Util/Math/CuboidSlicing.hpp"

using namespace medyan;
using namespace mathfunc;

namespace {

bool planeCuboidSlicingResultEqual(const PlaneCuboidSlicingResult<double>& r1, const PlaneCuboidSlicingResult<double>& r2, double eps) {
    if(abs(r1.volumeIn - r2.volumeIn) > eps) return false;
    size_t s = r1.areaIn.size();
    for(size_t i = 0; i < r1.areaIn.size(); ++i) {
        if(abs(r1.areaIn[i] - r2.areaIn[i]) > eps) return false;
    }
    return true;
}

PlaneCuboidSlicingResult<double> planeUnitCubeSliceByIntersection(double x, double y, double z) {
    auto normal = normalizedVector(Vec3{1.0/x, 1.0/y, 1.0/z});
    auto point = Vec3{x, 0, 0};
    return planeUnitCubeSlice(point, normal);
}

} // namespace

TEST_CASE("Unit cube slicing transition continuity", "[Cuboid Slicing]") {
    /**************************************************************************
    Test continuity under different transitions under the simplified
    conditions, i.e. normal has all non-neg components (no normal flip is
    required, however, the reverse might secretly use a whole flip).
    
    Cube slicing result is obtained by calculating intersections on the axes,
    and the formula are different depending on the positions of the
    intersections.
    **************************************************************************/

    double inEps = 1e-6;
    double outEps = 1e-5;

    // From 3 points inside, to 2 points inside
    {
        auto r1 = planeUnitCubeSliceByIntersection(1 - inEps, 0.5, 0.5);
        auto r2 = planeUnitCubeSliceByIntersection(1 + inEps, 0.5, 0.5);
        CHECK(planeCuboidSlicingResultEqual(r1, r2, outEps));
    }
    {
        auto r1 = planeUnitCubeSliceByIntersection(0.5, 1 - inEps, 0.5);
        auto r2 = planeUnitCubeSliceByIntersection(0.5, 1 + inEps, 0.5);
        CHECK(planeCuboidSlicingResultEqual(r1, r2, outEps));
    }
    {
        auto r1 = planeUnitCubeSliceByIntersection(0.5, 0.5, 1 - inEps);
        auto r2 = planeUnitCubeSliceByIntersection(0.5, 0.5, 1 + inEps);
        CHECK(planeCuboidSlicingResultEqual(r1, r2, outEps));
    }

    // From 3 points inside, to 1 point inside
    {
        auto r1 = planeUnitCubeSliceByIntersection(0.5, 1 - inEps, 1 - inEps);
        auto r2 = planeUnitCubeSliceByIntersection(0.5, 1 + inEps, 1 + inEps);
        CHECK(planeCuboidSlicingResultEqual(r1, r2, outEps));
    }
    {
        auto r1 = planeUnitCubeSliceByIntersection(1 - inEps, 0.5, 1 - inEps);
        auto r2 = planeUnitCubeSliceByIntersection(1 + inEps, 0.5, 1 + inEps);
        CHECK(planeCuboidSlicingResultEqual(r1, r2, outEps));
    }
    {
        auto r1 = planeUnitCubeSliceByIntersection(1 - inEps, 1 - inEps, 0.5);
        auto r2 = planeUnitCubeSliceByIntersection(1 + inEps, 1 + inEps, 0.5);
        CHECK(planeCuboidSlicingResultEqual(r1, r2, outEps));
    }

    // Transitions between the scenario with 1 point inside
    {
        auto r1 = planeUnitCubeSliceByIntersection(0.5, 2 - inEps, 2 - inEps);
        auto r2 = planeUnitCubeSliceByIntersection(0.5, 2 + inEps, 2 + inEps);
        CHECK(planeCuboidSlicingResultEqual(r1, r2, outEps));
    }
    {
        auto r1 = planeUnitCubeSliceByIntersection(2 - inEps, 0.5, 2 - inEps);
        auto r2 = planeUnitCubeSliceByIntersection(2 + inEps, 0.5, 2 + inEps);
        CHECK(planeCuboidSlicingResultEqual(r1, r2, outEps));
    }
    {
        auto r1 = planeUnitCubeSliceByIntersection(2 - inEps, 2 - inEps, 0.5);
        auto r2 = planeUnitCubeSliceByIntersection(2 + inEps, 2 + inEps, 0.5);
        CHECK(planeCuboidSlicingResultEqual(r1, r2, outEps));
    }

    // From 3 points inside, to no point inside
    {
        auto r1 = planeUnitCubeSliceByIntersection(1 - inEps, 1 - inEps, 1 - inEps);
        auto r2 = planeUnitCubeSliceByIntersection(1 + inEps, 1 + inEps, 1 + inEps);
        CHECK(planeCuboidSlicingResultEqual(r1, r2, outEps));
    }

    // Check reverse condition
    {
        auto r1 = planeUnitCubeSliceByIntersection(1.5 - inEps, 1.5 - inEps, 1.5 - inEps);
        auto r2 = planeUnitCubeSliceByIntersection(1.5 + inEps, 1.5 + inEps, 1.5 + inEps);
        CHECK(planeCuboidSlicingResultEqual(r1, r2, outEps));
        CHECK_FALSE(planeCuboidSlicingResultEqual(r1, r2, inEps * 1e-5)); // Check actually changed
    }
    
}

TEST_CASE("Unit cube flipping", "[Cuboid Slicing]") {
    /**************************************************************************
    Test the flipping of normal vector
    **************************************************************************/

    const double nVal = 1.0 / sqrt(3);
    const double pVal = 1.0 / 6.0;

    auto r1 = planeUnitCubeSlice({pVal, pVal, pVal}, {nVal, nVal, nVal});

    // Flip x
    {
        std::array<size_t, 6> comp {{1, 0, 2, 3, 4, 5}};
        auto r2 = planeUnitCubeSlice({{1-pVal, pVal, pVal}}, {{-nVal, nVal, nVal}});
        CHECK(r1.volumeIn == Approx(r2.volumeIn));
        for(size_t i = 0; i < 6; ++i) {
            CHECK(r1.areaIn[i] == Approx(r2.areaIn[comp[i]]));
        }
    }

    // Flip y
    {
        std::array<size_t, 6> comp {{0, 1, 3, 2, 4, 5}};
        auto r2 = planeUnitCubeSlice({{pVal, 1-pVal, pVal}}, {{nVal, -nVal, nVal}});
        CHECK(r1.volumeIn == Approx(r2.volumeIn));
        for(size_t i = 0; i < 6; ++i) {
            CHECK(r1.areaIn[i] == Approx(r2.areaIn[comp[i]]));
        }
    }

    // Flip z
    {
        std::array<size_t, 6> comp {{0, 1, 2, 3, 5, 4}};
        auto r2 = planeUnitCubeSlice({{pVal, pVal, 1-pVal}}, {{nVal, nVal, -nVal}});
        CHECK(r1.volumeIn == Approx(r2.volumeIn));
        for(size_t i = 0; i < 6; ++i) {
            CHECK(r1.areaIn[i] == Approx(r2.areaIn[comp[i]]));
        }
    }

    // Flip yz
    {
        std::array<size_t, 6> comp {{0, 1, 3, 2, 5, 4}};
        auto r2 = planeUnitCubeSlice({{pVal, 1-pVal, 1-pVal}}, {{nVal, -nVal, -nVal}});
        CHECK(r1.volumeIn == Approx(r2.volumeIn));
        for(size_t i = 0; i < 6; ++i) {
            CHECK(r1.areaIn[i] == Approx(r2.areaIn[comp[i]]));
        }
    }

    // Flip zx
    {
        std::array<size_t, 6> comp {{1, 0, 2, 3, 5, 4}};
        auto r2 = planeUnitCubeSlice({{1-pVal, pVal, 1-pVal}}, {{-nVal, nVal, -nVal}});
        CHECK(r1.volumeIn == Approx(r2.volumeIn));
        for(size_t i = 0; i < 6; ++i) {
            CHECK(r1.areaIn[i] == Approx(r2.areaIn[comp[i]]));
        }
    }

    // Flip xy
    {
        std::array<size_t, 6> comp {{1, 0, 3, 2, 4, 5}};
        auto r2 = planeUnitCubeSlice({{1-pVal, 1-pVal, pVal}}, {{-nVal, -nVal, nVal}});
        CHECK(r1.volumeIn == Approx(r2.volumeIn));
        for(size_t i = 0; i < 6; ++i) {
            CHECK(r1.areaIn[i] == Approx(r2.areaIn[comp[i]]));
        }
    }

    // Flip xyz
    {
        std::array<size_t, 6> comp {{1, 0, 3, 2, 5, 4}};
        auto r2 = planeUnitCubeSlice({{1-pVal, 1-pVal, 1-pVal}}, {{-nVal, -nVal, -nVal}});
        CHECK(r1.volumeIn == Approx(r2.volumeIn));
        for(size_t i = 0; i < 6; ++i) {
            CHECK(r1.areaIn[i] == Approx(r2.areaIn[comp[i]]));
        }
    }

}

TEST_CASE("Cube slicing translation and scaling", "[Cuboid Slicing]") {
    /**************************************************************************
    Test the slicing some cube with certain position and size
    **************************************************************************/

    const double nVal = 1.0 / sqrt(3);
    const double pVal = 1.0 / 6.0;
    const double r0Val = 10.0;
    const double s = 2.0;

    auto r = PlaneCubeSlicer() (
        Vec3{r0Val + s * pVal, r0Val + s * pVal, r0Val + s * pVal},
        Vec3{-nVal, -nVal, -nVal},
        Vec3{r0Val, r0Val, r0Val},
        s
    );

    const double exVolumeIn = 47.0 / 6.0;
    const std::array<double, 6> exAreaIn {{3.5, 4.0, 3.5, 4.0, 3.5, 4.0}};

    CHECK(r.volumeIn == Approx(exVolumeIn).epsilon(1e-5));
    for(size_t i = 0; i < 6; ++i) {
        CHECK(r.areaIn[i] == Approx(exAreaIn[i]).epsilon(1e-5));
    }

}

TEST_CASE("Cuboid slicing transition and scaling", "[Cuboid Slicing]") {
    /**************************************************************************
    Test slicing non-cube cuboid
    **************************************************************************/

    const std::array<double, 3> boxSize { 2.0, 4.0, 1.0 };
    const Vec3 normal { 0.0, sqrt(0.5), sqrt(0.5) };
    const Vec3 r0 { 100.0, 100.0, 100.0 };
    const Vec3 point { 100.0, 102.0, 100.0 };

    auto r = PlaneCuboidSlicer() (point, normal, r0, boxSize);

    const double exVolumeIn = 3.0;
    const std::array<double, 6> exAreaIn {{ 1.5, 1.5, 2.0, 0.0, 4.0, 2.0}};

    CHECK(r.volumeIn == Approx(exVolumeIn).epsilon(1e-5));
    for(size_t i = 0; i < 6; ++i) {
        CHECK(r.areaIn[i] == Approx(exAreaIn[i]).epsilon(1e-5));
    }
}

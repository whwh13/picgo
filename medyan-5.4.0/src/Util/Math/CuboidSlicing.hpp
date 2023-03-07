#ifndef MEDYAN_Util_Math_CuboidSlicing_Hpp
#define MEDYAN_Util_Math_CuboidSlicing_Hpp

#include <array>

#include "MathFunctions.h"

namespace medyan {
// The result produced by plane intersecting a cube
template< typename Float = double >
struct PlaneCuboidSlicingResult {
    Float volumeIn;

    std::array<Float, 6> areaIn;
    // in the order of
    // x = x_min, x = x_max,
    // y = y_min, y = y_max,
    // z = z_min, z = z_max

    PlaneCuboidSlicingResult& operator*=(Float a) {
        // Expand the volume and area
        Float a2 = a * a;
        volumeIn *= a2 * a;
        for(Float& eachAreaIn: areaIn) eachAreaIn *= a2;
        return *this;
    }

    PlaneCuboidSlicingResult& operator*=(const std::array<Float, 3>& a) {
        // Expand the volume and area with different aspect ratio
        const std::array<Float, 3> areaFac {{ a[1] * a[2], a[2] * a[0], a[0] * a[1] }};
        for(size_t idx = 0; idx < 3; ++idx) {
            areaIn[2*idx    ] *= areaFac[idx];
            areaIn[2*idx + 1] *= areaFac[idx];
        }
        volumeIn *= a[0] * a[1] * a[2];
        return *this;
    }

    PlaneCuboidSlicingResult& flip(unsigned short int flippedCoords) { // coords (0-7) in binary is zyx
        for(size_t inspect = 0; inspect < 3; ++inspect) {
            if((flippedCoords >> inspect) & 1) {
                std::swap(areaIn[inspect*2], areaIn[inspect*2+1]);
            }
        }
        return *this;
    }
    PlaneCuboidSlicingResult& reverse(Float a) { // a is cube size
        Float a2 = a * a;
        volumeIn = a2 * a - volumeIn;
        for(Float& eachAreaIn: areaIn) eachAreaIn = a2 - eachAreaIn;
        return *this;
    }

};

template< typename Float = double >
inline auto planeUnitCubeSlice( // Cube [0, 1] x [0, 1] x [0, 1]
    const medyan::Vec< 3, Float >& point,  // A point on the plane
    const medyan::Vec< 3, Float >& normal  // Unit normal of the plane pointing outwards
) {
    PlaneCuboidSlicingResult< Float > res{}; // zero-initiate

    // First consider the case where normal has only non-neg components
    auto flippedPoint = point;
    auto flippedNormal = normal;
    unsigned short int flippedCoords = 0;
    for(size_t idx = 0; idx < 3; ++idx) {
        if(flippedNormal[idx] < 0) {
            flippedCoords += (1 << idx);
            flippedPoint[idx] = 1 - flippedPoint[idx];
            flippedNormal[idx] = -flippedNormal[idx];
        }
    }

    Float signedDist = medyan::dot(flippedPoint, flippedNormal);
    Float signedDistMax = medyan::dot(medyan::Vec< 3, Float >{ 1,1,1 }, flippedNormal);

    bool reverse = false;

    if(signedDist <= 0) {
        // do nothing because everything is zero
    } else if(signedDist >= signedDistMax) {
        res.volumeIn = 1;
        res.areaIn = {{1, 1, 1, 1, 1, 1}};
    } else {
        if(signedDist > 0.5 * signedDistMax) {
            // Doing a full flip, plus swapping inside/outside,
            // which leaves the normal unchanged.
            reverse = true;
            flippedCoords ^= ((1 << 3) - 1); // full flip
            for(size_t idx = 0; idx < 3; ++idx) {
                flippedPoint[idx] = 1 - flippedPoint[idx];
            }

            signedDist = signedDistMax - signedDist;
        }

        // Now the signedDist should be within [0, 0.5 signedDistMax]

        // Find the intersections on the axis. The values should be within [0, +inf]
        Float x = signedDist / flippedNormal[0];
        Float y = signedDist / flippedNormal[1];
        Float z = signedDist / flippedNormal[2];

        size_t numGT1 = 0;
        if(x > 1) ++numGT1;
        if(y > 1) ++numGT1;
        if(z > 1) ++numGT1;

        switch(numGT1) {
        case 0:
            {
                res.volumeIn = x * y * z / 6;
                res.areaIn[0] = y * z / 2;
                res.areaIn[2] = z * x / 2;
                res.areaIn[4] = x * y / 2;
            }
            break;

        case 1:
            if(x > 1) {
                Float r = 1 - 1 / x;
                Float y1 = y * r;
                Float z1 = z * r;
                res.volumeIn = (x == std::numeric_limits<Float>::infinity()?
                    y * z / 2:
                    (x * y * z - (x - 1) * y1 * z1) / 6
                );
                res.areaIn[0] = y * z / 2;  res.areaIn[1] = y1 * z1 / 2;
                res.areaIn[2] = (z + z1) / 2;
                res.areaIn[4] = (y + y1) / 2;
            } else if(y > 1) {
                Float r = 1 - 1 / y;
                Float z1 = z * r;
                Float x1 = x * r;
                res.volumeIn = (y == std::numeric_limits<Float>::infinity()?
                    z * x / 2:
                    (x * y * z - (y - 1) * z1 * x1) / 6
                );
                res.areaIn[2] = z * x / 2;  res.areaIn[3] = z1 * x1 / 2;
                res.areaIn[4] = (x + x1) / 2;
                res.areaIn[0] = (z + z1) / 2;
            } else if(z > 1) {
                Float r = 1 - 1 / z;
                Float x1 = x * r;
                Float y1 = y * r;
                res.volumeIn = (z == std::numeric_limits<Float>::infinity()?
                    x * y / 2:
                    (x * y * z - (z - 1) * x1 * y1) / 6
                );
                res.areaIn[4] = x * y / 2;  res.areaIn[5] = x1 * y1 / 2;
                res.areaIn[0] = (y + y1) / 2;
                res.areaIn[2] = (x + x1) / 2;
            }
            break;

        case 2:
            if(x <= 1) {
                Float r1 = 1 - 1 / y;
                Float r2 = 1 - 1 / z;
                Float r3 = r1 + r2 - 1;

                Float y1 = y * r2;
                Float z1 = z * r1;

                Float x1 = x * r1;
                Float x2 = x * r2;
                if(r3 <= 0) {
                    // Neither y nor z can be inf
                    res.volumeIn = (x * y * z - y1 * x2 * (z - 1) - z1 * x1 * (y - 1)) / 6;
                    res.areaIn[0] = (y * z - y1 * (z - 1) - z1 * (y - 1)) / 2;
                    res.areaIn[2] = (x + x2) / 2;   res.areaIn[3] = z1 * x1 / 2;
                    res.areaIn[4] = (x + x1) / 2;   res.areaIn[5] = y1 * x2 / 2;
                } else {
                    Float x3 = x * r3;
                    res.volumeIn = (x + x3) / 2;
                    res.areaIn[0] = 1;
                    res.areaIn[2] = (x + x2) / 2;   res.areaIn[3] = (x1 + x3) / 2;
                    res.areaIn[4] = (x + x1) / 2;   res.areaIn[5] = (x2 + x3) / 2;
                }
            } else if(y <= 1) {
                Float r1 = 1 - 1 / z;
                Float r2 = 1 - 1 / x;
                Float r3 = r1 + r2 - 1;

                Float z1 = z * r2;
                Float x1 = x * r1;

                Float y1 = y * r1;
                Float y2 = y * r2;
                if(r3 <= 0) {
                    // Neither z nor x can be inf
                    res.volumeIn = (x * y * z - z1 * y2 * (x - 1) - x1 * y1 * (z - 1)) / 6;
                    res.areaIn[2] = (z * x - z1 * (x - 1) - x1 * (z - 1)) / 2;
                    res.areaIn[4] = (y + y2) / 2;   res.areaIn[5] = x1 * y1 / 2;
                    res.areaIn[0] = (y + y1) / 2;   res.areaIn[1] = z1 * y2 / 2;
                } else {
                    Float y3 = y * r3;
                    res.volumeIn = (y + y3) / 2;
                    res.areaIn[2] = 1;
                    res.areaIn[4] = (y + y2) / 2;   res.areaIn[5] = (y1 + y3) / 2;
                    res.areaIn[0] = (y + y1) / 2;   res.areaIn[1] = (y2 + y3) / 2;
                }
            } else if(z <= 1) {
                Float r1 = 1 - 1 / x;
                Float r2 = 1 - 1 / y;
                Float r3 = r1 + r2 - 1;

                Float x1 = x * r2;
                Float y1 = y * r1;

                Float z1 = z * r1;
                Float z2 = z * r2;
                if(r3 <= 0) {
                    // Neither x nor y can be inf
                    res.volumeIn = (x * y * z - x1 * z2 * (y - 1) - y1 * z1 * (x - 1)) / 6;
                    res.areaIn[4] = (x * y - x1 * (y - 1) - y1 * (x - 1)) / 2;
                    res.areaIn[0] = (z + z2) / 2;   res.areaIn[1] = y1 * z1 / 2;
                    res.areaIn[2] = (z + z1) / 2;   res.areaIn[3] = x1 * z2 / 2;
                } else {
                    Float z3 = z * r3;
                    res.volumeIn = (z + z3) / 2;
                    res.areaIn[4] = 1;
                    res.areaIn[0] = (z + z2) / 2;   res.areaIn[1] = (z1 + z3) / 2;
                    res.areaIn[2] = (z + z1) / 2;   res.areaIn[3] = (z2 + z3) / 2;
                }
            }
            break;
        
        case 3:
            {
                // no x, y or z can be inf
                Float r1 = 1 - 1 / x;
                Float r2 = 1 - 1 / y;
                Float r3 = 1 - 1 / z;

                Float x1 = x * r2, x2 = x * r3;
                Float y1 = y * r3, y2 = y * r1;
                Float z1 = z * r1, z2 = z * r2;

                res.volumeIn = (x * y * z - y2 * z1 * (x - 1) - z2 * x1 * (y - 1) - x2 * y1 * (x - 1)) / 6;
                res.areaIn[0] = (y * z - y1 * (z - 1) - z2 * (y - 1)) / 2;  res.areaIn[1] = y2 * z1 / 2;
                res.areaIn[2] = (z * x - z1 * (x - 1) - x2 * (z - 1)) / 2;  res.areaIn[3] = z2 * x1 / 2;
                res.areaIn[4] = (x * y - x1 * (y - 1) - y2 * (x - 1)) / 2;  res.areaIn[5] = x2 * y1 / 2;
            }
            break;
        }
    }

    // Finally restore the flipping
    res.flip(flippedCoords);
    if(reverse) res.reverse(1);

    return res;
}

struct PlaneCubeSlicer {
    template< typename Float = double >
    auto operator() (
        const medyan::Vec< 3, Float >& point,  // A point on the plane
        const medyan::Vec< 3, Float >& normal, // Unit normal of the plane pointing outwards
        const medyan::Vec< 3, Float >& r0,     // (x_min, y_min, z_min) of the cube
        Float                            a       // Side length of the cube
    ) {
        PlaneCuboidSlicingResult< Float > res;
        res = planeUnitCubeSlice(
            (point - r0) * (1.0 / a),
            normal
        );
        res *= a;
        return res;
    }
};

struct PlaneCuboidSlicer {
    template< typename Float = double >
    auto operator() (
        const medyan::Vec< 3, Float >& point,   // A point on the plane
        const medyan::Vec< 3, Float >& normal,  // Unit normal of the plane pointing outwards
        const medyan::Vec< 3, Float >& r0,      // (x_min, y_min, z_min) of the cuboid
        const std::array<Float, 3>&      a        // Edge length of the cuboid
    ) {
        PlaneCuboidSlicingResult< Float > res;
        const medyan::Vec< 3, Float > pointInUnitCube {
            (point[0] - r0[0]) / a[0],
            (point[1] - r0[1]) / a[1],
            (point[2] - r0[2]) / a[2]
        };
        const auto normalInUnitCube = medyan::normalizedVector(medyan::Vec< 3, Float > {
            normal[0] * a[0],
            normal[1] * a[1],
            normal[2] * a[2]
        });
        res = planeUnitCubeSlice(pointInUnitCube, normalInUnitCube);
        res *= a;
        return res;
    }
};

} // namespace medyan

#endif

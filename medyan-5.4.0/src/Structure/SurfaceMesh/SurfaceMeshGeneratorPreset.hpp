#ifndef MEDYAN_Structure_SurfaceMesh_SurfaceMeshGeneratorPreset_hpp
#define MEDYAN_Structure_SurfaceMesh_SurfaceMeshGeneratorPreset_hpp

#include <algorithm> // min, min_element
#include <sstream>
#include <stdexcept> // runtime_error
#include <string>
#include <vector>

#include "Structure/SurfaceMesh/SurfaceMeshGenerator.hpp"
#include "Util/Io/Log.hpp"

namespace medyan::mesh_gen {

template< typename Float = double >
inline typename MarchingTetrahedraGenerator<Float>::Result
generateEllipsoid(const Vec< 3, Float >& center, const Vec< 3, Float >& semiAxes) {
    using namespace mathfunc;

    constexpr Float  maxCubeSize               = 60.0;
    constexpr size_t minCubesEachSemiDirection = 20;

    const Float cubeSize = std::min(
        *std::min_element(semiAxes.begin(), semiAxes.end()) / minCubesEachSemiDirection,
        maxCubeSize
    );
    const Vec< 3, Float > boxOrigin {
        center[0] - semiAxes[0] - cubeSize,
        center[1] - semiAxes[1] - cubeSize,
        center[2] - semiAxes[2] - cubeSize
    };
    const std::array< std::size_t, 3 > numCubes {
        size_t(2 * semiAxes[0] / cubeSize) + 2,
        size_t(2 * semiAxes[1] / cubeSize) + 2,
        size_t(2 * semiAxes[2] / cubeSize) + 2
    };

    return MarchingTetrahedraGenerator< Float >(cubeSize, boxOrigin, numCubes)(
        [&](const Vec< 3, Float >& v) {
            const auto x = v - center;
            return x[0] * x[0] / (semiAxes[0] * semiAxes[0]) +
                   x[1] * x[1] / (semiAxes[1] * semiAxes[1]) +
                   x[2] * x[2] / (semiAxes[2] * semiAxes[2]) - 1;
        }
    );
}

template< typename Float = double >
inline auto generatePlane(
    const Vec< 3, Float >& pointOnPlane, const Vec< 3, Float >& normal,
    const Vec< 3, Float >& boundingBoxOrigin, const Vec< 3, Float >& boundingBoxSize
) {
    using namespace mathfunc;

    constexpr Float cubeSize = 40.0;

    const std::array< std::size_t, 3 > numCubes {
        size_t(boundingBoxSize[0] / cubeSize),
        size_t(boundingBoxSize[1] / cubeSize),
        size_t(boundingBoxSize[2] / cubeSize)
    };
    // Overall size might be smaller than specified bounding box size

    const auto un = normalizedVector(normal);

    return MarchingTetrahedraGenerator< Float >(cubeSize, boundingBoxOrigin, numCubes)(
        [&](const Vec< 3, Float >& v) {
            return dot(v - pointOnPlane, un);
        }
    );
}

template< typename Float = double >
inline typename MarchingTetrahedraGenerator<Float>::Result
generateMeshViaParams(const std::vector< std::string >& param) {
    if(param.empty()) {
        LOG(ERROR) << "Mesh initialization parameters should not be empty";
        throw std::runtime_error("Empty mesh init param");
    }

    if(param[0] == "ELLIPSOID") {
        if(param.size() < 7) {
            log::error("Mesh ellipsoid initializer has 6 parameters");
            log::info("Mesh ellipsoid initializer: ELLIPSOID center_x center_y center_z a b c");
            throw std::runtime_error("Insufficient parameters in mesh ellipsoid initializer");
        }

        Float value[6];
        for(size_t i = 0; i < 6; ++i) {
            std::istringstream iss(param[i+1]);
            iss >> value[i];
        }
        return generateEllipsoid< Float >({value[0], value[1], value[2]}, {value[3], value[4], value[5]});
    }

    if(param[0] == "PLANE") {
        if(param.size() < 13) {
            log::error("Mesh plane initializer has 12 parameters");
            log::info("Mesh plane initializer: PLANE center_xyz... normal_xyz... box_origin_xyz... box_size_abc...");
            throw std::runtime_error("Insufficient parameters in mesh plane initializer");
        }

        Float value[12];
        for(size_t i = 0; i < 12; ++i) {
            std::istringstream iss(param[i+1]);
            iss >> value[i];
        }
        return generatePlane< Float >(
            {value[0], value[1], value[2]},
            {value[3], value[4], value[5]},
            {value[6], value[7], value[8]},
            {value[9], value[10], value[11]}
        );
    }

    log::error("Unrecognized mesh initializer");
    throw std::runtime_error("Unrecognized mesh initializer");
}

} // namespace medyan::mesh_gen

#endif

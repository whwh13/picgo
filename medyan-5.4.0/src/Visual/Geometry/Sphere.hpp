#ifndef MEDYAN_Visual_Geometry_Sphere_Hpp
#define MEDYAN_Visual_Geometry_Sphere_Hpp

#include <array>
#include <cstdint> // uint_fast8_t
#include <tuple>
#include <vector>

#include "Util/Math/Vec.hpp"

namespace medyan::visual {

// This function transforms a bead to GL_TRIANGLES compatible vertex coord list
// and index list.
template< typename Float = double >
struct SphereUv {
    using CoordType = Vec< 3, Float >;

    // ShapeCache stores the coordinates on a unit sphere, as well as the
    // triangle indices
    struct ShapeCache {
        int longitudeSegs;
        int latitudeSegs;
        std::vector< CoordType > vertices; // poles at index size() - 2 and size() - 1
        std::vector< std::array< int, 3 > > triInd;
    };

    Float radius; // of circumscribing sphere
    int   longitudeSegs; // at least 3, about twice of latitudeSegs
    int   latitudeSegs; // at least 2

    ShapeCache makeCache() const {
        ShapeCache res { longitudeSegs, latitudeSegs };

        const int numVertices = longitudeSegs * (latitudeSegs - 1) + 2; // 2 means polar vertices
        const int numTriangles = longitudeSegs * (2 * latitudeSegs - 2);

        res.vertices.resize(numVertices);
        res.triInd.resize(numTriangles);

        const auto dPhi = 2 * M_PI / longitudeSegs;
        const auto dTheta = M_PI / latitudeSegs;

        // Calculate vertices
        for(std::uint_fast8_t nTheta = 1; nTheta < latitudeSegs; ++nTheta) {
            for(std::uint_fast8_t nPhi = 0; nPhi < longitudeSegs; ++nPhi) {
                res.vertices[(nTheta - 1) * longitudeSegs + nPhi] = CoordType {
                    static_cast<Float>(std::sin(nTheta * dTheta) * std::cos(nPhi * dPhi)),
                    static_cast<Float>(std::sin(nTheta * dTheta) * std::sin(nPhi * dPhi)),
                    static_cast<Float>(std::cos(nTheta * dTheta))
                };
            }
        }
        res.vertices[numVertices - 2] = CoordType { (Float)0.0, (Float)0.0, (Float)1.0    }; // Nouth pole
        res.vertices[numVertices - 1] = CoordType { (Float)0.0, (Float)0.0, (Float)(-1.0) }; // Sorth pole

        // Register triangles (GL_TRIANGLES)
        // from North pole to South pole
        for(std::uint_fast8_t nTheta = 0; nTheta < latitudeSegs; ++nTheta) {
            const int numPrev = (nTheta == 0 ? 0 : (2 * nTheta - 1) * longitudeSegs);
            for(std::uint_fast8_t nPhi = 0; nPhi < longitudeSegs; ++nPhi) {
                if(nTheta == 0) {
                    res.triInd[numPrev + nPhi][0] = numVertices - 2; // North pole
                    res.triInd[numPrev + nPhi][1] = nTheta * longitudeSegs + nPhi;
                    res.triInd[numPrev + nPhi][2] = nTheta * longitudeSegs + (nPhi + 1) % longitudeSegs;
                }
                else if(nTheta == latitudeSegs - 1) {
                    res.triInd[numPrev + nPhi][0] = (nTheta - 1) * longitudeSegs + nPhi;
                    res.triInd[numPrev + nPhi][1] = numVertices - 1; // South pole
                    res.triInd[numPrev + nPhi][2] = (nTheta - 1) * longitudeSegs + (nPhi + 1) % longitudeSegs;
                }
                else {
                    res.triInd[numPrev + 2 * nPhi][0]     = (nTheta - 1) * longitudeSegs + nPhi;
                    res.triInd[numPrev + 2 * nPhi][1]     =  nTheta      * longitudeSegs + nPhi;
                    res.triInd[numPrev + 2 * nPhi][2]     = (nTheta - 1) * longitudeSegs + (nPhi + 1) % longitudeSegs;
                    res.triInd[numPrev + 2 * nPhi + 1][0] = (nTheta - 1) * longitudeSegs + (nPhi + 1) % longitudeSegs;
                    res.triInd[numPrev + 2 * nPhi + 1][1] =  nTheta      * longitudeSegs + nPhi;
                    res.triInd[numPrev + 2 * nPhi + 1][2] =  nTheta      * longitudeSegs + (nPhi + 1) % longitudeSegs;
                }
            }
        }

        return res;
    } // End makeCache()

    // This function transforms a point in space to the mesh of a ball
    // Parameters
    //   - coords:  the container where the coordinates of beads can be found
    //   - cache:   the shape cache
    auto generate(const CoordType& coord, const ShapeCache& cache) const {
        using namespace mathfunc;

        constexpr CoordType e0 = { 1.0, 0.0, 0.0 };
        constexpr CoordType e1 = { 0.0, 1.0, 0.0 };
        const CoordType e2 = cross(e0, e1);

        const int numVertices = cache.vertices.size();

        // Results
        std::vector< CoordType > vertices(numVertices); // poles at index nv - 2 and nv - 1

        // Calculate vertices
        for(int i = 0; i < numVertices; ++i) {
            vertices[i] = coord + radius * (
                e0 * cache.vertices[i][0] +
                e1 * cache.vertices[i][1] +
                e2 * cache.vertices[i][2]
            );
        }

        return std::make_tuple(vertices, cache.triInd);
    }

    // Finds an estimate of number of triangles returned by "generate".
    int estimateNumTriangles() const {
        return longitudeSegs * (2 * latitudeSegs - 2);
    }

};

} // namespace medyan::visual

#endif

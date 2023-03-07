#ifndef MEDYAN_Structure_SurfaceMesh_SurfaceMeshGenerator_hpp
#define MEDYAN_Structure_SurfaceMesh_SurfaceMeshGenerator_hpp

#include <algorithm> // max, min
#include <array>
#include <cstdint> // uint_fast8_t
#include <type_traits> // integral_constant, enable_if
#include <vector>

#include "Util/Environment.hpp"
#include "Util/Math/Vec.hpp"

namespace medyan::mesh_gen {

namespace indexer {

using uif8 = std::uint_fast8_t;

//-------------------------------------------------------------------------
// Indexing in a tetrahedron:
//
// The vertices of the tetrahedra are v0, v1, v2, v3, which must satisfy
//     r01 x r12 . r23 > 0 (specifically, the cuboid volume)
//
// The edges in a tetrahedra are in the following order:
//     01, 02, 03, 12, 13, 23
//-------------------------------------------------------------------------

namespace internal {

    template< uif8 v0, uif8 v1 > struct TetraEdgeIndexFromVertexIndex;
    template<> struct TetraEdgeIndexFromVertexIndex< 0, 1 > : std::integral_constant< uif8, 0 > {};
    template<> struct TetraEdgeIndexFromVertexIndex< 0, 2 > : std::integral_constant< uif8, 1 > {};
    template<> struct TetraEdgeIndexFromVertexIndex< 0, 3 > : std::integral_constant< uif8, 2 > {};
    template<> struct TetraEdgeIndexFromVertexIndex< 1, 2 > : std::integral_constant< uif8, 3 > {};
    template<> struct TetraEdgeIndexFromVertexIndex< 1, 3 > : std::integral_constant< uif8, 4 > {};
    template<> struct TetraEdgeIndexFromVertexIndex< 2, 3 > : std::integral_constant< uif8, 5 > {};
    template< uif8 v0, uif8 v1 >
    constexpr uif8 tetraEdgeIndexFromVertexIndex = TetraEdgeIndexFromVertexIndex< v0, v1 >::value;

} // namespace internal

struct TetraTriangleIntersectionIndex {
    uif8 vertexIndex[3][2];
    uif8 edgeIndex[3];
};
template< uif8 v00, uif8 v01, uif8 v10, uif8 v11, uif8 v20, uif8 v21 >
constexpr auto getTetraTriangleIntersectionIndex() {
    constexpr uif8 e0 = internal::tetraEdgeIndexFromVertexIndex< v00, v01 >;
    constexpr uif8 e1 = internal::tetraEdgeIndexFromVertexIndex< v10, v11 >;
    constexpr uif8 e2 = internal::tetraEdgeIndexFromVertexIndex< v20, v21 >;
    return TetraTriangleIntersectionIndex {
        { {v00, v01}, {v10, v11}, {v20, v21} },
        { e0, e1, e2 }
    };
}

struct TetraQuadIntersectionIndex {
    uif8 vertexIndex[4][2];
    uif8 edgeIndex[4];
};
template< uif8 v00, uif8 v01, uif8 v10, uif8 v11, uif8 v20, uif8 v21, uif8 v30, uif8 v31 >
constexpr auto getTetraQuadIntersectionIndex() {
    constexpr uif8 e0 = internal::tetraEdgeIndexFromVertexIndex< v00, v01 >;
    constexpr uif8 e1 = internal::tetraEdgeIndexFromVertexIndex< v10, v11 >;
    constexpr uif8 e2 = internal::tetraEdgeIndexFromVertexIndex< v20, v21 >;
    constexpr uif8 e3 = internal::tetraEdgeIndexFromVertexIndex< v30, v31 >;
    return TetraQuadIntersectionIndex {
        { {v00, v01}, {v10, v11}, {v20, v21}, {v30, v31} },
        { e0, e1, e2, e3 }
    };
}

} // namespace indexer

template< typename Float = double >
class MarchingTetrahedraGenerator {
public:
    using coordinate_type = medyan::Vec< 3, Float >;
    using small_size_t = std::uint_fast8_t;

    // This struct can be replaced by std::optional since C++17
    struct TetraEdgeData {
        int  vertexIdxInMesh;
        bool hasIntersection = false;
    };

    // The result type
    // Note:
    //   - The terms "vertex" and "triangle" refer to the elements in the
    //     resulting surface mesh, not the tetrahedra system.
    struct Result {
        std::vector< coordinate_type > vertexCoordinateList;
        std::vector< std::array< int, 3 > > triangleList; // Each triangle consists of vertex indices
    };

    MarchingTetrahedraGenerator(
        Float cubeSize,
        const coordinate_type& boundingBoxOrigin,
        const std::array< std::size_t, 3 >& numCubes
    ) : _cuboidSize{ cubeSize, cubeSize, cubeSize },
        _numCuboids(numCubes),
        _boundingBoxOrigin(boundingBoxOrigin)
    { }

    // The function that generates the meshwork
    // Template parameters:
    //   - FieldFunc: a functor that takes a point and returns given how far it is outside the surface
    template< typename FieldFunc >
    Result operator()(FieldFunc&& func) const {
        using std::size_t;

        // Temporary data
        std::vector< Float > distValue(_getVertexListSize()); // using vertex indexing
        std::vector< TetraEdgeData > edgeData(_getEdgeListSize()); // using edge indexing

        // Result
        Result res;

        // Helper function: computing intersections for a tetrahedron
        const auto calcTetra = [&](size_t nx, size_t ny, size_t nz, small_size_t tetIdx) {
            using namespace std;

            // switch sign of 4 switches and decide which procedure to follow (0.0 count as positive)
            // all pos / neg: nothing
            // one pos: gen triangle pointing pos
            // one neg: gen triangle pointing neg
            // 2 pos 2 neg: gen 2 triangles

            const array< size_t, 4 > vertexIndices {
                _getVertexIdxInTetra(nx, ny, nz, tetIdx, 0),
                _getVertexIdxInTetra(nx, ny, nz, tetIdx, 1),
                _getVertexIdxInTetra(nx, ny, nz, tetIdx, 2),
                _getVertexIdxInTetra(nx, ny, nz, tetIdx, 3)
            };
            const array< Float, 4 > vertexValues {
                distValue[vertexIndices[0]],
                distValue[vertexIndices[1]],
                distValue[vertexIndices[2]],
                distValue[vertexIndices[3]]
            };
            const array< Vec< 3, Float >, 4 > vertexCoords {
                _vertexCoordinateInTetra(nx, ny, nz, tetIdx, 0),
                _vertexCoordinateInTetra(nx, ny, nz, tetIdx, 1),
                _vertexCoordinateInTetra(nx, ny, nz, tetIdx, 2),
                _vertexCoordinateInTetra(nx, ny, nz, tetIdx, 3)
            };
            const array< size_t, 6 > edgeIndices {
                _getEdgeIdxInTetra(nx, ny, nz, tetIdx, 0),
                _getEdgeIdxInTetra(nx, ny, nz, tetIdx, 1),
                _getEdgeIdxInTetra(nx, ny, nz, tetIdx, 2),
                _getEdgeIdxInTetra(nx, ny, nz, tetIdx, 3),
                _getEdgeIdxInTetra(nx, ny, nz, tetIdx, 4),
                _getEdgeIdxInTetra(nx, ny, nz, tetIdx, 5)
            };

            const auto newIntersect = [&](small_size_t eIdx, small_size_t v0Idx, small_size_t v1Idx) {
                if(!edgeData[edgeIndices[eIdx]].hasIntersection) {
                    // Compute position
                    Float pos = vertexValues[v0Idx] / (vertexValues[v0Idx] - vertexValues[v1Idx]);
                    pos = std::min(1 - _minPositionShift, std::max(_minPositionShift, pos));
                    const auto newCoord = vertexCoords[v0Idx] * (1 - pos) + vertexCoords[v1Idx] * pos;

                    // Add new vertex
                    res.vertexCoordinateList.push_back(newCoord);

                    // Update edge data
                    edgeData[edgeIndices[eIdx]].hasIntersection = true;
                    edgeData[edgeIndices[eIdx]].vertexIdxInMesh = res.vertexCoordinateList.size() - 1;
                } // else do nothing
            };
            const auto newTriangle = [&](indexer::TetraTriangleIntersectionIndex idx) {
                // Create intersections
                newIntersect(idx.edgeIndex[0], idx.vertexIndex[0][0], idx.vertexIndex[0][1]);
                newIntersect(idx.edgeIndex[1], idx.vertexIndex[1][0], idx.vertexIndex[1][1]);
                newIntersect(idx.edgeIndex[2], idx.vertexIndex[2][0], idx.vertexIndex[2][1]);

                // Create a new triangle
                res.triangleList.push_back({
                    edgeData[edgeIndices[idx.edgeIndex[0]]].vertexIdxInMesh,
                    edgeData[edgeIndices[idx.edgeIndex[1]]].vertexIdxInMesh,
                    edgeData[edgeIndices[idx.edgeIndex[2]]].vertexIdxInMesh
                });
            };
            const auto newQuad = [&](indexer::TetraQuadIntersectionIndex idx) {
                // Create intersections
                newIntersect(idx.edgeIndex[0], idx.vertexIndex[0][0], idx.vertexIndex[0][1]);
                newIntersect(idx.edgeIndex[1], idx.vertexIndex[1][0], idx.vertexIndex[1][1]);
                newIntersect(idx.edgeIndex[2], idx.vertexIndex[2][0], idx.vertexIndex[2][1]);
                newIntersect(idx.edgeIndex[3], idx.vertexIndex[3][0], idx.vertexIndex[3][1]);

                // Create a new quad by adding 2 triangles with vertices on edges 012 and edges 230
                //
                //     e3 ------- e0
                //     |        / |
                //     |     /    |
                //     |  /       |
                //     e2 ------- e1
                res.triangleList.push_back({
                    edgeData[edgeIndices[idx.edgeIndex[0]]].vertexIdxInMesh,
                    edgeData[edgeIndices[idx.edgeIndex[1]]].vertexIdxInMesh,
                    edgeData[edgeIndices[idx.edgeIndex[2]]].vertexIdxInMesh
                });
                res.triangleList.push_back({
                    edgeData[edgeIndices[idx.edgeIndex[2]]].vertexIdxInMesh,
                    edgeData[edgeIndices[idx.edgeIndex[3]]].vertexIdxInMesh,
                    edgeData[edgeIndices[idx.edgeIndex[0]]].vertexIdxInMesh
                });
            };

            small_size_t cond = 0; // 0bxxxx <-- each bit is 1 if sign is pos, 0 if sign is neg
            for(small_size_t i = 0; i < 4; ++i) {
                cond <<= 1;
                cond |= (vertexValues[i] >= 0.0 ? 1 : 0);
            }

            switch(cond) {
            case 0b0000:
            case 0b1111:
                break;

            case 0b0001:
                // intersects 03, 13, 23, triangle (03, 13, 23)
                newTriangle(indexer::getTetraTriangleIntersectionIndex<
                    0, 3, 1, 3, 2, 3
                >());
                break;
            case 0b1110:
                // intersects 23, 13, 03, triangle (23, 13, 03)
                newTriangle(indexer::getTetraTriangleIntersectionIndex<
                    2, 3, 1, 3, 0, 3
                >());
                break;

            case 0b0010:
                // intersects 12, 02, 23, triangle (12, 02, 23)
                newTriangle(indexer::getTetraTriangleIntersectionIndex<
                    1, 2, 0, 2, 2, 3
                >());
                break;
            case 0b1101:
                // intersects 23, 02, 12, triangle (23, 02, 12)
                newTriangle(indexer::getTetraTriangleIntersectionIndex<
                    2, 3, 0, 2, 1, 2
                >());
                break;

            case 0b0100:
                // intersects 12, 13, 01, triangle (12, 13, 01)
                newTriangle(indexer::getTetraTriangleIntersectionIndex<
                    1, 2, 1, 3, 0, 1
                >());
                break;
            case 0b1011:
                // intersects 01, 13, 12, triangle (01, 13, 12)
                newTriangle(indexer::getTetraTriangleIntersectionIndex<
                    0, 1, 1, 3, 1, 2
                >());
                break;

            case 0b1000:
                // intersects 03, 02, 01, triangle (03, 02, 01)
                newTriangle(indexer::getTetraTriangleIntersectionIndex<
                    0, 3, 0, 2, 0, 1
                >());
                break;
            case 0b0111:
                // intersects 01, 02, 03, triangle (01, 02, 03)
                newTriangle(indexer::getTetraTriangleIntersectionIndex<
                    0, 1, 0, 2, 0, 3
                >());
                break;

            case 0b0011:
                // intersects 02, 03, 13, 12, 2 triangles
                newQuad(indexer::getTetraQuadIntersectionIndex<
                    0, 2, 0, 3, 1, 3, 1, 2
                >());
                break;
            case 0b1100:
                // intersects 12, 13, 03, 02, 2 triangles
                newQuad(indexer::getTetraQuadIntersectionIndex<
                    1, 2, 1, 3, 0, 3, 0, 2
                >());
                break;

            case 0b0101:
                // intersects 01, 12, 23, 03, 2 triangles
                newQuad(indexer::getTetraQuadIntersectionIndex<
                    0, 1, 1, 2, 2, 3, 0, 3
                >());
                break;
            case 0b1010:
                // intersects 03, 23, 12, 01, 2 triangles
                newQuad(indexer::getTetraQuadIntersectionIndex<
                    0, 3, 2, 3, 1, 2, 0, 1
                >());
                break;

            case 0b0110:
                // intersects 01, 02, 23, 13, 2 triangles
                newQuad(indexer::getTetraQuadIntersectionIndex<
                    0, 1, 0, 2, 2, 3, 1, 3
                >());
                break;
            case 0b1001:
                // intersects 13, 23, 02, 01, 2 triangles
                newQuad(indexer::getTetraQuadIntersectionIndex<
                    1, 3, 2, 3, 0, 2, 0, 1
                >());
                break;

            } // switch(cond)
        }; // End of function void calcTetra(...)

        // Precompute phase
        for(size_t nx = 0; nx <= _numCuboids[0]; ++nx)
            for(size_t ny = 0; ny <= _numCuboids[1]; ++ny)
                for(size_t nz = 0; nz <= _numCuboids[2]; ++nz)
                    distValue[_getVertexIdx(nx, ny, nz)] = func(_vertexCoordinate(nx, ny, nz));

        // Generating triangles
        for(size_t nx = 0; nx < _numCuboids[0]; ++nx)
            for(size_t ny = 0; ny < _numCuboids[1]; ++ny)
                for(size_t nz = 0; nz < _numCuboids[2]; ++nz)
                    for(small_size_t tetIdx = 0; tetIdx < _numTetrahedraPerCuboid; ++tetIdx)
                        calcTetra(nx, ny, nz, tetIdx);

        return res;

    } // Result operator()(...)

private:
    //-------------------------------------------------------------------------
    // Indexing in cuboid:
    //
    // Number of vertices: (nx + 1)(ny + 1)(nz + 1)
    //
    // Tetrahedra in a cuboid is ordered in the following order (labeled using r's):
    //     ijk, (i+j)(-i)(k+i), jki, (j+k)(-j)(i+j), kij, (k+i)(-k)(j+k)
    //
    // Edges in a cuboid is simply indexed by 0bxyz(edge direction) - 1
    //-------------------------------------------------------------------------
    static constexpr small_size_t _numTetrahedraPerCuboid = 6;
    static constexpr small_size_t _numTetraEdgesPerCuboid = 7; // index = 0bxyz - 1

    // Local vertex index [local tetra idx (6)][vertex idx (4)]
#ifdef COMPILER_GCC
    // Compiler bug: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=92625
    static constexpr std::uint_fast16_t _tetraVertexLocalIndex[6][4] {
#else
    static constexpr small_size_t _tetraVertexLocalIndex[6][4] {
#endif
        {0b000, 0b100, 0b110, 0b111},
        {0b000, 0b110, 0b010, 0b111},
        {0b000, 0b010, 0b011, 0b111},
        {0b000, 0b011, 0b001, 0b111},
        {0b000, 0b001, 0b101, 0b111},
        {0b000, 0b101, 0b100, 0b111}
    };
    // Local edge index [local tetra idx (6)][edge idx (6)]
    // Value:  0b  xyz                  xyz
    //             ^ starting position  ^ edge direction (positive)
    // Real index in cuboid = 0bxyz(direction) - 1
#ifdef COMPILER_GCC
    // Compiler bug: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=92625
    static constexpr std::uint_fast16_t _tetraEdgeLocalIndex[6][6] {
#else
    static constexpr small_size_t _tetraEdgeLocalIndex[6][6] {
#endif
        {0b000'100, 0b000'110, 0b000'111, 0b100'010, 0b100'011, 0b110'001},
        {0b000'110, 0b000'010, 0b000'111, 0b010'100, 0b110'001, 0b010'101},
        {0b000'010, 0b000'011, 0b000'111, 0b010'001, 0b010'101, 0b011'100},
        {0b000'011, 0b000'001, 0b000'111, 0b001'010, 0b011'100, 0b001'110},
        {0b000'001, 0b000'101, 0b000'111, 0b001'100, 0b001'110, 0b101'010},
        {0b000'101, 0b000'100, 0b000'111, 0b100'001, 0b101'010, 0b100'011}
    };

    // Parameters
    Float                           _minPositionShift = 1e-5; // The position will have a minimal shift from either end
    coordinate_type                 _cuboidSize;
    std::array< std::size_t, 3 >    _numCuboids;
    coordinate_type                 _boundingBoxOrigin;

    // Indexers
    auto _getVertexListSize() const {
        return (_numCuboids[0] + 1) * (_numCuboids[1] + 1) * (_numCuboids[2] + 1);
    }
    auto _getVertexIdx(std::size_t nx, std::size_t ny, std::size_t nz) const {
        std::size_t res = nx;
        res *= _numCuboids[1] + 1; res += ny;
        res *= _numCuboids[2] + 1; res += nz;
        return res;
    }
    auto _getVertexIdxInTetra(std::size_t nx, std::size_t ny, std::size_t nz, small_size_t tetIdx, small_size_t vtxIdx) const {
        const small_size_t i = _tetraVertexLocalIndex[tetIdx][vtxIdx];
        return _getVertexIdx(
            nx + ((i >> 2) & 1),
            ny + ((i >> 1) & 1),
            nz + ( i       & 1)
        );
    }
    auto _getCuboidListSize() const {
        return (_numCuboids[0] + 1) * (_numCuboids[1] + 1) * (_numCuboids[2] + 1);
    }
    auto _getCuboidIdx(std::size_t nx, std::size_t ny, std::size_t nz) const {
        std::size_t res = nx;
        res *= _numCuboids[1] + 1; res += ny;
        res *= _numCuboids[2] + 1; res += nz;
        return res;
    }
    auto _getEdgeListSize() const {
        return _numTetraEdgesPerCuboid * _getCuboidListSize();
    }
    auto _getEdgeIdx(std::size_t nx, std::size_t ny, std::size_t nz, small_size_t edgeIdx) const {
        return _numTetraEdgesPerCuboid * _getCuboidIdx(nx, ny, nz) + edgeIdx;
    }
    auto _getEdgeIdxInTetra(std::size_t nx, std::size_t ny, std::size_t nz, small_size_t tetIdx, small_size_t edgeIdx) const {
        const small_size_t i = _tetraEdgeLocalIndex[tetIdx][edgeIdx];
        const small_size_t idxInCuboid = (i & 0b111) - 1;
        return _getEdgeIdx(
            nx + ((i >> 5) & 1),
            ny + ((i >> 4) & 1),
            nz + ((i >> 3) & 1),
            idxInCuboid
        );
    }

    // Coordinates
    auto _vertexCoordinate(std::size_t nx, std::size_t ny, std::size_t nz) const {
        return medyan::Vec< 3, Float > {
            _boundingBoxOrigin[0] + _cuboidSize[0] * nx,
            _boundingBoxOrigin[1] + _cuboidSize[1] * ny,
            _boundingBoxOrigin[2] + _cuboidSize[2] * nz
        };
    }
    auto _vertexCoordinateInTetra(std::size_t nx, std::size_t ny, std::size_t nz, small_size_t tetIdx, small_size_t vtxIdx) const {
        const small_size_t i = _tetraVertexLocalIndex[tetIdx][vtxIdx];
        return _vertexCoordinate(
            nx + ((i >> 2) & 1),
            ny + ((i >> 1) & 1),
            nz + ( i       & 1)
        );
    }
};

} // namespace medyan::mesh_gen

#endif

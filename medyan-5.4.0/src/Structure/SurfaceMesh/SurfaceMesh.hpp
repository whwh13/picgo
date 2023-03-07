#ifndef MEDYAN_SurfaceMesh_hpp
#define MEDYAN_SurfaceMesh_hpp

#include <algorithm>
#include <array>
#include <cstddef> // ptrdiff_t
#include <cstdint> // uint_fast8_t
#include <iterator>
#include <optional>
#include <sstream>
#include <stdexcept> // runtime_error
#include <string>
#include <type_traits>
#include <utility> // move, pair
#include <vector>

#include "common.h" // Index, Size
#include "Util/Io/Log.hpp"

namespace medyan {

/******************************************************************************
The data structure for an orientable, manifold 2d triangular meshwork in 3d
space. The data structure only provides topological relationship, and is
completed by the Attribute class.

The connectivity is similar to CGAL library, which is halfedge based.

target:     The vertex that the halfedge points to.
opposite:   The halfedge on the same edge but in the opposite direction.
face:       The triangle that the halfedge is circling ccw.
next:       The next halfedge in the current face circling ccw.
prev:       The previous halfedge in the current face circling ccw.
edge:       The undirected edge of this halfedge.

All the other elements must have at least one index pointing to a halfedge.
******************************************************************************/

// The DeletableVector is a thin wrapper for std::vector.
// This container is designed for types with assignment operators.
// When an element is removed, instead of doing vector::erase,
// it essentially swaps the element with the last one, and pops the vector.
template< typename T > class DeletableVector {

private:
    std::vector< T > value_;

public:

    using iterator       = typename std::vector< T >::iterator;
    using const_iterator = typename std::vector< T >::const_iterator;
    iterator       begin() noexcept       { return value_.begin(); }
    const_iterator begin() const noexcept { return value_.begin(); }
    iterator       end() noexcept       { return value_.end(); }
    const_iterator end() const noexcept { return value_.end(); }

    // Insert a new element. Returns the new index.
    Index insert() {
        value_.emplace_back();
        return value_.size() - 1;
    }

    //-------------------------------------------------------------------------
    // Remove an element from the container.
    // Returns the new index of the previous last element.
    // Might change the position of certain elements.
    // If index for deletion is out of range, the behavior is undefined.
    //
    // Uses swap-delete algorithm:
    //   - If the item is the last in the vector, pop it;
    //   - Otherwise, swap the item with the last one and pop back. (As a result,
    //     the indices pointing to the last element will be INVALIDATED, so the
    //     caller must manually retarget all indices to the last element.)
    struct EraseResult {
        Size newSize;
        std::optional< Index > newIndex;

        template< typename... IndexT >
        auto track(IndexT... oldIndices) const {
            using IT = std::array< Index, sizeof...(IndexT) >;
            IT oi = {{ oldIndices ... }};
            IT res;
            for(Index i = 0; i < oi.size(); ++i) {
                if(oi[i] < newSize) res[i] = oi[i];
                else res[i] = *newIndex; // In this case newIndex must have value
            }
            return res;
        }
    };
    template< typename Retargeter > // The retargeter must implement operation()(from, to)
    auto erase(Index index, Retargeter&& r) {
        EraseResult res {};

        const Index lastIndex = value_.size() - 1;
        if(index == lastIndex) {
            value_.pop_back();
        } else {
            // Move value from lastIndex to index
            value_[index] = std::move(value_[lastIndex]);
            value_.pop_back();
            r(lastIndex, index);
            res.newIndex = index;
        }

        res.newSize = value_.size();
        return res;
    }

    //-------------------------------------------------------------------------
    // Remove several items at once.
    // This function is needed because deleting one element might invalidate
    // other indices to be removed, so sequential one-by-one deletion is not
    // safe.
    //
    // Returns the new indices of previous elements with indices greater than
    // or equal to the final size.
    // The i-th value indicates the new index of previously (finalSize + i)-th
    // component.
    //
    // Invalidates any indices bigger than the final size.
    //
    // If there are any index out of range, or there are repeated indices, the
    // behavior is undefined.
    //
    // The algorithm works as follows:
    //   - Computes the final size
    //   - Removes to-be-deleted items with indices larger than the final size
    //   - Moves the not-to-be-deleted items with indices larger than the final
    //     size to the to-be-deleted items with indices smaller than the final
    //     size.
    //   - Adjust the size to the final size
    template< Size n >
    struct EraseNResult {
        Index newSize;
        std::array< std::optional< Index >, n > newIndices;

        template< typename... IndexT >
        auto track(IndexT... oldIndices) const {
            using IT = std::array< Index, sizeof...(IndexT) >;
            IT oi = {{ oldIndices ... }};
            IT res;
            for(Index i = 0; i < oi.size(); ++i) {
                if(oi[i] < newSize) res[i] = oi[i];
                else res[i] = *newIndices[oi[i] - newSize]; // In this case the element must have value
            }
            return res;
        }
    };
    template< std::size_t n, typename Retargeter > // The retargeter must implement operator()(from, to)
    auto erase(const std::array< Index, n >& indices, Retargeter&& r) {
        using namespace std;
        EraseNResult<n> res {};

        // isDeleted[i]: whether value_[finalSize + i] should be deleted. initialized to false
        std::array< bool, n > isDeleted {};
        const Size currentSize = value_.size();
        const Size finalSize = currentSize - n;

        // Mark to-be-deleted items with bigger indices as deleted
        for(Index i = 0; i < n; ++i) {
            if(indices[i] >= finalSize) {
                isDeleted[indices[i] - finalSize] = true;
            }
        }
        // Move the not-to-be-deleted items with bigger indices to the to-be-deleted items with small indices
        for(Index indAfterFinal = 0, i = 0; indAfterFinal < n; ++indAfterFinal) {
            if(!isDeleted[indAfterFinal]) {
                while(i < n && indices[i] >= finalSize) ++i; // Find (including current i) the next i with small index
                if(i < n) {
                    // Found. This should always be satisfied.
                    value_[indices[i]] = std::move(value_[finalSize + indAfterFinal]);
                    r(finalSize + indAfterFinal, indices[i]);
                    res.newIndices[indAfterFinal] = indices[i];
                }
                ++i;
            }
        }

        // Remove garbage
        value_.resize(finalSize);

        res.newSize = finalSize;
        return res;
    }

    Size size() const noexcept { return value_.size(); }

    T&       operator[](Index index)       { return value_[index]; }
    const T& operator[](Index index) const { return value_[index]; }

    // This function should only be called during initialization/finalization
    std::vector< T >& getValue() { return value_; }

};


struct HalfEdgeMeshConnection {

    enum class PolygonType { triangle, border };

    // Forward declarations
    struct VertexConnection;
    struct HalfEdgeConnection;
    struct EdgeConnection;
    struct TriangleConnection;
    struct BorderConnection;

    // Wrapper around Index type, to enhance safety
    template< typename Element >
    struct IndexType {
        Index index = 0;

        // Modifier
        IndexType& operator++() { ++index; return *this; }

        // Comparison
        bool operator==(IndexType rhs) const { return index == rhs.index; }
        bool operator!=(IndexType rhs) const { return !(*this == rhs); }
        bool operator<(Index     rhs) const { return index < rhs; }
        bool operator<(IndexType rhs) const { return index < rhs.index; }
    };
    using VertexIndex   = IndexType< VertexConnection >;
    using HalfEdgeIndex = IndexType< HalfEdgeConnection >;
    using EdgeIndex     = IndexType< EdgeConnection >;
    using TriangleIndex = IndexType< TriangleConnection >;
    using BorderIndex   = IndexType< BorderConnection >;

    // The elements should be movable.
    struct VertexConnection {
        // Only one HalfEdge targeting the vertex is needed.
        HalfEdgeIndex halfEdgeIndex;

        // Number of neighbors
        Size degree;

        // 0: vertex is inside, 1: vertex is on the border, >=2: pathological
        std::uint_fast8_t numTargetingBorderHalfEdges;
    };
    struct HalfEdgeConnection {

        PolygonType polygonType;
        Index       polygonIndex;
        VertexIndex   targetVertexIndex;
        HalfEdgeIndex oppositeHalfEdgeIndex;
        HalfEdgeIndex nextHalfEdgeIndex;
        HalfEdgeIndex prevHalfEdgeIndex;
        EdgeIndex     edgeIndex;
    };
    struct EdgeConnection {
        // Only one HalfEdge is needed.
        HalfEdgeIndex halfEdgeIndex;

        // 0: edge is inside, 1: edge is on the border, 2: pathological
        std::uint_fast8_t numBorderHalfEdges;
    };
    // A triangle is a closed polygon which has exactly 3 half edges.
    struct TriangleConnection {
        // Only one HalfEdge is needed.
        HalfEdgeIndex halfEdgeIndex;
    };
    // A border is a closed polygon which might be non-planar, and should have
    // more than 2 half edges.
    struct BorderConnection {
        // Only one half edge is needed
        HalfEdgeIndex halfEdgeIndex;
    };
};

// The Attribute class must implement
//   - Type VertexAttribute
//     - void setIndex(size_t)
//   - Type EdgeAttribute
//     - void setIndex(size_t)
//   - Type HalfEdgeAttribute
//     - void setIndex(size_t)
//   - Type TriangleAttribute
//     - void setIndex(size_t)
//   - Type BorderAttribute
//     - void setIndex(size_t)
//   - Type MetaAttribute
//   - Type AttributeInitializerInfo
//   - void init(Mesh, const AttributeInitializerInfo&)
//   - AttributeInitializerInfo extract(Mesh)
template< typename Attribute > class HalfEdgeMesh {
public:

    using AttributeType = Attribute;
    using MeshType = HalfEdgeMesh;

    using VertexAttribute   = typename Attribute::VertexAttribute;
    using EdgeAttribute     = typename Attribute::EdgeAttribute;
    using HalfEdgeAttribute = typename Attribute::HalfEdgeAttribute;
    using TriangleAttribute = typename Attribute::TriangleAttribute;
    using BorderAttribute   = typename Attribute::BorderAttribute;
    using MetaAttribute     = typename Attribute::MetaAttribute;

    using PolygonType   = HalfEdgeMeshConnection::PolygonType;

    using VertexIndex   = HalfEdgeMeshConnection::VertexIndex;
    using HalfEdgeIndex = HalfEdgeMeshConnection::HalfEdgeIndex;
    using EdgeIndex     = HalfEdgeMeshConnection::EdgeIndex;
    using TriangleIndex = HalfEdgeMeshConnection::TriangleIndex;
    using BorderIndex   = HalfEdgeMeshConnection::BorderIndex;

    // The elements should be movable.
    struct Vertex : HalfEdgeMeshConnection::VertexConnection {
        using ConnectionType = HalfEdgeMeshConnection::VertexConnection;

        VertexAttribute attr;
    };
    struct HalfEdge : HalfEdgeMeshConnection::HalfEdgeConnection {
        using ConnectionType = HalfEdgeMeshConnection::HalfEdgeConnection;

        HalfEdgeAttribute attr;
    };
    struct Edge : HalfEdgeMeshConnection::EdgeConnection {
        using ConnectionType = HalfEdgeMeshConnection::EdgeConnection;

        EdgeAttribute attr;
    };
    struct Triangle : HalfEdgeMeshConnection::TriangleConnection {
        using ConnectionType = HalfEdgeMeshConnection::TriangleConnection;

        TriangleAttribute attr;
    };
    struct Border : HalfEdgeMeshConnection::BorderConnection {
        using ConnectionType = HalfEdgeMeshConnection::BorderConnection;

        BorderAttribute attr;
    };


    // Auxiliary static factories for indices.
    static auto vertexIndex  (Index index) { return VertexIndex  {index}; }
    static auto halfEdgeIndex(Index index) { return HalfEdgeIndex{index}; }
    static auto edgeIndex    (Index index) { return EdgeIndex    {index}; }
    static auto triangleIndex(Index index) { return TriangleIndex{index}; }
    static auto borderIndex  (Index index) { return BorderIndex  {index}; }

private:

    DeletableVector<Triangle> triangles_;     // collection of triangles
    DeletableVector<HalfEdge> halfEdges_;     // collection of halfedges
    DeletableVector<Edge>     edges_;         // collection of edges
    DeletableVector<Vertex>   vertices_;      // collection of vertices
    DeletableVector<Border>   borders_;       // collection of borders

    MetaAttribute meta_ {};

    int genus_ = 0; // Genus of the surface. Currently it is not tracked.

    // Element accessor
    template< typename Element, std::enable_if_t<std::is_same<Element, Triangle>::value, void>* = nullptr>
    auto& getElements_() { return triangles_; }
    template< typename Element, std::enable_if_t<std::is_same<Element, HalfEdge>::value, void>* = nullptr>
    auto& getElements_() { return halfEdges_; }
    template< typename Element, std::enable_if_t<std::is_same<Element, Edge>::value, void>* = nullptr>
    auto& getElements_() { return edges_; }
    template< typename Element, std::enable_if_t<std::is_same<Element, Vertex>::value, void>* = nullptr>
    auto& getElements_() { return vertices_; }
    template< typename Element, std::enable_if_t<std::is_same<Element, Border>::value, void>* = nullptr>
    auto& getElements_() { return borders_; }

    auto& element_(TriangleIndex ti ) { return triangles_[ti.index]; }
    auto& element_(HalfEdgeIndex hei) { return halfEdges_[hei.index]; }
    auto& element_(EdgeIndex     ei ) { return edges_    [ei.index]; }
    auto& element_(VertexIndex   vi ) { return vertices_ [vi.index]; }
    auto& element_(BorderIndex   bi ) { return borders_  [bi.index]; }
    // The const version is public

    // Meshwork registration helper
    void registerTriangle_(TriangleIndex ti, HalfEdgeIndex hei0, HalfEdgeIndex hei1, HalfEdgeIndex hei2) {
        element_(ti).halfEdgeIndex = hei0;
        element_(hei0).nextHalfEdgeIndex = hei1;
        element_(hei0).prevHalfEdgeIndex = hei2;
        element_(hei0).polygonIndex = ti.index;
        element_(hei1).nextHalfEdgeIndex = hei2;
        element_(hei1).prevHalfEdgeIndex = hei0;
        element_(hei1).polygonIndex = ti.index;
        element_(hei2).nextHalfEdgeIndex = hei0;
        element_(hei2).prevHalfEdgeIndex = hei1;
        element_(hei2).polygonIndex = ti.index;
    }
    void registerEdge_(EdgeIndex ei, HalfEdgeIndex hei0, HalfEdgeIndex hei1) {
        element_(ei).halfEdgeIndex = hei0;
        element_(ei).numBorderHalfEdges =
            static_cast<std::uint_fast8_t>(element_(hei0).polygonType == PolygonType::border) +
            static_cast<std::uint_fast8_t>(element_(hei1).polygonType == PolygonType::border);
        element_(hei0).oppositeHalfEdgeIndex = hei1;
        element_(hei0).edgeIndex = ei;
        element_(hei1).oppositeHalfEdgeIndex = hei0;
        element_(hei1).edgeIndex = ei;
    }

    template< typename Mod >
    auto newVertex_(Mod&& mod) {
        const VertexIndex index { vertices_.insert() };
        mod.newVertex(*this, index);
        return index;
    }
    template< typename Mod >
    auto newEdge_(Mod&& mod) {
        const EdgeIndex index { edges_.insert() };
        mod.newEdge(*this, index);
        return index;
    }
    template< typename Mod >
    auto newHalfEdge_(Mod&& mod) {
        const HalfEdgeIndex index { halfEdges_.insert() };
        mod.newHalfEdge(*this, index);
        return index;
    }
    template< typename Mod >
    auto newTriangle_(Mod&& mod) {
        const TriangleIndex index { triangles_.insert() };
        mod.newTriangle(*this, index);
        return index;
    }
    template< typename Mod >
    auto newBorder_(Mod&& mod) {
        const BorderIndex index { borders_.insert() };
        mod.newBorder(*this, index);
        return index;
    }

    template< typename Context >
    void retargetElement_(Context& sys, VertexIndex from, VertexIndex to) {
        // Need to update all stored indices/reference/pointer to the vertex.
        forEachHalfEdgeTargetingVertex(to, [&](HalfEdgeIndex hei) {
            element_(hei).targetVertexIndex = to;
        });
        element_(to).attr.setIndex(sys, to.index);
    }
    template< typename Context >
    void retargetElement_(Context& sys, HalfEdgeIndex from, HalfEdgeIndex to) {
        element_(opposite(to)).oppositeHalfEdgeIndex = to;

        switch(element_(to).polygonType) {
        case PolygonType::triangle:
            if(element_(triangle(to)).halfEdgeIndex == from)
                element_(triangle(to)).halfEdgeIndex = to;
            break;
        case PolygonType::border:
            if(borders_[polygon(to)].halfEdgeIndex == from)
                borders_[polygon(to)].halfEdgeIndex = to;
            break;
        }

        if(element_(target(to)).halfEdgeIndex == from)
            element_(target(to)).halfEdgeIndex = to;
        if(element_(edge(to)).halfEdgeIndex == from)
            element_(edge(to)).halfEdgeIndex = to;
        element_(next(to)).prevHalfEdgeIndex = to;
        element_(prev(to)).nextHalfEdgeIndex = to;
        element_(to).attr.setIndex(sys, to.index);
    }
    template< typename Context >
    void retargetElement_(Context& sys, EdgeIndex from, EdgeIndex to) {
        forEachHalfEdgeInEdge(to, [this, to](HalfEdgeIndex hei) {
            element_(hei).edgeIndex = to;
        });
        element_(to).attr.setIndex(sys, to.index);
    }
    template< typename Context >
    void retargetElement_(Context& sys, TriangleIndex from, TriangleIndex to) {
        forEachHalfEdgeInTriangle(to, [this, to](HalfEdgeIndex hei) {
            element_(hei).polygonIndex = to.index;
        });
        element_(to).attr.setIndex(sys, to.index);
    }
    template< typename Context >
    void retargetElement_(Context& sys, BorderIndex from, BorderIndex to) {
        forEachHalfEdgeInPolygon(element(to), [this, to](HalfEdgeIndex hei) {
            element_(hei).polygonIndex = to.index;
        });
        element_(to).attr.setIndex(sys, to.index);
    }

    template< typename Element, typename Mod >
    auto removeElement_(Mod&& mod, Index index) {
        using IT = HalfEdgeMeshConnection::IndexType< typename Element::ConnectionType >;
        mod.removeElement(*this, IT {index});
        return getElements_<Element>().erase(
            index,
            [&, this](Index from, Index to) {
                retargetElement_(mod.context(), IT {from}, IT {to});
            }
        );
    }
    template< typename Element, Size n, typename Mod >
    auto removeElements_(Mod&& mod, const std::array< Index, n >& indices) {
        using IT = HalfEdgeMeshConnection::IndexType< typename Element::ConnectionType >;
        for(Index i : indices) mod.removeElement(*this, IT {i});
        return getElements_<Element>().erase(
            indices,
            [&, this](Index from, Index to) {
                retargetElement_(mod.context(), IT {from}, IT {to});
            }
        );
    }

    template< typename Element, typename Mod >
    void clearElement_(Mod&& mod) {
        using IT = HalfEdgeMeshConnection::IndexType< typename Element::ConnectionType >;
        auto& elements = getElements_< Element >();
        for(IT i {0}; i < elements.size(); ++i)
            mod.removeElement(*this, i);
        elements.getValue().clear();
    }

public:
    template< typename Mod >
    void clear(Mod&& mod) {
        clearElement_<Vertex>  (mod);
        clearElement_<HalfEdge>(mod);
        clearElement_<Edge>    (mod);
        clearElement_<Triangle>(mod);
        clearElement_<Border>  (mod);
    }


    // Constructors
    HalfEdgeMesh() = default;

    struct VertexTriangleInitializer {
        struct Info {
            Size numVertices;
            std::vector< std::array< int, 3 > > triangleVertexIndexList;
            typename Attribute::AttributeInitializerInfo attributeInitializerInfo;

            // Optional information
            std::optional< std::vector< std::vector< int > > > borderVertexIndexList;
        };

        struct PathologicalTopologyError : public std::runtime_error {
            PathologicalTopologyError(const std::string& what_arg) : std::runtime_error("Pathological Topology: " + what_arg) {}
        };

        template< typename AttributeModifier >
        void init(
            AttributeModifier&&      mod,
            HalfEdgeMesh&            mesh,
            Size                     numVertices,
            const std::vector< std::array< int, 3 > >& triangleVertexIndexList,
            const typename Attribute::AttributeInitializerInfo& attributeInitializerInfo
        ) const {
            const auto numTriangles = triangleVertexIndexList.size();
            mesh.vertices_.getValue().resize(numVertices);
            mesh.triangles_.getValue().resize(numTriangles);

            const auto estimatedNumHalfEdges = 3 * numTriangles;  // Might be more than this number with borders.
            mesh.halfEdges_.getValue().reserve(estimatedNumHalfEdges);
            mesh.edges_.getValue().reserve(estimatedNumHalfEdges / 2);

            struct VertexAdditionalInfo {
                bool hasTargetingHalfEdge = false;
                std::vector< HalfEdgeIndex > leavingHalfEdgeIndices;
            };
            std::vector< VertexAdditionalInfo > vai(numVertices);

            struct HalfEdgeAdditionalInfo {
                bool isAtBorder;
                bool oppositeBorderCreated = false;
            };
            std::vector< HalfEdgeAdditionalInfo > hai;
            hai.reserve(estimatedNumHalfEdges);

            // Reset targeting border half edge counter
            for(auto& v : mesh.vertices_) {
                v.numTargetingBorderHalfEdges = 0;
            }

            // Build topological information by inserting new half edges
            // The newly added half edges are always in triangle, but might be at the border,
            // if no opposite half edge is registered.
            for(Index ti = 0; ti < numTriangles; ++ti) {
                const auto& t = triangleVertexIndexList[ti];
                mesh.triangles_[ti].halfEdgeIndex.index = mesh.halfEdges_.size(); // The next inserted halfedge index

                for(Index i = 0; i < 3; ++i) {

                    // Insert a new half edge
                    const HalfEdgeIndex hei { mesh.halfEdges_.insert() };
                    HalfEdge& he = mesh.element_(hei);
                    hai.push_back({true});
                    he.polygonType = PolygonType::triangle;
                    he.polygonIndex = ti;
                    he.targetVertexIndex.index = t[i];
                    he.nextHalfEdgeIndex.index = (i == 2 ? hei.index - 2 : hei.index + 1);
                    he.prevHalfEdgeIndex.index = (i == 0 ? hei.index + 2 : hei.index - 1);

                    ++ mesh.element_(he.targetVertexIndex).numTargetingBorderHalfEdges;

                    // Remembering this edge in the vertices.
                    const int leftVertexIndex = t[i == 0 ? 2 : i - 1];
                    vai[leftVertexIndex].leavingHalfEdgeIndices.push_back(hei);

                    // Search in the target vertex, whether there's an opposite halfedge leaving
                    {
                        const auto findRes = std::find_if(
                            vai[t[i]].leavingHalfEdgeIndices.begin(),
                            vai[t[i]].leavingHalfEdgeIndices.end(),
                            [&mesh, leftVertexIndex](HalfEdgeIndex leavingHalfEdgeIndex) {
                                return leftVertexIndex == mesh.target(leavingHalfEdgeIndex).index;
                            }
                        );
                        if(findRes == vai[t[i]].leavingHalfEdgeIndices.end()) {
                            // opposite not found
                            const auto newEdgeIndex = mesh.edges_.insert();
                            mesh.edges_[newEdgeIndex].halfEdgeIndex = hei;
                            mesh.edges_[newEdgeIndex].numBorderHalfEdges = 1;
                            he.edgeIndex.index = newEdgeIndex;
                        } else {
                            // opposite found
                            hai[hei.index].isAtBorder = false;
                            he.oppositeHalfEdgeIndex = *findRes;
                            he.edgeIndex = mesh.edge(he.oppositeHalfEdgeIndex);

                            hai[he.oppositeHalfEdgeIndex.index].isAtBorder = false;
                            mesh.element_(he.oppositeHalfEdgeIndex).oppositeHalfEdgeIndex = hei;

                            mesh.element_(he.edgeIndex).numBorderHalfEdges = 0;
                            -- mesh.element_(he.targetVertexIndex).numTargetingBorderHalfEdges;
                            -- mesh.element_(mesh.target(he.oppositeHalfEdgeIndex)).numTargetingBorderHalfEdges;
                        }
                    }

                    // TODO: search current vertex and make sure no overlapping half edge exists

                    // Set vertex half edge index
                    if(!vai[t[i]].hasTargetingHalfEdge) {
                        mesh.vertices_[t[i]].halfEdgeIndex = hei;
                        vai[t[i]].hasTargetingHalfEdge = true;
                    }
                } // end loop halfedges
            } // end loop triangles

            // Registering vertex degrees and check for vertex topology
            std::vector< int > unusedVertices;
            std::vector< std::pair< int, std::uint_fast8_t > > pathoVertices;
            for(Index vi = 0; vi < numVertices; ++vi) {
                mesh.vertices_[vi].degree = vai[vi].leavingHalfEdgeIndices.size();
                if(mesh.vertices_[vi].degree == 0) unusedVertices.emplace_back(vi);
                if(mesh.vertices_[vi].numTargetingBorderHalfEdges > 1)
                    pathoVertices.emplace_back(vi, mesh.vertices_[vi].numTargetingBorderHalfEdges);
            }
            if(!unusedVertices.empty() || !pathoVertices.empty()) {
                if(!unusedVertices.empty()) {
                    std::ostringstream oss;
                    oss << "The following vertices are unused:";
                    for(auto x : unusedVertices) oss << ' ' << x;
                    LOG(ERROR) << oss.str();
                }
                if(!pathoVertices.empty()) {
                    std::ostringstream oss;
                    oss << "The following vertices have more than 1 targeting border half edges:";
                    for(const auto& x : pathoVertices) oss << ' ' << x.first << ':' << +x.second;
                    LOG(ERROR) << oss.str();
                }

                throw PathologicalTopologyError("There are unused vertices or vertices on multiple borders");
            }

            // Make borders
            {
                // Method of finding the next half edge inside border
                const auto findNextIn = [&mesh, &hai](HalfEdgeIndex hei_b_cur) -> HalfEdgeIndex {
                    auto hei_in = hei_b_cur;
                    do {
                        hei_in = mesh.prev(mesh.opposite(hei_in));
                    } while(!hai[hei_in.index].isAtBorder);
                    return hei_in;
                };

                // Method of associating a new border half edge with vertices, edges and the border
                const auto registerBorderHalfEdge = [&mesh](HalfEdgeIndex hei_b, HalfEdgeIndex hei_in, BorderIndex bi) {
                    mesh.element_(hei_b).polygonType = PolygonType::border;
                    mesh.element_(hei_b).polygonIndex = bi.index;
                    mesh.element_(hei_b).targetVertexIndex = mesh.target(mesh.prev(hei_in));
                    mesh.registerEdge_(mesh.edge(hei_in), hei_in, hei_b);
                };

                // Method of associating 2 consequtive half edges
                const auto connectBorderHalfEdge = [&mesh](HalfEdgeIndex hei_b_cur, HalfEdgeIndex hei_b_last) {
                    mesh.element_(hei_b_cur).prevHalfEdgeIndex = hei_b_last;
                    mesh.element_(hei_b_last).nextHalfEdgeIndex = hei_b_cur;
                };

                const Size currentNumHalfEdges = mesh.halfEdges_.size();
                for(Index hei = 0; hei < currentNumHalfEdges; ++hei) {
                    if(hai[hei].isAtBorder && ! hai[hei].oppositeBorderCreated) {
                        HalfEdgeIndex hei_in_cur { hei };

                        // Create first border half edge
                        const HalfEdgeIndex hei_b_first { mesh.halfEdges_.insert() };
                        hai[hei_in_cur.index].oppositeBorderCreated = true;
                        auto hei_b_cur = hei_b_first;
                        auto hei_b_last = hei_b_first;

                        // Create a border object
                        const BorderIndex bi { mesh.borders_.insert() };
                        mesh.element_(bi).halfEdgeIndex = hei_b_first;
                        registerBorderHalfEdge(hei_b_first, hei_in_cur, bi);

                        // Sequentially create border half edges
                        while( (hei_in_cur = findNextIn(hei_b_cur)).index != hei ) {
                            hei_b_cur.index = mesh.halfEdges_.insert();
                            hai[hei_in_cur.index].oppositeBorderCreated = true;

                            registerBorderHalfEdge(hei_b_cur, hei_in_cur, bi);
                            connectBorderHalfEdge(hei_b_cur, hei_b_last);

                            hei_b_last = hei_b_cur;
                        }
                        connectBorderHalfEdge(hei_b_first, hei_b_last);

                    }
                } // End looping half edges for making border
            } // End making border

            // Initialize attributes
            mod.init(mesh, attributeInitializerInfo);

        } // End of function void init(...)

        template< typename Context >
        Info extract(const Context& sys, const HalfEdgeMesh& mesh) const {
            Info info;
            info.numVertices = mesh.vertices_.size();
            const auto numTriangles = mesh.triangles_.size();
            info.triangleVertexIndexList.resize(numTriangles);

            for(TriangleIndex ti {}; ti < numTriangles; ++ti) {
                Index i = 0;
                mesh.forEachHalfEdgeInTriangle(ti, [&](HalfEdgeIndex hei) {
                    info.triangleVertexIndexList[ti.index][i++] = mesh.target(hei).index;
                });
            }

            info.attributeInitializerInfo = Attribute::extract(sys, mesh);

            // Extract border lists
            info.borderVertexIndexList.emplace(mesh.numBorders());
            for(Index bi = 0; bi < mesh.numBorders(); ++bi) {
                mesh.forEachHalfEdgeInPolygon(
                    mesh.borders_[bi],
                    [&](HalfEdgeIndex hei) {
                        (*info.borderVertexIndexList)[bi].push_back(mesh.target(hei).index);
                    }
                );
            }

            return info;
        }
    };

    // Initialize the meshwork using triangle vertex index lists. Throws on error.
    template< typename Initializer, typename Mod, typename... Args >
    void init(Mod&& mod, Args&&... args) {
        clear(mod); // Clear all the current topology
        Initializer().init(mod, *this, std::forward<Args>(args)...);
    }
    template< typename Initializer, typename Context >
    auto extract(const Context& sys) const {
        return Initializer().extract(sys, *this);
    }

    bool isClosed()const noexcept { return numBorders() == 0; }

    // Data accessors
    auto numVertices()  const noexcept { return vertices_.size(); }
    auto numHalfEdges() const noexcept { return halfEdges_.size(); }
    auto numEdges()     const noexcept { return edges_.size(); }
    auto numTriangles() const noexcept { return triangles_.size(); }
    auto numBorders()   const noexcept { return borders_.size(); }
    const auto& getTriangles() const { return triangles_; }
    const auto& getHalfEdges() const { return halfEdges_; }
    const auto& getEdges()     const { return edges_; }
    const auto& getVertices()  const { return vertices_; }
    const auto& getBorders()   const { return borders_; }

    const auto& element(TriangleIndex ti ) const { return triangles_[ti.index]; }
    const auto& element(HalfEdgeIndex hei) const { return halfEdges_[hei.index]; }
    const auto& element(EdgeIndex     ei ) const { return edges_    [ei.index]; }
    const auto& element(VertexIndex   vi ) const { return vertices_ [vi.index]; }
    const auto& element(BorderIndex   bi ) const { return borders_  [bi.index]; }

    // Attribute accessor
    template< typename TheIndexType >
    auto      & attribute(TheIndexType i)       { return element_(i).attr; }
    template< typename TheIndexType >
    const auto& attribute(TheIndexType i) const { return element(i).attr; }

    MetaAttribute&       metaAttribute()       noexcept { return meta_; }
    const MetaAttribute& metaAttribute() const noexcept { return meta_; }

    // Meshwork traverse
    auto polygonType(HalfEdgeIndex hei) const { return element(hei).polygonType; }
    auto opposite(HalfEdgeIndex hei) const { return element(hei).oppositeHalfEdgeIndex; }
    auto next(HalfEdgeIndex hei) const { return element(hei).nextHalfEdgeIndex; }
    auto prev(HalfEdgeIndex hei) const { return element(hei).prevHalfEdgeIndex; }
    auto triangle(HalfEdgeIndex hei) const { return TriangleIndex{ polygon(hei) }; }
    auto border(HalfEdgeIndex hei) const { return BorderIndex{ polygon(hei) }; }
    Index polygon(HalfEdgeIndex hei) const { return element(hei).polygonIndex; }
    auto target(HalfEdgeIndex hei) const { return element(hei).targetVertexIndex; }
    auto edge(HalfEdgeIndex hei) const { return element(hei).edgeIndex; }

    auto halfEdge(EdgeIndex ei) const { return element(ei).halfEdgeIndex; }
    auto halfEdge(TriangleIndex ti) const { return element(ti).halfEdgeIndex; }

    Size degree(VertexIndex vi) const { return element(vi).degree; }

    bool isVertexOnBorder(const Vertex& v) const { return v.numTargetingBorderHalfEdges >= 1; }
    bool isVertexOnBorder(VertexIndex vi) const { return isVertexOnBorder(element(vi)); }
    bool isEdgeOnBorder(EdgeIndex ei) const { return element(ei).numBorderHalfEdges >= 1; }
    bool isInTriangle(HalfEdgeIndex hei) const { return polygonType(hei) == PolygonType::triangle; }

    // Mesh neighbor iterators
    template< typename Func > void forEachHalfEdgeTargetingVertex(const Vertex& v, Func&& func) const {
        // Counter-clockwise iterating
        auto hei0 = v.halfEdgeIndex;
        auto hei = hei0;
        do {
            func(hei);
            hei = prev(opposite(hei));
        } while(hei != hei0);
    }
    template< typename Func > void forEachHalfEdgeTargetingVertex(VertexIndex vi, Func&& func) const {
        forEachHalfEdgeTargetingVertex(element(vi), std::forward<Func>(func));
    }
    template< typename Polygon, typename Func >
    void forEachHalfEdgeInPolygon(const Polygon& p, Func&& func) const {
        auto hei0 = p.halfEdgeIndex;
        auto hei = hei0;
        do {
            func(hei);
            hei = next(hei);
        } while(hei != hei0);
    }
    template< typename Polygon, typename Func >
    void forEachHalfEdgeInPolygon(Index pi, Func&& func) const {
        forEachHalfEdgeInPolygon(getElements_< Polygon >()[pi], std::forward<Func>(func));
    }
    template< typename Func > void forEachHalfEdgeInTriangle(TriangleIndex ti, Func&& func) const {
        forEachHalfEdgeInPolygon< Triangle >(element(ti), std::forward<Func>(func));
    }
    template< typename Func > void forEachHalfEdgeInEdge(const Edge& e, Func&& func) const {
        auto hei0 = e.halfEdgeIndex;
        func(hei0);
        func(opposite(hei0));
    }
    template< typename Func > void forEachHalfEdgeInEdge(EdgeIndex ei, Func&& func) const {
        forEachHalfEdgeInEdge(element(ei), std::forward<Func>(func));
    }

    //-------------------------------------------------------------------------
    // The following are basic mesh topology operators
    //-------------------------------------------------------------------------

    // template class AttributeModify should have following functions:
    // - Context& context()
    // - void newVertex(...)
    // - void newEdge(...)
    // - void newHalfEdge(...)
    // - void newTriangle(...)
    // - void newBorder(...)
    // - void removeElement(...)

    // Vertex insertion on edge.
    struct VertexInsertionOnEdge {
        static constexpr int deltaNumVertex = 1;

        struct ChangeInfo {

            VertexIndex viNew;

            // tiOld[0] -> tiNew[0] and tiNew[1]
            // tiOld[1] -> tiNew[2] and tiNew[3]
            // Note:
            //   - New triangles are arranged in the ccw direction
            //   - If the insertion happens on the border, the border-side info
            //     is undefined.
            std::array< TriangleIndex, 2 > tiOld;
            std::array< TriangleIndex, 4 > tiNew;
        };


        template< typename AttributeModify >
        auto operator()(HalfEdgeMesh& mesh, EdgeIndex edgeIndex, AttributeModify&& mod) const {

            ChangeInfo res {};

            // Get index of current elements
            const auto ohei       = mesh.halfEdge(edgeIndex);
            const auto ohei_n     = mesh.next(ohei);
            const auto ohei_p     = mesh.prev(ohei);
            const auto ohei_o     = mesh.opposite(ohei);
            const auto ohei_on    = mesh.next(ohei_o);
            const auto ohei_op    = mesh.prev(ohei_o);

            const auto opi0       = mesh.polygon(ohei);
            const auto opt0       = mesh.polygonType(ohei);
            const bool ist0       = opt0 == PolygonType::triangle;
            const auto opi2       = mesh.polygon(ohei_o);
            const auto opt2       = mesh.polygonType(ohei_o);
            const bool ist2       = opt2 == PolygonType::triangle;

            const auto vi0        = mesh.target(ohei);
            const auto vi2        = mesh.target(ohei_o);

            if(!ist0 && !ist2) {
                LOG(ERROR) << "During vertex insertion on edge, neither side of the edge is a triangle.";
                throw std::runtime_error("Invalid edge condition at vertex insertion.");
            }

            // Create new elements
            const auto vi     = mesh.newVertex_(mod);
            const auto ei2    = mesh.newEdge_(mod); // New edge created by splitting
            const auto hei0_o = mesh.newHalfEdge_(mod); // Targeting new vertex, oppositing ohei
            const auto hei2_o = mesh.newHalfEdge_(mod); // Targeting new vertex, oppositing ohei_o

            mesh.element_(vi).degree = 2;
            mesh.element_(vi).numTargetingBorderHalfEdges = 0;

            if(ist0) {
                const auto oti0   = TriangleIndex { opi0 };
                const auto vi1    = mesh.target(ohei_n);

                const auto ei1    = mesh.newEdge_(mod); // New edge cutting t0
                const auto hei1   = mesh.newHalfEdge_(mod); // Leaving new vertex
                const auto hei1_o = mesh.newHalfEdge_(mod); // Targeting new vertex
                const auto ti1    = mesh.newTriangle_(mod);

                mesh.element_(hei1_o).targetVertexIndex = vi;
                mesh.element_(hei1_o).polygonType = PolygonType::triangle;
                mesh.element_(hei1).targetVertexIndex = vi1;
                mesh.element_(hei1).polygonType = PolygonType::triangle;

                ++mesh.element_(vi).degree;
                ++mesh.element_(vi1).degree;

                mesh.registerTriangle_(oti0, ohei,   ohei_n,  hei1_o);
                mesh.registerTriangle_(ti1,  hei1,   ohei_p,  hei2_o);
                mesh.registerEdge_(ei1, hei1, hei1_o);

                res.tiOld[0] = oti0;
                res.tiNew[0] = oti0;
                res.tiNew[1] = ti1;
            } else {
                mesh.element_(ohei_p).nextHalfEdgeIndex = hei2_o;
                mesh.element_(hei2_o).nextHalfEdgeIndex = ohei;
                mesh.element_(ohei).prevHalfEdgeIndex = hei2_o;
                mesh.element_(hei2_o).prevHalfEdgeIndex = ohei_p;

                mesh.element_(hei2_o).polygonIndex = opi0;

                ++mesh.element_(vi).numTargetingBorderHalfEdges;
            }

            if(ist2) {
                const auto oti2   = TriangleIndex { opi2 };
                const auto vi3    = mesh.target(ohei_on);

                const auto ei3    = mesh.newEdge_(mod); // New edge cutting t2
                const auto hei3   = mesh.newHalfEdge_(mod); // Leaving new vertex
                const auto hei3_o = mesh.newHalfEdge_(mod); // Targeting new vertex
                const auto ti3    = mesh.newTriangle_(mod);

                mesh.element_(hei3_o).targetVertexIndex = vi;
                mesh.element_(hei3_o).polygonType = PolygonType::triangle;
                mesh.element_(hei3).targetVertexIndex = vi3;
                mesh.element_(hei3).polygonType = PolygonType::triangle;

                ++mesh.element_(vi).degree;
                ++mesh.element_(vi3).degree;

                mesh.registerTriangle_(oti2, ohei_o, ohei_on, hei3_o);
                mesh.registerTriangle_(ti3,  hei3,   ohei_op, hei0_o);
                mesh.registerEdge_(ei3, hei3, hei3_o);

                res.tiOld[1] = oti2;
                res.tiNew[2] = oti2;
                res.tiNew[3] = ti3;
            } else {
                mesh.element_(ohei_op).nextHalfEdgeIndex = hei0_o;
                mesh.element_(hei0_o).nextHalfEdgeIndex = ohei_o;
                mesh.element_(ohei_o).prevHalfEdgeIndex = hei0_o;
                mesh.element_(hei0_o).prevHalfEdgeIndex = ohei_op;

                mesh.element_(hei0_o).polygonIndex = opi2;

                ++mesh.element_(vi).numTargetingBorderHalfEdges;
            }

            // Adjust vertex and new half edges
            mesh.element_(vi).halfEdgeIndex = hei0_o;

            mesh.element_(hei0_o).targetVertexIndex = vi;
            mesh.element_(hei2_o).targetVertexIndex = vi;
            mesh.element_(hei0_o).polygonType = opt2;
            mesh.element_(hei2_o).polygonType = opt0;

            // Adjust edge
            mesh.registerEdge_(edgeIndex, ohei,   hei0_o);
            mesh.registerEdge_(ei2,       ohei_o, hei2_o);

            res.viNew = vi;
            return res;
        }
    };
    // Edge collapse
    // Note:
    //   - If the edge is not border, at most one vertex can lie on the border
    struct HalfEdgeCollapse {
        static constexpr int deltaNumVertex = -1;

        struct ChangeInfo {
            VertexIndex viTo;
        };

        // The target of the halfedge ohei will be preserved
        // Notice that halfedge index (not edge index) is used in this function.
        template< typename AttributeModifier >
        auto operator()(HalfEdgeMesh& mesh, HalfEdgeIndex ohei, AttributeModifier&& mod) const {

            // Preconditions should be handled by the caller

            ChangeInfo res {};

            const auto collapseInnerEdge = [&](
                MeshType::EdgeIndex     oei,
                MeshType::HalfEdgeIndex ohei,
                MeshType::HalfEdgeIndex ohei_o
            ) {

                const auto ohei_n  = mesh.next(ohei);
                const auto ohei_p  = mesh.prev(ohei);
                const auto ohei_on = mesh.next(ohei_o);
                const auto ohei_op = mesh.prev(ohei_o);

                const auto ot0 = mesh.triangle(ohei);
                const auto ot1 = mesh.triangle(ohei_o);
                const auto ov0 = mesh.target(ohei); // Will collapse to this vertex
                const auto ov1 = mesh.target(ohei_o); // Will be removed
                const auto oei1 = mesh.edge(ohei_n); // Will collapse to this edge
                const auto oei2 = mesh.edge(ohei_p); // Will be removed
                const auto oei3 = mesh.edge(ohei_op); // Will collapse on this edge
                const auto oei4 = mesh.edge(ohei_on); // Will be removed

                // Retarget all halfedges pointing v1 to v0
                const auto hei_begin = mesh.opposite(ohei_on); // Changed halfedge begin
                const auto hei_end = mesh.opposite(ohei_n);    // Changed halfedge end
                for(auto hei1 = hei_begin; hei1 != ohei_p; hei1 = mesh.opposite(mesh.next(hei1))) {
                    mesh.element_(hei1).targetVertexIndex = ov0;
                }
                mesh.element_(ov0).halfEdgeIndex = hei_begin;
                mesh.element_(mesh.target(ohei_n)).halfEdgeIndex = mesh.opposite(ohei_p);
                mesh.element_(mesh.target(ohei_on)).halfEdgeIndex = mesh.opposite(ohei_op);

                // Collapse edges
                mesh.registerEdge_(oei1, mesh.opposite(ohei_n), mesh.opposite(ohei_p));
                mesh.registerEdge_(oei3, mesh.opposite(ohei_op), mesh.opposite(ohei_on));

                // Adjust vertex degrees
                mesh.element_(ov0).degree += mesh.element_(ov1).degree - 4;
                mesh.element_(ov0).numTargetingBorderHalfEdges += mesh.element_(ov1).numTargetingBorderHalfEdges;
                --mesh.element_(mesh.target(ohei_n)).degree;
                --mesh.element_(mesh.target(ohei_on)).degree;

                // Remove elements
                const auto rmv = mesh.template removeElement_<Vertex>(mod, ov1.index);
                mesh.template removeElements_<Edge, 3>(mod, {oei.index, oei2.index, oei4.index});
                mesh.template removeElements_<HalfEdge, 6>(mod, {
                    ohei.index,   ohei_n.index,  ohei_p.index,
                    ohei_o.index, ohei_on.index, ohei_op.index
                });
                mesh.template removeElements_<Triangle, 2>(mod, {ot0.index, ot1.index});

                // Record change info
                res.viTo.index = rmv.track(ov0.index)[0];
            };

            const auto collapseBorderEdge = [&](
                MeshType::EdgeIndex     oei,
                MeshType::HalfEdgeIndex ohei_t,
                MeshType::HalfEdgeIndex ohei_b
            ) {
                const auto ohei_tn = mesh.next(ohei_t);
                const auto ohei_tp = mesh.prev(ohei_t);
                const auto ohei_bn = mesh.next(ohei_b);
                const auto ohei_bp = mesh.prev(ohei_b);

                const auto ot0 = mesh.triangle(ohei_t);
                const auto ob1 = BorderIndex{ mesh.polygon(ohei_b) };

                // Edge collapse cannot be done, if the border has only 3 edges
                if(mesh.next(ohei_bn) == ohei_bp) {
                    log::error("In edge collapse, border {} has only 3 edges.", mesh.polygon(ohei_b));
                    throw std::runtime_error("Invalid connection in edge collapse.");
                }

                const auto ov0 = mesh.target(ohei_t); // Will collapse to this vertex
                const auto ov1 = mesh.target(ohei_b); // Will be removed
                const auto oei1 = mesh.edge(ohei_tn); // Will collapse to this edge
                const auto oei2 = mesh.edge(ohei_tp); // Will be removed

                // Retarget all halfedges pointing v1 to v0
                const auto hei_begin = mesh.opposite(ohei_bn); // Changed halfedge begin
                for(auto hei1 = hei_begin; hei1 != ohei_tp; hei1 = mesh.opposite(mesh.next(hei1))) {
                    mesh.element_(hei1).targetVertexIndex = ov0;
                }
                mesh.element_(ov0).halfEdgeIndex = mesh.opposite(ohei_tn);
                mesh.element_(mesh.target(ohei_tn)).halfEdgeIndex = mesh.opposite(ohei_tp);

                // Collapse edges
                mesh.registerEdge_(oei1, mesh.opposite(ohei_tn), mesh.opposite(ohei_tp));

                // Connect borders
                mesh.element_(ob1).halfEdgeIndex = ohei_bn;
                mesh.element_(ohei_bn).prevHalfEdgeIndex = ohei_bp;
                mesh.element_(ohei_bp).nextHalfEdgeIndex = ohei_bn;

                // Adjust vertex degrees
                mesh.element_(ov0).degree += mesh.element_(ov1).degree - 3;
                --mesh.element_(mesh.target(ohei_tn)).degree;

                // Remove elements
                const auto rmv = mesh.template removeElement_<Vertex>(mod, ov1.index);
                mesh.template removeElements_<Edge, 2>(mod, {oei.index, oei2.index});
                mesh.template removeElements_<HalfEdge, 4>(mod, {
                    ohei_t.index, ohei_tn.index, ohei_tp.index,
                    ohei_b.index
                });
                mesh.template removeElement_<Triangle>(mod, ot0.index);

                // Record change info
                res.viTo.index = rmv.track(ov0.index)[0];
            };

            // Get index of current elements
            const auto oei     = mesh.edge(ohei);
            const auto ohei_o  = mesh.opposite(ohei);
            const auto ist     = mesh.polygonType(ohei)   == PolygonType::triangle;
            const auto ist_o   = mesh.polygonType(ohei_o) == PolygonType::triangle;

            if(ist) {
                if(ist_o) {
                    collapseInnerEdge(oei, ohei, ohei_o);
                }
                else {
                    collapseBorderEdge(oei, ohei, ohei_o);
                }
            }
            else {
                if(ist_o) {
                    collapseBorderEdge(oei, ohei_o, ohei);
                }
                else {
                    log::error("During edge collapse, neither side of the edge is a triangle.");
                    throw std::runtime_error("Invalid edge condition at edge collapse.");
                }
            }

            return res;
        }

    };
    // Edge flip
    // Note: Border edges cannot be flipped
    struct EdgeFlip {
        static constexpr int deltaNumVertex = 0;

        void operator()(HalfEdgeMesh& mesh, EdgeIndex edgeIndex) const {

            // Preconditions (like topology) should be handled by the caller.

            // Get index of current elements
            const auto ohei = mesh.halfEdge(edgeIndex);
            const auto ohei_n = mesh.next(ohei);
            const auto ohei_p = mesh.prev(ohei);
            const auto ohei_o = mesh.opposite(ohei);
            const auto ohei_on = mesh.next(ohei_o);
            const auto ohei_op = mesh.prev(ohei_o);
            const auto ov0 = mesh.target(ohei);
            const auto ov1 = mesh.target(ohei_n);
            const auto ov2 = mesh.target(ohei_o);
            const auto ov3 = mesh.target(ohei_on);
            const auto ot0 = mesh.triangle(ohei);
            const auto ot1 = mesh.triangle(ohei_o);

            // Retarget vertices
            mesh.element_(ohei).targetVertexIndex = ov1;
            mesh.element_(ohei_o).targetVertexIndex = ov3;
            mesh.element_(ov0).halfEdgeIndex = ohei_op;
            mesh.element_(ov2).halfEdgeIndex = ohei_p;

            --mesh.element_(ov0).degree;
            --mesh.element_(ov2).degree;
            ++mesh.element_(ov1).degree;
            ++mesh.element_(ov3).degree;

            // Remake triangles
            mesh.registerTriangle_(ot0, ohei, ohei_p, ohei_on);
            mesh.registerTriangle_(ot1, ohei_o, ohei_op, ohei_n);

        }

    };

}; // End of class HalfEdgeMesh

} // namespace medyan

#endif

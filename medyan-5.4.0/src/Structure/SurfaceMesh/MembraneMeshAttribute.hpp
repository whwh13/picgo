#ifndef MEDYAN_Structure_SurfaceMesh_MembraneMeshAttribute_hpp
#define MEDYAN_Structure_SurfaceMesh_MembraneMeshAttribute_hpp

#include <algorithm> // max
#include <array>
#include <limits> // numeric_limits
#include <memory> // unique_ptr
#include <stdexcept> // logic_error
#include <tuple>
#include <vector>

#include "MathFunctions.h"
#include "Structure/SurfaceMesh/AdaptiveMeshAttribute.hpp"
#include "Structure/SurfaceMesh/Edge.hpp"
#include "Structure/SurfaceMesh/GeometricMeshAttribute.hpp"
#include "Structure/SurfaceMesh/HalfEdge.hpp"
#include "Structure/SurfaceMesh/SurfaceMesh.hpp"
#include "Structure/SurfaceMesh/Triangle.hpp"
#include "Structure/SurfaceMesh/Vertex.hpp"
#include "Structure/SurfaceMesh/Types.hpp"
#include "SysParams.h"
#include "Util/Io/Log.hpp"
#include "Util/StableVector.hpp"

namespace medyan {

// Forward declarations
class Membrane;



/******************************************************************************
Implements the attributes of the meshwork used by the membrane, mainly
geometric attribute and adaptive mesh attributes.

Geometric attributes are mainly used in force field computations.
Adaptive attributes provides additional attributes mainly used in adaptive mesh
algorithm.

This struct can connect the topological meshwork with the actual system objects.
******************************************************************************/
struct MembraneMeshAttribute {

    using MeshType = HalfEdgeMesh< MembraneMeshAttribute >;

    struct VertexAttribute {
        using CoordinateType      = Vertex::CoordinateType;

        medyan::StableVector<Vertex>::Index vertexSysIndex {};

        GVertex gVertex;
        AdaptiveMeshAttribute::VertexAttribute aVertex;

        Size  cachedDegree;
        Index cachedCoordIndex;

        template< typename Context > CoordinateType&       getCoordinate(Context&       sys) const { return sys.vertices[vertexSysIndex].coord; }
        template< typename Context > const CoordinateType& getCoordinate(const Context& sys) const { return sys.vertices[vertexSysIndex].coord; }
        template< typename Context > Vertex&       vertex(Context&       sys) const { return sys.vertices[vertexSysIndex]; }
        template< typename Context > const Vertex& vertex(const Context& sys) const { return sys.vertices[vertexSysIndex]; }

        const GVertex& getGVertex() const { return gVertex; }
        GVertex&       getGVertex()       { return gVertex; }

        template< typename Context > void setIndex(Context& sys, Index index) {
            vertex(sys).setTopoIndex(index);
        }
    };
    struct EdgeAttribute {
        StableVectorIndex<Edge> edgeSysIndex {};

        GEdge gEdge;
        AdaptiveMeshAttribute::EdgeAttribute aEdge;

        // The index of x coordinate of vertices in the vectorized dof array
        // [target of half edge, opposite target]
        Index                               cachedCoordIndex[2];
        HalfEdgeMeshConnection::PolygonType cachedPolygonType[2];
        // [polygon of half edge, opposite polygon]
        Index                               cachedPolygonIndex[2];

        template< typename Context > Edge&       edge(Context&       sys) const { return sys.edges[edgeSysIndex]; }
        template< typename Context > const Edge& edge(const Context& sys) const { return sys.edges[edgeSysIndex]; }

        const GEdge& getGEdge() const { return gEdge; }
        GEdge&       getGEdge()       { return gEdge; }

        template< typename Context > void setIndex(Context& sys, Index index) {
            edge(sys).setTopoIndex(index);
        }
    };
    struct HalfEdgeAttribute {
        std::unique_ptr< medyan::HalfEdge > halfEdge;

        GHalfEdge gHalfEdge;
        AdaptiveMeshAttribute::HalfEdgeAttribute aHalfEdge;

        // The index of x coordinate of vertices in the vectorized dof array
        // [source, target, next target]
        Index cachedCoordIndex[3];

        const GHalfEdge& getGHalfEdge() const { return gHalfEdge; }
        GHalfEdge&       getGHalfEdge()       { return gHalfEdge; }

        template< typename Context > void setIndex(Context& sys, Index index) {}
    };
    struct TriangleAttribute {
        StableVectorIndex<Triangle> triangleSysIndex;

        GTriangle gTriangle;
        AdaptiveMeshAttribute::TriangleAttribute aTriangle;

        // The index of x coordinate of vertices in the vectorized dof array
        // [target of half edge, next target, prev target]
        Index                                 cachedCoordIndex[3];
        // [half edge, next, prev]
        HalfEdgeMeshConnection::HalfEdgeIndex cachedHalfEdgeIndex[3];
        // [edge of half edge, next edge, prev edge]
        HalfEdgeMeshConnection::EdgeIndex     cachedEdgeIndex[3];

        template< typename Context > Triangle&       triangle(Context&       sys) const { return sys.triangles[triangleSysIndex]; }
        template< typename Context > const Triangle& triangle(const Context& sys) const { return sys.triangles[triangleSysIndex]; }

        const GTriangle& getGTriangle() const { return gTriangle; }
        GTriangle&       getGTriangle()       { return gTriangle; }

        template< typename Context> void setIndex(Context& sys, Index index) {
            triangle(sys).setTopoIndex(index);
        }
    };
    struct BorderAttribute {
        // Whether this border touches a reservoir.
        //
        // If the border is connected to a lipid reservoir, a surface tension
        // would be applied at this border, and the value is supplied as a
        // membrane mechanical parameter.
        // The surface tension must be the same for all borders of a membrane,
        // to make energy minimization possible.
        //
        // If the border is not connected to a lipid reservoir, then it simply
        // represents the free boundary of the membrane.
        bool reservoir = false;

        template< typename Context > void setIndex(Context& sys, Index index) {}
    };
    struct MetaAttribute {
        StableVectorIndex<Membrane> membraneSysIndex {};

        MembraneMeshVertexSystem vertexSystem = MembraneMeshVertexSystem::material;
        bool hasLipidReservoir = false;

        // About mech and chem
        //-----------------------------
        bool isMechParamsSet = false;
        bool isChemParamsSet = false;

        // Chemistry information
        //-----------------------------
        MembraneMeshChemistryInfo chemInfo;

        // Cache related
        //-----------------------------
        bool indexCacheForFFValid = false;

        Size vertexMaxDegree;

        // A vertex has undetermined number of neighbors, so the cache structure needs to be determined at run-time.
        //
        // Notes:
        //   - The outer half edge targets the corresponding "neighbor vertex", and is the "prev" of the targeting half edge.
        //   - The polygon corresponds to the targeting half edge.
        std::vector< Index > cachedVertexTopo;
        Index cachedVertexTopoSize() const { return vertexMaxDegree * 5; }
        Index cachedVertexOffsetNeighborCoord(Index idx) const { return cachedVertexTopoSize() * idx; }
        Index cachedVertexOffsetTargetingHE  (Index idx) const { return cachedVertexTopoSize() * idx + vertexMaxDegree; }
        Index cachedVertexOffsetLeavingHE    (Index idx) const { return cachedVertexTopoSize() * idx + vertexMaxDegree * 2; }
        Index cachedVertexOffsetOuterHE      (Index idx) const { return cachedVertexTopoSize() * idx + vertexMaxDegree * 3; }
        Index cachedVertexOffsetPolygon      (Index idx) const { return cachedVertexTopoSize() * idx + vertexMaxDegree * 4; }
    };

    using CoordinateType = typename VertexAttribute::CoordinateType;

    struct AttributeInitializerInfo {
        std::vector< CoordinateType > vertexCoordinateList;
    };

    // Extraction can be done multiple times without allocating/deallocating
    template< typename Context >
    static auto extract(const Context& sys, const MeshType& mesh) {
        AttributeInitializerInfo info;

        const auto numVertices = mesh.getVertices().size();

        info.vertexCoordinateList.reserve(numVertices);
        for(Index i = 0; i < numVertices; ++i) {
            info.vertexCoordinateList.push_back(
                static_cast<CoordinateType>(sys.vertices[mesh.attribute(HalfEdgeMeshConnection::VertexIndex{i}).vertexSysIndex].coord));
        }

        return info;
    }

};


//-------------------------------------------------------------------------
// Basic mesh operations
//-------------------------------------------------------------------------

template< typename Context, typename ContextFunc >
struct ElementAttributeModifier {
    Context* ps = nullptr;
    ContextFunc sysFunc {};

    ElementAttributeModifier(Context* ps) : ps(ps) {}

    Context& context() const { return *ps; }

    void newVertex(MembraneMeshAttribute::MeshType& mesh, HalfEdgeMeshConnection::VertexIndex i) {
        mesh.attribute(i).vertexSysIndex = sysFunc.template emplaceTrackable<Vertex>(*ps, MembraneMeshAttribute::CoordinateType{}, mesh.metaAttribute().membraneSysIndex, i.index);
    }
    void newEdge(MembraneMeshAttribute::MeshType& mesh, HalfEdgeMeshConnection::EdgeIndex e) {
        mesh.attribute(e).edgeSysIndex = sysFunc.template emplaceTrackable<Edge>(*ps, mesh.metaAttribute().membraneSysIndex, e.index);
    }
    void newHalfEdge(MembraneMeshAttribute::MeshType& mesh, HalfEdgeMeshConnection::HalfEdgeIndex he) {
        mesh.attribute(he).halfEdge = std::make_unique< medyan::HalfEdge >();
    }
    void newTriangle(MembraneMeshAttribute::MeshType& mesh, HalfEdgeMeshConnection::TriangleIndex t) {
        mesh.attribute(t).triangleSysIndex = sysFunc.template emplaceTrackable<Triangle>(*ps, mesh.metaAttribute().membraneSysIndex, t.index);
    }
    void newBorder(MembraneMeshAttribute::MeshType& mesh, HalfEdgeMeshConnection::BorderIndex) {
        // Do nothing
    }

    void removeElement(MembraneMeshAttribute::MeshType& mesh, HalfEdgeMeshConnection::VertexIndex i) {
        sysFunc.template removeTrackable<Vertex>(*ps, mesh.attribute(i).vertexSysIndex);
    }
    void removeElement(MembraneMeshAttribute::MeshType& mesh, HalfEdgeMeshConnection::EdgeIndex i) {
        sysFunc.template removeTrackable<Edge>(*ps, mesh.attribute(i).edgeSysIndex);
    }
    void removeElement(MembraneMeshAttribute::MeshType& mesh, HalfEdgeMeshConnection::HalfEdgeIndex i) {
        // Do nothing
    }
    void removeElement(MembraneMeshAttribute::MeshType& mesh, HalfEdgeMeshConnection::TriangleIndex i) {
        sysFunc.template removeTrackable<Triangle>(*ps, mesh.attribute(i).triangleSysIndex);
    }
    void removeElement(MembraneMeshAttribute::MeshType& mesh, HalfEdgeMeshConnection::BorderIndex i) {
        // Do nothing
    }

    // Mesh attribute initializing and extracting
    // These operations do not follow the RAII idiom.
    // Initialization should be executed only once, as it allocates resources
    void init(MembraneMeshAttribute::MeshType& mesh, const MembraneMeshAttribute::AttributeInitializerInfo& info) {
        const MembraneMeshAttribute::MetaAttribute& meta = mesh.metaAttribute();
        for(Index i = 0; i < mesh.getVertices().size(); ++i) {
            mesh.attribute(HalfEdgeMeshConnection::VertexIndex{i}).vertexSysIndex
                = sysFunc.template emplaceTrackable<Vertex>(*ps, info.vertexCoordinateList[i], meta.membraneSysIndex, i);
        }
        for(Index i = 0; i < mesh.getHalfEdges().size(); ++i) {
            mesh.attribute(HalfEdgeMeshConnection::HalfEdgeIndex{i}).halfEdge
                = std::make_unique< medyan::HalfEdge >();
        }
        for(Index i = 0; i < mesh.getEdges().size(); ++i) {
            mesh.attribute(HalfEdgeMeshConnection::EdgeIndex{i}).edgeSysIndex
                = sysFunc.template emplaceTrackable<Edge>(*ps, meta.membraneSysIndex, i);
        }
        for(Index i = 0; i < mesh.getTriangles().size(); ++i) {
            mesh.attribute(HalfEdgeMeshConnection::TriangleIndex{i}).triangleSysIndex
                = sysFunc.template emplaceTrackable<Triangle>(*ps, meta.membraneSysIndex, i);
        }
    }
};

} // namespace medyan

#endif

// WARNING: THIS FILE IS CURRENTLY UNUSABLE
#error "File is not ready to use."

#ifndef MEDYAN_Structure_SurfaceMesh_MembraneMeshCheck_hpp
#define MEDYAN_Structure_SurfaceMesh_MembraneMeshCheck_hpp

#include <unordered_map>
#include <unordered_set>
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Structure/SurfaceMesh/MeshTriangleQuality.hpp"
#include "Util/Io/Log.hpp"

namespace medyan::membrane_mesh_check {

using MeshType = Membrane::MeshType;

struct MembraneMeshInfoDump {
    void addInfo1Ring(
        std::unordered_set<size_t>& vs, std::unordered_set<size_t>& es,
        const MeshType& mesh, size_t vi
    ) const {
        mesh.forEachHalfEdgeTargetingVertex(vi, [&](size_t hei) {
            const size_t hei_p = mesh.prev(hei);
            const size_t vp = mesh.target(hei_p);
            es.insert(mesh.edge(hei));
            es.insert(mesh.edge(hei_p));
            vs.insert(vp);
        });
        vs.insert(vi);
    }
    void operator()(const MeshType& mesh, size_t vi, size_t ring) const {
        std::unordered_set<size_t> vs, es, cvs;
        cvs.insert(vi);
        for(size_t curRing = 1; curRing <= ring; ++curRing) {
            for(auto i : cvs) addInfo1Ring(vs, es, mesh, i);
            std::unordered_set<size_t> vs_next;
            for (auto i : vs) if (cvs.find(i) == cvs.end()) {
                vs_next.insert(i);
            }
            cvs = std::move(vs_next);
        }

        // Output
        std::unordered_map<size_t, size_t> vri;
        size_t index = 0; // 1 based index
        for (auto i : vs) {
            std::cout << mesh.getVertexAttribute(i).getCoordinate() << std::endl;
            vri[i] = (++index);
        }
        for (auto i : es) {
            const size_t hei = mesh.getEdges()[i].halfEdgeIndex;
            std::cout << vri[mesh.target(hei)] << ' ' << vri[mesh.target(mesh.opposite(hei))] << std::endl;
        }

    }
};
struct MembraneMeshTopologyCheck {
    size_t minDegree;
    size_t maxDegree;
    size_t genus = 0;

    bool operator()(const MeshType& mesh, bool report = false) const {
        bool res = true; // Whether the topology is consistent

        const auto numVertices = mesh.numVertices();
        const auto numTriangles = mesh.numTriangles();
        const auto numEdges = mesh.numEdges();
        const auto numHalfEdges = mesh.numHalfEdges();
        const auto numBorders = mesh.numBorders();

        if(numEdges * 2 != numHalfEdges) {
            res = false;
            if(report) {
                LOG(ERROR) << "Incorrect number of edges (" << numEdges << ") vs half edges (" << numHalfEdges << ").";
            }
        }
        if(!numBorders && (numEdges * 2 != numTriangles * 3 || numHalfEdges != numEdges * 2)) {
            res = false;
            if(report) {
                LOG(ERROR) << "Incorrect number of elements of a closed surface";
            }
        }

        // Check vertices
        for(size_t i = 0; i < numVertices; ++i) {
            // Check degree
            const size_t degree = mesh.degree(i);
            size_t degreeActual = 0;
            mesh.forEachHalfEdgeTargetingVertex(i, [&](size_t hei) {
                ++degreeActual;
            });
            if(degreeActual != degree) {
                res = false;
                if(report)
                    LOG(ERROR) << "Inconsistent degree at vertex " << i << ": " << degreeActual
                        << " (expected " << degree << ")";
            }
            if(degree < minDegree) {
                res = false;
                if(report)
                    LOG(ERROR) << "Vertex " << i << " has too small degree " << degree;
            }
            if(degree > maxDegree) {
                res = false;
                if(report)
                    LOG(ERROR) << "Vertex " << i << " has too large degree " << degree;
            }
            // Check targeting half edge index
            const size_t hei = mesh.getVertices()[i].halfEdgeIndex;
            const size_t vi = mesh.target(hei);
            if(vi != i) {
                res = false;
                if(report) {
                    LOG(ERROR) << "Inconsistent targeting half edge at vertex " << i;
                    LOG(INFO) << "Half edge " << hei << " points to vertex " << vi;
                }
            }
        }

        // Check edges
        for(size_t i = 0; i < numEdges; ++i) {
            // Check half edge index
            const size_t hei = mesh.getEdges()[i].halfEdgeIndex;
            const size_t ei = mesh.edge(hei);
            if(ei != i) {
                res = false;
                if(report) {
                    LOG(ERROR) << "Inconsistent half edge index at edge " << i;
                    LOG(INFO) << "Half edge " << hei << " has edge " << ei;
                }
            }
        }

        // Check triangles
        for(size_t i = 0; i < numTriangles; ++i) {
            // Check half edge index
            const size_t hei = mesh.getTriangles()[i].halfEdgeIndex;
            const size_t ti = mesh.triangle(hei);
            if(ti != i) {
                res = false;
                if(report) {
                    LOG(ERROR) << "Inconsistent half edge index at triangle " << i;
                    LOG(INFO) << "Half edge " << hei << " has triangle " << ti;
                }
            }
        }

        // Check half edges
        for(size_t i = 0; i < numHalfEdges; ++i) {
            // Check opposite
            const size_t hei_o = mesh.opposite(i);
            if(mesh.opposite(hei_o) != i) {
                res = false;
                if(report) {
                    LOG(ERROR) << "Inconsistent opposite half edges: "
                        << i << " -> " << hei_o << " -> " << mesh.opposite(hei_o);
                }
            }

            // Check next/prev
            const size_t next = mesh.next(i);
            const size_t prev = mesh.prev(i);
            if(i != mesh.prev(next)) {
                res = false;
                if(report) {
                    LOG(ERROR) << "Next of half edge " << i << " is " << next
                        << ", but its prev is " << mesh.prev(next);
                }
            }
            if(i != mesh.next(prev))  {
                res = false;
                if(report) {
                    LOG(ERROR) << "Prev of half edge " << i << " is " << prev
                        << ", but its next is " << mesh.next(prev);
                }
            }
        }

        if(!res && report)
            LOG(INFO)
                << "Triangles: " << numTriangles
                << "Vertices: " << numVertices
                << "Edges: " << numEdges
                << "Half edges: " << numHalfEdges;

        return res;
    }

};

struct MembraneMeshDihedralCheck {
    double cosDihedralError;
    double cosDihedralWarning;

    bool operator()(const MeshType& mesh, bool report = false) const {
        bool res = true;

        const size_t numEdges = mesh.getEdges().size();
        // Requires triangle unit normals
        for(size_t i = 0; i < numEdges; ++i) {
            const size_t hei = mesh.getEdges()[i].halfEdgeIndex;
            const size_t hei_o = mesh.opposite(hei);

            using PolygonType = MeshType::PolygonType;
            if(
                mesh.getHalfEdges()[hei].polygonType   == PolygonType::triangle &&
                mesh.getHalfEdges()[hei_o].polygonType == PolygonType::triangle
            ) {
                const size_t t0 = mesh.triangle(hei);
                const size_t t1 = mesh.triangle(mesh.opposite(hei));
                const auto& un0 = mesh.getTriangleAttribute(t0).gTriangle.unitNormal;
                const auto& un1 = mesh.getTriangleAttribute(t1).gTriangle.unitNormal;
                const auto cosDihedral = dot(un0, un1);
                if(cosDihedral < cosDihedralError) {
                    res = false;
                    if(report)
                        LOG(ERROR) << "Dihedral on edge " << i << " is too low: " << cosDihedral;
                } else if(cosDihedral < cosDihedralWarning) {
                    if(report)
                        LOG(WARNING) << "Dihedral on edge " << i << " is low: " << cosDihedral;
                }
            }
        }

        return res;
    }
};

template< TriangleQualityCriteria c >
struct MembraneMeshQualityCheck {
    using TriangleQualityType = TriangleQuality< c >;
    double qualityError;
    double qualityWarning;

    bool operator()(const MeshType& mesh, bool report = false) const {
        using namespace mathfunc;

        bool res = true;

        const size_t numTriangles = mesh.getTriangles().size();
        for(size_t i = 0; i < numTriangles; ++i) {
            const size_t hei = mesh.getTriangles()[i].halfEdgeIndex;
            const size_t v0 = mesh.target(hei);
            const size_t v1 = mesh.target(mesh.next(hei));
            const size_t v2 = mesh.target(mesh.prev(hei));
            const Vec3 c0 (mesh.getVertexAttribute(v0).getCoordinate());
            const Vec3 c1 (mesh.getVertexAttribute(v1).getCoordinate());
            const Vec3 c2 (mesh.getVertexAttribute(v2).getCoordinate());

            const auto q = TriangleQualityType{}(c0, c1, c2);

            if(TriangleQualityType::worse(q, qualityError)) {
                res = false;
                if(report)
                    LOG(ERROR) << "Quality of triangle " << i << " is too low: " << q;
            } else if(TriangleQualityType::worse(q, qualityWarning)) {
                if(report)
                    LOG(WARNING) << "Quality of triangle " << i << " is low: " << q;
            }
        }

        return res;
    }

    bool sizeQuality(const MeshType& mesh, bool report = false) const {
        bool res = true;

        const size_t numTriangles = mesh.getTriangles().size();
        for(size_t i = 0; i < numTriangles; ++i) {
            const size_t hei = mesh.getTriangles()[i].halfEdgeIndex;
            const size_t e0 = mesh.edge(hei);
            const size_t e1 = mesh.edge(mesh.next(hei));
            const size_t e2 = mesh.edge(mesh.prev(hei));
            const auto l0 = mesh.getEdgeAttribute(e0).aEdge.eqLength;
            const auto l1 = mesh.getEdgeAttribute(e1).aEdge.eqLength;
            const auto l2 = mesh.getEdgeAttribute(e2).aEdge.eqLength;

            const auto q = TriangleQualityType{}(l0, l1, l2);

            if(TriangleQualityType::worse(q, qualityError)) {
                res = false;
                if(report)
                    LOG(ERROR) << "Quality of size of triangle " << i << " is too low: " << q;
            } else if(TriangleQualityType::worse(q, qualityWarning)) {
                if(report)
                    LOG(WARNING) << "Quality of size of triangle " << i << " is low: " << q;
            }
        }

        return res;
    }
};

template< TriangleQualityCriteria c >
struct MembraneMeshQualityReport {
    using TriangleQualityType = TriangleQuality< c >;

    void operator()(const MeshType& mesh) const {
        using namespace mathfunc;

        const size_t numEdges = mesh.getEdges().size();
        const size_t numTriangles = mesh.getTriangles().size();

        // Check triangle quality
        double worstQuality = TriangleQualityType::best;
        double avgQuality = 0.0;
        for(size_t i = 0; i < numTriangles; ++i) {
            const size_t hei = mesh.getTriangles()[i].halfEdgeIndex;
            const size_t v0 = mesh.target(hei);
            const size_t v1 = mesh.target(mesh.next(hei));
            const size_t v2 = mesh.target(mesh.prev(hei));
            const Vec3 c0 (mesh.getVertexAttribute(v0).getCoordinate());
            const Vec3 c1 (mesh.getVertexAttribute(v1).getCoordinate());
            const Vec3 c2 (mesh.getVertexAttribute(v2).getCoordinate());

            const auto q = TriangleQualityType{}(c0, c1, c2);

            if(TriangleQualityType::worse(q, worstQuality)) {
                worstQuality = q;
            }
            avgQuality += q;
        }
        avgQuality /= numTriangles;
        LOG(INFO) << "Quality of triangles: Avg " << avgQuality << " Worst " << worstQuality;

        // Check dihedral angle
        double minDotNormal = 1.0;
        double avgDotNormal = 0.0;
        for(size_t i = 0; i < numEdges; ++i) {
            const size_t hei = mesh.getEdges()[i].halfEdgeIndex;
            const size_t hei_o = mesh.opposite(hei);

            using PolygonType = MeshType::PolygonType;
            if(
                mesh.getHalfEdges()[hei].polygonType   == PolygonType::triangle &&
                mesh.getHalfEdges()[hei_o].polygonType == PolygonType::triangle
            ) {
                const size_t t0 = mesh.triangle(hei);
                const size_t t1 = mesh.triangle(mesh.opposite(hei));
                const auto& un0 = mesh.getTriangleAttribute(t0).gTriangle.unitNormal;
                const auto& un1 = mesh.getTriangleAttribute(t1).gTriangle.unitNormal;
                const auto cosDihedral = dot(un0, un1);

                if(cosDihedral < minDotNormal) {
                    minDotNormal = cosDihedral;
                }
                avgDotNormal += cosDihedral;
            }
        }
        avgDotNormal /= numEdges;
        LOG(INFO) << "Cosine of dihedral angles: Avg " << avgDotNormal << " Worst " << minDotNormal;
    }

    void sizeQuality(const MeshType& mesh) const {
        const size_t numTriangles = mesh.getTriangles().size();

        // Check triangle quality
        double worstQuality = TriangleQualityType::best;
        double avgQuality = 0.0;
        for(size_t i = 0; i < numTriangles; ++i) {
            const size_t hei = mesh.getTriangles()[i].halfEdgeIndex;
            const size_t e0 = mesh.edge(hei);
            const size_t e1 = mesh.edge(mesh.next(hei));
            const size_t e2 = mesh.edge(mesh.prev(hei));
            const auto l0 = mesh.getEdgeAttribute(e0).aEdge.eqLength;
            const auto l1 = mesh.getEdgeAttribute(e1).aEdge.eqLength;
            const auto l2 = mesh.getEdgeAttribute(e2).aEdge.eqLength;

            const auto q = TriangleQualityType{}(l0, l1, l2);

            if(TriangleQualityType::worse(q, worstQuality)) {
                worstQuality = q;
            }
            avgQuality += q;
        }
        avgQuality /= numTriangles;
        LOG(INFO) << "EXPECTED quality of triangles: Avg " << avgQuality << " Worst " << worstQuality;
    }

};

} // namespace medyan::membrane_mesh_check

#endif

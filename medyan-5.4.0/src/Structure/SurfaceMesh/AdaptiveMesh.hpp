/*

Adaptive mesh algorithm

Implementation inspired by
"An Adaptive Mesh Algorithm for Evolving Surfaces: Simulations of Drop Breakup and Coalescence" (2001)
by Vittorio Cristini, Jerzy Blawzdziewicz and Michael Loewenberg,
"About surface remeshing" (2000) by Pascal J. Frey,
"Geometric surface mesh optimization" (1998) by Pascal J. Frey and Houman Borouchaki

Performs mesh relaxation and topological transformation to improve mesh size
quality and shape quality, while maintaining the geometrical accuracy.

The algorithm works as follows.

Loop
    Update size measures
    Resample vertex to fulfill size quality by inserting/deleting vertex
    Relocate vertex to optimize triangle quality
Until all criteria are met or maximum iterations reached

*/

#ifndef MEDYAN_AdaptiveMesh_hpp
#define MEDYAN_AdaptiveMesh_hpp

#include <algorithm> // max min
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <vector>

#include "MathFunctions.h"
#include "Structure/SubSystem.h"
#include "Structure/SubSystemFunc.hpp"
#include "Structure/SurfaceMesh/AdaptiveMeshGeometryManager.hpp"
#include "Structure/SurfaceMesh/AdaptiveMeshVertexRelocation.hpp"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Structure/SurfaceMesh/MembraneMeshModifier.hpp"
#include "Structure/SurfaceMesh/MeshTriangleQuality.hpp"

namespace medyan::adaptive_mesh {

// Implementation
//-------------------------------------

template< typename Mesh, TriangleQualityCriteria c > class EdgeFlipManager {
public:
    using TriangleQualityType = TriangleQuality< c >;
    using CoordinateType = typename Mesh::AttributeType::CoordinateType;

    enum class State {
        Success,
        InvalidTopo,
        NonCoplanar,
        BadQuality
    };

private:
    size_t minDegree_;
    size_t maxDegree_;
    double minDotNormal_; // To assess whether triangles are coplanar.

public:

    // Constructor
    // Parameters
    //   - minDegree: minimum number of neighbors of any vertex
    //   - maxDegree: maximum number of neighbors of any vertex
    //   - minDotNormal: minimum dot product required between neighboring triangles. Range (-1, 1)
    EdgeFlipManager(size_t minDegree, size_t maxDegree, double minDotNormal) :
        minDegree_(minDegree),
        maxDegree_(maxDegree),
        minDotNormal_(minDotNormal)
    {}

    // Returns whether the edge is flipped.
    // Requires
    //   - vertex degrees
    //   - triangle unit normal
    State tryFlip(SubSystem& sys, Mesh& mesh, typename Mesh::EdgeIndex ei) const {
        using namespace mathfunc;

        const auto hei = mesh.halfEdge(ei);
        const auto hei_o = mesh.opposite(hei);
        const auto hei_n = mesh.next(hei);
        const auto hei_on = mesh.next(hei_o);

        const auto vi0 = mesh.target(hei);
        const auto vi1 = mesh.target(hei_n);
        const auto vi2 = mesh.target(hei_o);
        const auto vi3 = mesh.target(hei_on);
        // Currently the edge connects v0 and v2.
        // If the edge flips, the connection would be between v1 and v3.

        // Check if topo constraint is satisfied.
        if(mesh.isEdgeOnBorder(ei)) return State::InvalidTopo;

        const auto ti0 = mesh.triangle(hei);
        const auto ti1 = mesh.triangle(hei_o);

        if(
            mesh.degree(vi0) <= minDegree_ ||
            mesh.degree(vi2) <= minDegree_ ||
            mesh.degree(vi1) >= maxDegree_ ||
            mesh.degree(vi3) >= maxDegree_
        ) return State::InvalidTopo;

        // Check if the current triangles are coplanar.
        if(dot(
            mesh.attribute(ti0).gTriangle.unitNormal,
            mesh.attribute(ti1).gTriangle.unitNormal
        ) < minDotNormal_) return State::NonCoplanar;

        // Check if the target triangles are coplanar.
        const CoordinateType c0 (mesh.attribute(vi0).getCoordinate(sys));
        const CoordinateType c1 (mesh.attribute(vi1).getCoordinate(sys));
        const CoordinateType c2 (mesh.attribute(vi2).getCoordinate(sys));
        const CoordinateType c3 (mesh.attribute(vi3).getCoordinate(sys));
        const auto n013 = cross(c1 - c0, c3 - c0);
        const auto mag_n013 = magnitude(n013);
        const auto n231 = cross(c3 - c2, c1 - c2);
        const auto mag_n231 = magnitude(n231);
        if(mag_n013 == 0.0 || mag_n231 == 0.0) return State::BadQuality; // Degenerate case
        if(dot(n013, n231) < minDotNormal_ * mag_n013 * mag_n231) return State::NonCoplanar;

        // Check whether triangle quality can be improved.
        const auto q012 = TriangleQualityType{}(c0, c1, c2);
        const auto q230 = TriangleQualityType{}(c2, c3, c0);
        const auto qBefore = TriangleQualityType::worseOne(q012, q230);
        const auto q013 = TriangleQualityType{}(c0, c1, c3);
        const auto q231 = TriangleQualityType{}(c2, c3, c1);
        const auto qAfter = TriangleQualityType::worseOne(q013, q231);
        if( !TriangleQualityType::better(qAfter, qBefore) ) return State::BadQuality;

        // All checks complete. Do the flip.
        medyan::flipEdge(sys, mesh, ei);

        // Set attributes
        for(auto ti : {ti0, ti1}) {
            medyan::adaptiveComputeTriangleNormal(sys, mesh, ti);
            mesh.forEachHalfEdgeInTriangle(ti, [&](auto nhei) {
                medyan::adaptiveComputeAngle(sys, mesh, nhei);
            });
        }

        for(auto vi : {vi0, vi1, vi2, vi3}) {
            medyan::adaptiveComputeVertexNormal(mesh, vi);
        }

        // Does not change the edge preferrable length

        return State::Success;
    }
};

enum class EdgeSplitVertexInsertionMethod {
    MidPoint,
    AvgCurv     // Mid point snapped to the sphere of curvature, averaged
};
template< EdgeSplitVertexInsertionMethod > struct EdgeSplitVertexInsertion;
template<> struct EdgeSplitVertexInsertion< EdgeSplitVertexInsertionMethod::MidPoint > {
    Index v0, v1;
    template< typename Mesh >
    auto coordinate(const SubSystem& sys, const Mesh& mesh) const {
        using CoordinateType = typename Mesh::AttributeType::CoordinateType;

        const auto c0 = mesh.attribute(v0).getCoordinate(sys);
        const auto c1 = mesh.attribute(v1).getCoordinate(sys);
        return static_cast<CoordinateType>((c0 + c1) * 0.5);
    }
};
template<> struct EdgeSplitVertexInsertion< EdgeSplitVertexInsertionMethod::AvgCurv > {
    static constexpr double maxRadiusDistanceRatio = 1e6;
    static constexpr double maxRadiusDistanceRatio2 = maxRadiusDistanceRatio * maxRadiusDistanceRatio;

    Index v0, v1;

    // Requires
    //   - Vertex unit normal
    template< typename Mesh >
    auto coordinate(const SubSystem& sys, const Mesh& mesh) const {
        using namespace mathfunc;
        using CoordinateType = typename Mesh::AttributeType::CoordinateType;

        const CoordinateType c0 (mesh.attribute(typename Mesh::VertexIndex {v0}).getCoordinate(sys));
        const CoordinateType c1 (mesh.attribute(typename Mesh::VertexIndex {v1}).getCoordinate(sys));
        const auto& un0 = mesh.attribute(typename Mesh::VertexIndex {v0}).aVertex.unitNormal;
        const auto& un1 = mesh.attribute(typename Mesh::VertexIndex {v1}).aVertex.unitNormal;

        const auto r = c1 - c0;
        const auto mag2_r = magnitude2(r);

        // Compute radius of curvature (for both vertices)
        // negative: convex; positive: concave
        const auto r0 = mag2_r / (2 * dot(un0, r));
        const auto r1 = -mag2_r / (2 * dot(un1, r));

        CoordinateType res0, res1;

        if(std::abs(r0 * r0 / mag2_r) > maxRadiusDistanceRatio2) {
            res0 = 0.5 * (c0 + c1);
        } else {
            // Compute vector from center of sphere to mid point
            const auto ro0 = 0.5 * r - r0 * un0;
            const auto mag_ro0 = magnitude(ro0);

            if(mag_ro0 == 0.0) {
                throw std::runtime_error("Unit normal is parallel to edge.");
            }

            res0 = c0 + r0 * un0 + ro0 * (std::abs(r0) / mag_ro0);
        }

        if(std::abs(r1 * r1 / mag2_r) > maxRadiusDistanceRatio2) {
            res1 = 0.5 * (c0 + c1);
        } else {
            // Compute vector from center of sphere to mid point
            const auto ro1 = -0.5 * r - r1 * un1;
            const auto mag_ro1 = magnitude(ro1);

            if(mag_ro1 == 0.0) {
                throw std::runtime_error("Unit normal is parallel to edge.");
            }

            res1 = c1 + r1 * un1 + ro1 * (std::abs(r1) / mag_ro1);
        }

        return static_cast<CoordinateType>(0.5 * (res0 + res1));
    }
};

template<
    typename Mesh,
    TriangleQualityCriteria c,
    EdgeSplitVertexInsertionMethod m
> class EdgeSplitManager {
public:
    using EdgeFlipManagerType          = EdgeFlipManager< Mesh, c >;
    using EdgeSplitVertexInsertionType = EdgeSplitVertexInsertion< m >;
    using CoordinateType               = typename Mesh::AttributeType::CoordinateType;
    enum class State {
        success,
        invalidTopo,
        notLongestEdge
    };

private:
    size_t maxDegree_;

public:

    // Constructor
    EdgeSplitManager(size_t maxDegree) : maxDegree_(maxDegree) {}

    // Returns whether a new vertex is inserted.
    // Requires
    //   - Vertex degree
    State trySplit(SubSystem& sys, Mesh& mesh, typename Mesh::EdgeIndex ei, const EdgeFlipManagerType& efm) const {
        using namespace mathfunc;

        const auto hei = mesh.halfEdge(ei);
        const auto hei_o = mesh.opposite(hei);
        const auto hei_n = mesh.next(hei);
        const auto hei_p = mesh.prev(hei);
        const auto hei_on = mesh.next(hei_o);
        const auto hei_op = mesh.prev(hei_o);

        const auto vi0 = mesh.target(hei);
        const auto vi1 = mesh.target(hei_n);
        const auto vi2 = mesh.target(hei_o);
        const auto vi3 = mesh.target(hei_on);

        const auto ei0 = mesh.edge(hei_n); // v0 - v1
        const auto ei1 = mesh.edge(hei_p); // v1 - v2
        const auto ei2 = mesh.edge(hei_on); // v2 - v3
        const auto ei3 = mesh.edge(hei_op); // v3 - v0

        // Check topology constraints
        // Currently does not support insertion on border edges, but we may also implement that in the future.
        if(mesh.isEdgeOnBorder(ei)) return State::invalidTopo;

        // A new vertex with degree 4 will always be introduced
        if(
            mesh.degree(vi1) >= maxDegree_ ||
            mesh.degree(vi3) >= maxDegree_
        ) return State::invalidTopo;

        // Check whether the current edge is the longest in the triangle
        const CoordinateType c0 (mesh.attribute(vi0).getCoordinate(sys));
        const CoordinateType c1 (mesh.attribute(vi1).getCoordinate(sys));
        const CoordinateType c2 (mesh.attribute(vi2).getCoordinate(sys));
        const CoordinateType c3 (mesh.attribute(vi3).getCoordinate(sys));
        const auto l2_e = distance2(c0, c2);
        const auto l2_01 = distance2(c0, c1);
        const auto l2_12 = distance2(c1, c2);
        const auto l2_23 = distance2(c2, c3);
        const auto l2_30 = distance2(c3, c0);
        if( !(
            l2_e >= l2_01 && l2_e >= l2_12 &&
            l2_e >= l2_23 && l2_e >= l2_30
        )) return State::notLongestEdge;

        // All checks passed. Do the splitting.
        const auto eqLength = mesh.attribute(ei).aEdge.eqLength;
        const auto change = medyan::insertVertexOnEdge<medyan::SubSystemFunc>(
            sys, mesh, ei,
            EdgeSplitVertexInsertionType { vi0.index, vi2.index }.coordinate(sys, mesh)
        );

        // Update local geometries for adaptive remeshing algorithm
        mesh.forEachHalfEdgeTargetingVertex(change.viNew, [&](auto nhei) {
            const auto nti = mesh.triangle(nhei);
            const auto nei = mesh.edge(nhei);

            medyan::adaptiveComputeTriangleNormal(sys, mesh, nti);
            mesh.forEachHalfEdgeInTriangle(nti, [&](auto nnhei) {
                medyan::adaptiveComputeAngle(sys, mesh, nnhei);
            });

            // Set preferrable length of edges to be the same as before
            mesh.attribute(nei).aEdge.eqLength = eqLength;
        });

        medyan::adaptiveComputeVertexNormal(mesh, change.viNew);
        mesh.forEachHalfEdgeTargetingVertex(change.viNew, [&](auto nhei) {
            const auto nvi = mesh.target(mesh.opposite(nhei));

            medyan::adaptiveComputeVertexNormal(mesh, nvi);
        });

        // Propose edge flipping on surrounding quad edges
        efm.tryFlip(sys, mesh, ei0);
        efm.tryFlip(sys, mesh, ei1);
        efm.tryFlip(sys, mesh, ei2);
        efm.tryFlip(sys, mesh, ei3);

        return State::success;

    }
};

template< typename Mesh, TriangleQualityCriteria c > class EdgeCollapseManager {
public:
    using TriangleQualityType = TriangleQuality< c >;
    using CoordinateType      = typename Mesh::AttributeType::CoordinateType;

    enum class State {
        success,
        invalidTopo,
        notSuitable,
        badQuality
    };

private:
    size_t minDegree_;
    size_t maxDegree_;
    double minQualityImprovement_; // If smaller than 1, then some degradation is allowed.
    double minDotNormal_; // Coplanarness requirement after collapse

    struct PrequalifyResult_ {
        bool   suitable = true; // If this is false, then other values might be undefined.
        double qualityImproved; // The improvement of the worst quality after collapsing.
    };
    // Prequalify the collapse, when the edge is not on the border
    // hei is the direction of collapsing (source gets removed, and target is preserved)
    auto prequalify_(const SubSystem& sys, const Mesh& mesh, typename Mesh::HalfEdgeIndex hei) const {
        using namespace mathfunc;

        PrequalifyResult_ res;

        const auto hei_o = mesh.opposite(hei); // targeting vi1
        const auto hei_n = mesh.next(hei);
        const auto hei_ono = mesh.opposite(mesh.next(hei_o)); // targeting vi1
        const auto hei_opo = mesh.opposite(mesh.prev(hei_o));

        const auto vi0 = mesh.target(hei); // preserved
        const auto vi1 = mesh.target(hei_o); // to be removed

        if(mesh.isVertexOnBorder(vi1)) {
            // Removing border vertex would change the shape of the border.
            res.suitable = false;
            return res;
        }
        if(mesh.attribute(vi1).vertex(sys).getAttachmentRefCount() > 0) {
            // Removing a vertex with attachment is not allowed.
            res.suitable = false;
            return res;
        }

        const CoordinateType c0 (mesh.attribute(vi0).getCoordinate(sys));
        const CoordinateType c1 (mesh.attribute(vi1).getCoordinate(sys));

        const auto ti0 = mesh.triangle(hei);
        const auto ti1 = mesh.triangle(hei_o);

        double qBefore = TriangleQualityType::best;
        double qAfter  = TriangleQualityType::best;
        double minCosDihedral = 1.0; // Coplanar case (best)

        {
            auto chei = hei_o; // chei should always target vi1
            std::optional<Vec3> lastUnitNormal;
            if(mesh.polygonType(mesh.opposite(hei_n)) == Mesh::PolygonType::triangle) {
                lastUnitNormal = mesh.attribute(mesh.triangle(mesh.opposite(hei_n))).gTriangle.unitNormal;
            }
            do {
                const auto ti = mesh.triangle(chei); // valid because vi1 is not on the border
                const auto chei_po = mesh.opposite(mesh.prev(chei));
                const auto vn = mesh.target(mesh.next(chei));
                const auto vp = mesh.target(mesh.prev(chei));
                const CoordinateType cn (mesh.attribute(vn).getCoordinate(sys));
                const CoordinateType cp (mesh.attribute(vp).getCoordinate(sys));

                // Triangle quality before
                qBefore = TriangleQualityType::worseOne(
                    TriangleQualityType{}(cp, c1, cn),
                    qBefore
                );

                // Triangle quality after and Dihedral angle after
                if(ti != ti0 && ti != ti1) {
                    // Quality
                    qAfter = TriangleQualityType::worseOne(
                        TriangleQualityType{}(cp, c0, cn),
                        qAfter
                    );

                    // Dihedral angle
                    auto n_0np = cross(cn - c0, cp - c0);
                    const auto mag_n_0np = magnitude(n_0np);
                    if(mag_n_0np == 0.0) {
                        minCosDihedral = -1.0;
                        break;
                    } else {
                        n_0np *= (1.0 / mag_n_0np);
                        // Dihedral angle with last triangle
                        if(lastUnitNormal.has_value()) {
                            minCosDihedral = std::min(
                                minCosDihedral,
                                dot(n_0np, *lastUnitNormal)
                            );
                        }
                        lastUnitNormal = n_0np;
                        // Dihedral angle with outside triangle
                        if(mesh.polygonType(chei_po) == Mesh::PolygonType::triangle) {
                            minCosDihedral = std::min(
                                minCosDihedral,
                                dot(n_0np, mesh.attribute(mesh.triangle(chei_po)).gTriangle.unitNormal)
                            );
                        }
                        // Special dihedral angle
                        if(
                            chei == hei_ono &&
                            mesh.polygonType(hei_opo) == Mesh::PolygonType::triangle
                        ) {
                            minCosDihedral = std::min(
                                minCosDihedral,
                                dot(n_0np, mesh.attribute(mesh.triangle(hei_opo)).gTriangle.unitNormal)
                            );
                        }
                    }
                }

                // Change chei
                chei = mesh.prev(mesh.opposite(chei)); // counter clockwise around vi1
            } while(chei != hei_o);
        }

        if(minCosDihedral < minDotNormal_) {
            res.suitable = false;
            return res;
        }

        res.qualityImproved = TriangleQualityType::improvement(qBefore, qAfter);
        return res;
    }

public:

    // Constructor
    EdgeCollapseManager(size_t minDegree, size_t maxDegree, double minQualityImprovement, double minDotNormal) :
        minDegree_(minDegree),
        maxDegree_(maxDegree),
        minQualityImprovement_(minQualityImprovement),
        minDotNormal_(minDotNormal)
    {}

    // Returns whether the edge is collapsed
    // Requires
    //   - <None>
    State tryCollapse(SubSystem& sys, Mesh& mesh, typename Mesh::EdgeIndex ei) const {
        using namespace mathfunc;

        const auto hei = mesh.halfEdge(ei);
        const auto hei_o = mesh.opposite(hei);
        const auto hei_n = mesh.next(hei);
        const auto hei_on = mesh.next(hei_o);

        const auto vi0 = mesh.target(hei);
        const auto vi1 = mesh.target(hei_n);
        const auto vi2 = mesh.target(hei_o);
        const auto vi3 = mesh.target(hei_on);
        // Currently the edge connects v0 and v2.
        // If the edge collapses, v0 and v2 would become one point.

        // Check topology constraints
        // Currently does not allow collapsing of border edges, but we may also implement that in the future
        if(mesh.isEdgeOnBorder(ei)) return State::invalidTopo;

        if(
            mesh.degree(vi0) + mesh.degree(vi2) - 4 > maxDegree_ ||
            mesh.degree(vi0) + mesh.degree(vi2) - 4 < minDegree_ ||
            (mesh.isVertexOnBorder(vi0) && mesh.isVertexOnBorder(vi2)) ||
            mesh.degree(vi1) <= minDegree_ ||
            mesh.degree(vi3) <= minDegree_
        ) return State::invalidTopo;

        // Check triangle quality constraints
        // Calculate previous triangle qualities around a vertex
        // if v0 is removed
        const auto pr0 = prequalify_(sys, mesh, hei_o);
        // if v2 is removed
        const auto pr2 = prequalify_(sys, mesh, hei);

        // Choose the best result
        const auto worseThan = [](const PrequalifyResult_& pr0, const PrequalifyResult_& pr1) {
            return (
                (!pr0.suitable && pr1.suitable) ||
                (
                    pr0.suitable == pr1.suitable &&
                    pr0.qualityImproved < pr1.qualityImproved
                )
            );
        };
        const auto& [prChosen, heiChosen] = worseThan(pr0, pr2) ?
            std::tuple(pr2, hei) :
            std::tuple(pr0, hei_o);

        if(!prChosen.suitable) return State::notSuitable;
        if(prChosen.qualityImproved < minQualityImprovement_) return State::badQuality;

        // Do the collapse.
        const auto change = collapseHalfEdge<SubSystemFunc>(
            sys, mesh, heiChosen,
            mesh.attribute(mesh.target(heiChosen)).vertex(sys).coord
        );

        // Set attributes
        mesh.forEachHalfEdgeTargetingVertex(change.viTo, [&](auto nhei) {
            if(mesh.isInTriangle(nhei)) {
                const auto nti = mesh.triangle(nhei);
                medyan::adaptiveComputeTriangleNormal(sys, mesh, nti);
                mesh.forEachHalfEdgeInTriangle(nti, [&](auto nnhei) {
                    medyan::adaptiveComputeAngle(sys, mesh, nnhei);
                });
            }
        });

        medyan::adaptiveComputeVertexNormal(mesh, change.viTo);

        mesh.forEachHalfEdgeTargetingVertex(change.viTo, [&mesh](auto nhei) {
            const auto nvi = mesh.target(mesh.opposite(nhei));
            medyan::adaptiveComputeVertexNormal(mesh, nvi);
        });

        // Does not update edge preferred lengths

        return State::success;
    }
};

enum class SizeMeasureCriteria {
    curvature,
    border,
};
template< SizeMeasureCriteria > struct VertexSizeMeasure;
template<> struct VertexSizeMeasure< SizeMeasureCriteria::curvature > {
    double resolution; // size = res * min_radius_curvature
    double upperLimit; // maximum size

    // Requires
    //   - Vertex unit normal
    template< typename Mesh >
    auto vertexSize(SubSystem& sys, Mesh& mesh, typename Mesh::VertexIndex vi) const {
        using CoordinateType = typename Mesh::AttributeType::CoordinateType;

        double minRadiusCurvature = std::numeric_limits<double>::infinity();
        const auto& un = mesh.attribute(vi).aVertex.unitNormal;
        const CoordinateType ci (mesh.attribute(vi).getCoordinate(sys));
        mesh.forEachHalfEdgeTargetingVertex(vi, [&](auto hei) {
            const auto r = mesh.attribute(mesh.target(mesh.opposite(hei))).getCoordinate(sys) - ci;
            minRadiusCurvature = std::min(
                std::abs(0.5 * magnitude2(r) / dot(un, r)),
                minRadiusCurvature
            );
        });
        
        return std::min(resolution * minRadiusCurvature, upperLimit);
    }
};
template<> struct VertexSizeMeasure< SizeMeasureCriteria::border > {
    template< typename Mesh >
    double vertexSize(SubSystem& sys, Mesh& mesh, typename Mesh::VertexIndex vi) const {
        double ret = inf;
        if(mesh.isVertexOnBorder(vi)) {
            const auto& ci = mesh.attribute(vi).vertex(sys).coord;
            // Use the shortest border edge length.
            mesh.forEachHalfEdgeTargetingVertex(vi, [&](auto hei) {
                const auto hei_o = mesh.opposite(hei);
                if(!mesh.isInTriangle(hei) || !mesh.isInTriangle(hei_o)) {
                    const auto vn = mesh.target(hei_o);
                    ret = std::min(
                        ret,
                        distance(ci, mesh.attribute(vn).vertex(sys).coord)
                    );
                }
            });
        }
        return ret;
    }
};

template< SizeMeasureCriteria... > struct VertexSizeMeasureCombined;
template< SizeMeasureCriteria c, SizeMeasureCriteria... cs >
struct VertexSizeMeasureCombined< c, cs... > {
    template< typename Mesh >
    static auto vertexSize(SubSystem& sys, Mesh& mesh, typename Mesh::VertexIndex vi, const VertexSizeMeasure<c>& vsm, const VertexSizeMeasure<cs>&... vsms) {
        return std::min(vsm.vertexSize(sys, mesh, vi), VertexSizeMeasureCombined<cs...>::vertexSize(sys, mesh, vi, vsms...));
    }
};
template< SizeMeasureCriteria c >
struct VertexSizeMeasureCombined< c > {
    template< typename Mesh >
    static auto vertexSize(SubSystem& sys, Mesh& mesh, typename Mesh::VertexIndex vi, const VertexSizeMeasure<c>& vsm) {
        return vsm.vertexSize(sys, mesh, vi);
    }
};

template< typename Mesh > class SizeMeasureManager {
private:

    double _curvRes; // resolution used in radius curvature
    double _maxSize; // Hard upper bound of size
    Size   _diffuseIter; // Diffusion iterations used in gradation control

    template< SizeMeasureCriteria... cs >
    auto _vertexSize(SubSystem& sys, Mesh& mesh, typename Mesh::VertexIndex vi, const VertexSizeMeasure<cs>&... vsms) const {
        return VertexSizeMeasureCombined<cs...>::vertexSize(sys, mesh, vi, vsms...);
    }
    template< SizeMeasureCriteria... cs >
    void _updateVertexSize(SubSystem& sys, Mesh& mesh, const VertexSizeMeasure<cs>&... vsms) const {
        const Size numVertices = mesh.getVertices().size();
        for(Index i = 0; i < numVertices; ++i) {
            typename Mesh::VertexIndex vi {i};
            mesh.attribute(vi).aVertex.size = _vertexSize(sys, mesh, vi, vsms...);
        }
    }

    void _diffuseSize(Mesh& mesh) const {
        const Size numVertices = mesh.numVertices();

        // Initialize with max size
        for(Index i = 0; i < numVertices; ++i) {
            auto& av = mesh.attribute(typename Mesh::VertexIndex{i}).aVertex;
            av.size = std::min(av.size, _maxSize);
        }

        // Diffuse, with D * Delta t = 0.5, and uniformly weighted Laplace operator
        // l_new = l_old / 2 + (sum of neighbor l_old) / (2 * numNeighbors)
        // Currently, diffusion can only reduce local size measure.
        for(Index iter = 0; iter < _diffuseIter; ++iter) {
            for(Index i = 0; i < numVertices; ++i) {
                typename Mesh::VertexIndex vi {i};
                auto& av = mesh.attribute(vi).aVertex;
                const auto deg = mesh.degree(vi);
    
                double sumSizeNeighbor = 0.0;
                mesh.forEachHalfEdgeTargetingVertex(vi, [&](auto hei) {
                    sumSizeNeighbor += mesh.attribute(mesh.target(mesh.opposite(hei))).aVertex.size;
                });

                av.sizeAux = std::min(
                    0.5 * av.size + 0.5 * sumSizeNeighbor / deg,
                    _maxSize
                ); // capped by _maxSize
            }
            for(Index i = 0; i < numVertices; ++i) {
                auto& av = mesh.attribute(typename Mesh::VertexIndex{i}).aVertex;
                av.size = av.sizeAux;
            }
        }
    }

    void _updateEdgeEqLength(Mesh& mesh) const {
        const Size numEdges = mesh.numEdges();
        for(Index i = 0; i < numEdges; ++i) {
            typename Mesh::EdgeIndex ei{i};
            auto& l0 = mesh.attribute(ei).aEdge.eqLength;
            l0 = 0.0;
            mesh.forEachHalfEdgeInEdge(ei, [&](typename Mesh::HalfEdgeIndex hei) {
                l0 += 0.5 * mesh.attribute(mesh.target(hei)).aVertex.size;
            });
        }
    }

public:

    // Constructor
    SizeMeasureManager(double curvRes, double maxSize, Size diffuseIter) :
        _curvRes(curvRes), _maxSize(maxSize), _diffuseIter(diffuseIter) {}

    // Requires
    //   - Unit normal on each vertex
    void computeSizeMeasure(SubSystem& sys, Mesh& mesh) const {
        VertexSizeMeasure< SizeMeasureCriteria::curvature > vsmCurv {_curvRes, _maxSize};
        VertexSizeMeasure< SizeMeasureCriteria::border > vsmBorder;

        // Compute size on each vertex
        _updateVertexSize(sys, mesh, vsmCurv, vsmBorder);

        // Diffuse size on vertices
        _diffuseSize(mesh);

        // Compute preferred length of edges
        _updateEdgeEqLength(mesh);
    }
}; // End SizeMesasureManager


class MembraneMeshAdapter {
public:
    using MeshType            = Membrane::MeshType;
    using GeometryManagerType = GeometryManager< MeshType >;
    using CoordinateType      = MeshType::AttributeType::CoordinateType;

    static constexpr auto optimalVertexLocationMethod    = OptimalVertexLocationMethod::barycenter;
    static constexpr auto triangleQualityCriteria        = TriangleQualityCriteria::radiusRatio;
    static constexpr auto edgeSplitVertexInsertionMethod = EdgeSplitVertexInsertionMethod::AvgCurv;

private:
    SizeMeasureManager< MeshType > _sizeMeasureManager;
    DirectVertexRelocationManager< optimalVertexLocationMethod > _directVertexRelocationManager;

    EdgeFlipManager< MeshType, triangleQualityCriteria > _edgeFlipManager;
    EdgeSplitManager< MeshType, triangleQualityCriteria, edgeSplitVertexInsertionMethod > _edgeSplitManager;
    EdgeCollapseManager< MeshType, triangleQualityCriteria > edgeCollapseManager_;

    size_t _samplingAdjustmentMaxIter; // Maximum number of scans used in sampling.
    size_t _mainLoopSoftMaxIter; // Maximum iterations of the main loop if topology changes can be reduced to 0
    size_t _mainLoopHardMaxIter; // Maximum iterations of the main loop (hard cap)

    void computeSizeMeasures_(SubSystem& sys, MeshType& mesh) const {
        GeometryManagerType::computeAllTriangleNormals(sys, mesh);
        GeometryManagerType::computeAllAngles(sys, mesh);
        GeometryManagerType::computeAllVertexNormals(mesh);
        _sizeMeasureManager.computeSizeMeasure(sys, mesh);
    }
public:
    // Constructor
    MembraneMeshAdapter(const MeshAdapterSettings& param) :
        _sizeMeasureManager(param.curvatureResolution, param.maxSize, param.diffuseIter),
        _directVertexRelocationManager(
            param.relaxationMaxIterRelocation,
            param.relaxationMaxIterTotal
        ),
        _edgeFlipManager(param.minDegree, param.maxDegree, param.edgeFlipMinDotNormal),
        _edgeSplitManager(param.maxDegree),
        edgeCollapseManager_(
            param.minDegree,
            param.maxDegree,
            param.edgeCollapseMinQualityImprovement,
            param.edgeCollapseMinDotNormal
        ),
        _samplingAdjustmentMaxIter(param.samplingAdjustmentMaxIter),
        _mainLoopSoftMaxIter(param.mainLoopSoftMaxIter),
        _mainLoopHardMaxIter(param.mainLoopHardMaxIter)
    {}

    void adapt(SubSystem& sys, MeshType& mesh) const {
        using namespace mathfunc;

        size_t mainLoopIter = 0;
        while(true) {
            // This outer loop does the following in sequence:
            //
            // 1. Update the desired sizes of all the edges, based on the
            //    size criteria.
            // 2. Perform vertex insertion/deletion operations according to the
            //    desired size. The geometry should be updated.
            // 3. Relocate the vertices to local optimum iteratively. The
            //    geometry should be updated.

            if(!mesh.metaAttribute().isMechParamsSet) {
                // Before setting the mech params, we remove some sharp
                // features introduced by the mesh generation algorithm.
                meshSmoothing(sys, mesh, 0.01, 5);
            }
            computeSizeMeasures_(sys, mesh);

            bool sizeMeasureSatisfied = true;

            size_t countTopoModified;
            size_t iter = 0;
            do {
                // The inner loop traverses the meshwork multiple times to
                // check for elements that do not satisfy the size criteria,
                // and try to fix them by vertex insertion/deletion operations.

                countTopoModified = 0;
                for(MeshType::EdgeIndex ei {0}; ei < mesh.numEdges(); /* No increment here */) {
                    const auto hei0 = mesh.halfEdge(ei);
                    const auto v0 = mesh.target(hei0);
                    const auto v1 = mesh.target(mesh.opposite(hei0));

                    const CoordinateType c0 (mesh.attribute(v0).getCoordinate(sys));
                    const CoordinateType c1 (mesh.attribute(v1).getCoordinate(sys));
                    const double length2 = distance2(c0, c1);

                    const double eqLength = mesh.attribute(ei).aEdge.eqLength;
                    const double eqLength2 = eqLength * eqLength;

                    if(length2 > 2 * eqLength2) { // Too long
                        sizeMeasureSatisfied = false;
                        if(_edgeSplitManager.trySplit(sys, mesh, ei, _edgeFlipManager) == decltype(_edgeSplitManager)::State::success) {
                            // Edge splitting happened. Will check edge ei again next round
                            ++countTopoModified;
                        }
                        else
                            ++ei;
                    } else if(2 * length2 < eqLength2) { // Too short
                        sizeMeasureSatisfied = false;
                        if(edgeCollapseManager_.tryCollapse(sys, mesh, ei) == decltype(edgeCollapseManager_)::State::success) {
                            // Edge collapsing happened. The edge at ei will be different next round
                            ++countTopoModified;
                        }
                        else
                            ++ei;
                    } else { // Check passed
                        ++ei;
                    }
                }

                ++iter;
            } while(countTopoModified && iter < _samplingAdjustmentMaxIter); // If any topology was modified, will loop through all edges again.

            if(
                sizeMeasureSatisfied
                || (
                    countTopoModified == 0
                    && mainLoopIter >= _mainLoopSoftMaxIter
                )
                || mainLoopIter >= _mainLoopHardMaxIter
            ) break;

            _directVertexRelocationManager(sys, mesh, _edgeFlipManager);

            ++mainLoopIter;
        } // End loop TopoModifying-Relaxation

    } // End function adapt(...)

};

} // namespace medyan::adaptive_mesh

#endif

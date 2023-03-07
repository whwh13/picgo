#ifndef MEDYAN_Structure_SurfaceMesh_FuncMembraneMech_hpp
#define MEDYAN_Structure_SurfaceMesh_FuncMembraneMech_hpp

#include "Structure/SubSystem.h"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Structure/SurfaceMesh/FuncMembraneGeo.hpp"

namespace medyan {

//-------------------------------------------------------------------------
// Auxiliary functions for mechanical mesh attributes.
//-------------------------------------------------------------------------

// Set equilibrium area for a vertex.
// Precondition:
// - The equilibrium areas of neighboring triangles are computed.
inline void setVertexEqArea(
    SubSystem&                                   sys,
    MembraneMeshAttribute::MeshType&             mesh,
    MembraneMeshAttribute::MeshType::VertexIndex vi
) {
    using MT = MembraneMeshAttribute::MeshType;

    double totalEqAStar = 0;
    mesh.forEachHalfEdgeTargetingVertex(vi, [&](MT::HalfEdgeIndex hei) {
        if(mesh.polygonType(hei) == MT::PolygonType::triangle) {
            const auto ti = mesh.triangle(hei);
            totalEqAStar += mesh.attribute(ti).triangle(sys).mTriangle.eqArea;
        }
    });

    mesh.attribute(vi).vertex(sys).mVertex.eqArea = totalEqAStar / 3;
}




inline void initMechanicParams(SubSystem& sys, Membrane& membrane, const MembraneSetup& memSetup, const MembraneInit& memInit) {
    using MT = Membrane::MeshType;

    auto& mesh = membrane.getMesh();
    auto& mMembrane = membrane.mMembrane;

    // Validation.
    if(membrane.isClosed() && memInit.volumeOffset != 0) {
        log::error("In membrane {} mechanical initialization, the membrane is closed and the volume offset is not zero.", memInit.name);
        throw std::runtime_error("In membrane mechanical initialization, the membrane is closed and the volume offset is not zero.");
    }

    // Initialize MTriangle
    // Also calculate the total area and volume to set the equilibrium area
    // and volume for MMembrane
    double area = 0.0;
    double volume = memInit.volumeOffset;
    for(MT::TriangleIndex ti {0}; ti < mesh.numTriangles(); ++ti) {
        const auto theArea = medyan::area(sys, mesh, ti);

        area += theArea;
        volume += medyan::coneVolume(sys, mesh, ti);

        auto& mTriangle = mesh.attribute(ti).triangle(sys).mTriangle;
        mTriangle.kArea = memSetup.areaElasticity;
        mTriangle.eqArea = theArea * memInit.eqAreaFactor;
    }

    // Initialize MMembrane
    mMembrane.eqArea = area * memInit.eqAreaFactor;
    mMembrane.eqVolume = volume;
    mMembrane.volumeOffset = memInit.volumeOffset;
    mMembrane.kBending = memSetup.bendingElasticity;
    mMembrane.eqCurv = memSetup.eqMeanCurv;

    // Initialize MVertex
    for(MT::VertexIndex vi {0}; vi < mesh.numVertices(); ++vi) {
        auto& mVertex = mesh.attribute(vi).vertex(sys).mVertex;
        setVertexEqArea(sys, mesh, vi);
    }

    mesh.metaAttribute().isMechParamsSet = true;
} // void initMechanicParams(...)

} // namespace medyan

#endif

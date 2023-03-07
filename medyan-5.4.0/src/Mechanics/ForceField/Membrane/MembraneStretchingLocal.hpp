#ifndef MEDYAN_Mechanics_ForceField_Membrane_MembraneStretching_hpp
#define MEDYAN_Mechanics_ForceField_Membrane_MembraneStretching_hpp

#include <array>
#include <limits>
#include <type_traits>
#include <vector>

#include "Mechanics/ForceField/ForceField.h"
#include "Mechanics/ForceField/Membrane/MembraneStretchingImpl.hpp"
#include "Structure/SubSystem.h"
#include "Structure/SurfaceMesh/FuncMembraneGeo.hpp"
#include "Structure/SurfaceMesh/Membrane.hpp"

namespace medyan {

struct MembraneStretchingLocal : public ForceField {

    MembraneStretchingHarmonic impl;

    // (temp) holds the first index of coordinates of each vertex of triangles
    std::vector< std::array< Index, 3 >> vertexSet;

    // (temp) holds the mech params of triangles
    std::vector< double > kArea;
    std::vector< double > eqArea;

    virtual void vectorize(const FFCoordinateStartingIndex& si) override {
        using namespace std;
        using MT = Membrane::MeshType;

        vertexSet.clear();
        vertexSet.reserve(si.ps->triangles.size()); // Might be more than needed
        kArea.clear();
        kArea.reserve(si.ps->triangles.size());
        eqArea.clear();
        eqArea.reserve(si.ps->triangles.size());

        for(auto& m : si.ps->membranes) {
            // In area elasticity for each triangle, the triangles with any
            // vertex touching an reservoir connecting border will not be
            // included.

            const auto& mesh = m.getMesh();
            // Only applies for material surface coordinate system
            if(mesh.metaAttribute().vertexSystem == medyan::MembraneMeshVertexSystem::material) {

                for(const auto& t : mesh.getTriangles()) {
                    const auto vis = medyan::vertexIndices(mesh, t);

                    // Check if the vertices are reservoir touching.
                    const auto onReservoirBorder = [&](MT::VertexIndex vi) {
                        if(!mesh.isVertexOnBorder(vi)) return false;
                        bool ret = false;
                        mesh.forEachHalfEdgeTargetingVertex(vi, [&](MT::HalfEdgeIndex hei) {
                            if(!mesh.isInTriangle(hei)) {
                                if(mesh.attribute(mesh.border(hei)).reservoir) {
                                    ret = true;
                                }
                            }
                        });
                        return ret;
                    };

                    if(
                        !onReservoirBorder(vis[0]) &&
                        !onReservoirBorder(vis[1]) &&
                        !onReservoirBorder(vis[2])
                    ) {
                        vertexSet.push_back({
                            mesh.attribute(vis[0]).cachedCoordIndex,
                            mesh.attribute(vis[1]).cachedCoordIndex,
                            mesh.attribute(vis[2]).cachedCoordIndex,
                        });
                        kArea.push_back(t.attr.triangle(*si.ps).mTriangle.kArea);
                        eqArea.push_back(t.attr.triangle(*si.ps).mTriangle.eqArea);
                    }
                }
            }
        }
    }

    virtual FP computeEnergy(FP* coord) override {
        using namespace std;

        double en = 0;

        for(Index i = 0; i < vertexSet.size(); ++i) {
            const auto& vs = vertexSet[i];
            const auto enTriangle = impl.energy(
                medyan::area(
                    makeRefVec<3>(coord + vs[0]),
                    makeRefVec<3>(coord + vs[1]),
                    makeRefVec<3>(coord + vs[2])
                ),
                kArea[i],
                eqArea[i]
            );

            if(!isfinite(enTriangle)) {
                log::error("In {} energy calculation, triangle with index {} has energy {}.", getName(), i, enTriangle);

                return inffp;
            }
            else {
                en += enTriangle;
            }

        }

        return en;
    }

    virtual void computeForces(FP* coord, FP* force) override {
        using namespace std;

        for(Index i = 0; i < vertexSet.size(); ++i) {
            const auto& vs = vertexSet[i];
            const auto rv0 = makeRefVec<3>(coord + vs[0]);
            const auto rv1 = makeRefVec<3>(coord + vs[1]);
            const auto rv2 = makeRefVec<3>(coord + vs[2]);
            const auto [area, da0, da1, da2] = medyan::areaAndDerivative(rv0, rv1, rv2);

            impl.forces(force + vs[0], area, da0, kArea[i], eqArea[i]);
            impl.forces(force + vs[1], area, da1, kArea[i], eqArea[i]);
            impl.forces(force + vs[2], area, da2, kArea[i], eqArea[i]);

        }

    }

    virtual std::string getName() override { return "MembraneStretchingLocal"; }

};

} // namespace medyan

#endif

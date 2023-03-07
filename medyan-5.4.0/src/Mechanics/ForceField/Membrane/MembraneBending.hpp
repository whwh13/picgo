#ifndef MEDYAN_Mechancis_ForceField_Membrane_Bending_hpp
#define MEDYAN_Mechancis_ForceField_Membrane_Bending_hpp

#include <algorithm> // transform
#include <functional> // plus
#include <type_traits> // is_same
#include <vector>

#include "common.h" // floatingpoint
#include "Mechanics/ForceField/ForceField.h"
#include "Mechanics/ForceField/Membrane/MembraneBendingHelfrich.hpp"
#include "Mechanics/ForceField/Membrane/MembraneBendingTypes.hpp"
#include "Structure/SurfaceMesh/FuncMembraneGeo.hpp"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Structure/SurfaceMesh/Vertex.hpp"

namespace medyan {

/// Represents a Filament bending interaction
template< typename InteractionType >
struct MembraneBending : public ForceField {

    InteractionType impl;
    SubSystem* ps = nullptr;

    // Protein curvature mismatch parameters.
    bool curvatureMismatchEnableCurvatureSensing = true;
    bool curvatureMismatchEnableCurvatureGeneration = true;
    bool diminishBorderCurvatureMismatch = true;
    double molDiffuseDistance = 25;


    virtual std::string getName() override { return "MembraneBending"; }

    virtual void vectorize(const FFCoordinateStartingIndex& si, const SimulConfig& conf) override {
        ps = si.ps;
        curvatureMismatchEnableCurvatureSensing    = conf.mechParams.curvatureMismatchEnableCurvatureSensing;
        curvatureMismatchEnableCurvatureGeneration = conf.mechParams.curvatureMismatchEnableCurvatureGeneration;

        int numVertices = 0;
        // First pass: find number of vertices.
        for(auto& m : ps->membranes) {
            numVertices += m.getMesh().numVertices();
        }

        const int numCMSpecies = ps->proteinCurvatureMismatchParams.size();
        ps->allVertexCurvatureMismatchParams.resize(numCMSpecies, numVertices);
        // Second pass: allocate curvature mismatch structure.
        int vertexIndex = 0;
        for(auto& m : ps->membranes) {
            for(auto& v : m.getMesh().getVertices()) {
                for(int speciesIndex = 0; speciesIndex < numCMSpecies; ++speciesIndex) {
                    auto& param = ps->allVertexCurvatureMismatchParams(speciesIndex, vertexIndex);
                    auto& objVertex = v.attr.vertex(*ps);
                    auto  proteinArea = ps->proteinCurvatureMismatchParams[speciesIndex].area;

                    param.mol = objVertex.cVertex.species.findSpeciesByIndex(speciesIndex)->getN();
                }

                ++vertexIndex;
            }
        }

        // Diffuse mol to smoothen the distribution. Also diminish any border mols.
        // Run diffusion for unit time, with diffusion coefficient of molDiffuseDistance^2.
        // Determine the number of iterations needed to diffuse the molecules.
        double minArea = inf;
        for(auto& m : ps->membranes) {
            auto& mesh = m.getMesh();
            for(auto& v : mesh.getVertices()) {
                if(!mesh.isVertexOnBorder(v)) {
                    minArea = std::min(minArea, v.attr.gVertex.astar / 3);
                }
            }
        }
        // von Neumann stability condition. Using 1/4 instead of 1/2 to avoid numerical issues.
        const int molDiffuseNumIter = std::max<int>(
            3,
            static_cast<int>(std::ceil(4 * molDiffuseDistance * molDiffuseDistance / minArea))
        );
        if(molDiffuseNumIter > 50) {
            log::warn("Too many iterations ({}) for membrane bending diffusion. Consider decreasing molDiffuseDistance.", molDiffuseNumIter);
        }
        for(int cmIndex = 0; cmIndex < numCMSpecies; ++cmIndex) {
            // Temporary storage for average concentration.
            static std::vector<double> avgc;
            avgc.resize(numVertices);

            for(int iter = 0; iter < molDiffuseNumIter; ++iter) {
                for(auto& m : ps->membranes) {
                    auto& mesh = m.getMesh();
                    for(auto& v : mesh.getVertices()) {
                        const auto& objVertex = v.attr.vertex(*ps);
                        const auto& param = ps->allVertexCurvatureMismatchParams(cmIndex, objVertex.loopIndex);

                        const double curc = param.mol / objVertex.mVertex.eqArea;

                        // If this or any neighbor vertex is on border, set the concentration to zero.
                        bool nearBorder = mesh.isVertexOnBorder(v);

                        double lap = 0;
                        mesh.forEachHalfEdgeTargetingVertex(v, [&, this](auto& hei) {
                            const auto hei_o = mesh.opposite(hei);
                            const auto& vni = mesh.target(hei_o);
                            const auto& objvn = mesh.attribute(vni).vertex(*ps);
                            const auto& paramn = ps->allVertexCurvatureMismatchParams(cmIndex, objvn.loopIndex);
                            if(mesh.isVertexOnBorder(vni)) { nearBorder = true; }
                            const auto nc = paramn.mol / objvn.mVertex.eqArea;

                            // Find laplacian operator.
                            double sumCotTheta = 0;
                            if(mesh.isInTriangle(hei)) {
                                const auto hei_n = mesh.next(hei);
                                sumCotTheta += mesh.attribute(hei_n).gHalfEdge.cotTheta;
                            }
                            if(mesh.isInTriangle(hei_o)) {
                                const auto hei_on = mesh.next(hei_o);
                                sumCotTheta += mesh.attribute(hei_on).gHalfEdge.cotTheta;
                            }
                            lap += sumCotTheta * (nc - curc);
                        });
                        lap /= 2 * v.attr.gVertex.astar / 3;

                        avgc[objVertex.loopIndex] = (diminishBorderCurvatureMismatch && nearBorder)
                            ? 0
                            : curc + std::max(0.0, lap) * molDiffuseDistance * molDiffuseDistance / molDiffuseNumIter;
                    }
                }
                // Copy data back to CM params.
                for(auto& m : ps->membranes) {
                    auto& mesh = m.getMesh();
                    for(auto& v : mesh.getVertices()) {
                        auto& objVertex = v.attr.vertex(*ps);
                        auto& param = ps->allVertexCurvatureMismatchParams(cmIndex, objVertex.loopIndex);

                        param.mol = avgc[objVertex.loopIndex] * objVertex.mVertex.eqArea;
                    }
                }
            }
        }
    }

    virtual FP computeEnergy(FP *coord) override {
        using namespace std;

        auto& cmParams = ps->allVertexCurvatureMismatchParams;

        double en = 0;
        int vertexIndex = 0;

        for (auto& m : ps->membranes) {
            const auto &mesh = m.getMesh();

            // Note that eqCurv is not used for quadratic bending.
            const auto kBending = m.mMembrane.kBending;
            const auto eqCurv   = m.mMembrane.eqCurv;

            for (const auto &v : mesh.getVertices()) {
                double enVertex = 0;

                if (v.numTargetingBorderHalfEdges == 0) {

                    const auto area = (v.attr.gVertex.astar) / 3;

                    if constexpr(std::is_same_v< InteractionType, MembraneBendingHelfrich >) {
                        const auto curv = v.attr.gVertex.curv;

                        enVertex += impl.energy(area, curv, kBending, eqCurv);

                        for(int proteinIndex = 0; proteinIndex < cmParams.numProteins; ++proteinIndex) {
                            auto& vertexProteinParams = cmParams(proteinIndex, vertexIndex);
                            const auto& proteinParams = ps->proteinCurvatureMismatchParams[proteinIndex];

                            const auto eachEnergy = impl.energy(proteinParams.area, curv, proteinParams.kBending, proteinParams.eqCurv);
                            // The current curvature mismatch energy is stored as a side effect.
                            // The stored energy is for ONE protein only.
                            if(curvatureMismatchEnableCurvatureSensing) {
                                vertexProteinParams.energy = eachEnergy;
                            }
                            // However, the total energy includes contribution from all proteins of this kind.
                            if(curvatureMismatchEnableCurvatureGeneration) {
                                enVertex += vertexProteinParams.mol * eachEnergy;
                            }
                        }
                    }
                    else {
                        const auto curv2 = v.attr.gVertex.curv2;

                        enVertex += impl.energy(area, curv2, kBending);

                        // Quadratic bending does not support protein curvature mismatch.
                    }
                }
                else {
                    // In this case, v is a vertex on the membrane border.
                    // We do not compute bending energy for such vertices. Equivently, the curvature is always considered zero or the spontaneous curvature (if supported).
                    // However, for protein curvature mismatch, we need to compute the mismatch energy, where we assume the local curvature is always equal to the spontaneous curvature.

                    const auto area = v.attr.gVertex.astar / 3;
                    if constexpr(std::is_same_v< InteractionType, MembraneBendingHelfrich >) {
                        for(int proteinIndex = 0; proteinIndex < cmParams.numProteins; ++proteinIndex) {
                            auto& vertexProteinParams = cmParams(proteinIndex, vertexIndex);
                            const auto& proteinParams = ps->proteinCurvatureMismatchParams[proteinIndex];

                            const auto eachEnergy = impl.energy(proteinParams.area, eqCurv, proteinParams.kBending, proteinParams.eqCurv);
                            // The current curvature mismatch energy is stored as a side effect.
                            // The stored energy is for ONE protein only.
                            if(curvatureMismatchEnableCurvatureSensing) {
                                vertexProteinParams.energy = eachEnergy;
                            }
                        }
                    }
                }

                if (!isfinite(enVertex)) {
                    log::error("In {} energy calculation, a vertex has energy {}", getName(), enVertex);

                    return numeric_limits<double>::infinity();
                }
                else {
                    en += enVertex;
                }

                ++vertexIndex;
            }

        } // End for membrane

        return en;
    }

    virtual void computeForces(FP *coord, FP *force) override {
        using namespace std;
        using MT = Membrane::MeshType;

        auto& cmParams = ps->allVertexCurvatureMismatchParams;

        int vertexIndex = 0;

        for (auto& m : ps->membranes) {

            medyan::assertValidIndexCacheForFF(m.getMesh());

            const auto &mesh = m.getMesh();
            const auto &cvt = mesh.metaAttribute().cachedVertexTopo;

            // Note that eqCurv is not used for quadratic bending.
            const auto kBending = m.mMembrane.kBending;
            const auto eqCurv   = m.mMembrane.eqCurv;

            const Size numVertices = mesh.numVertices();
            for (Index i = 0; i < numVertices; ++i) {
                MT::VertexIndex vi {i};
                if (!mesh.isVertexOnBorder(vi)) {
                    const auto &va = mesh.attribute(vi);

                    const auto area = va.gVertex.astar / 3;
                    const auto curv = va.gVertex.curv;
                    const auto curv2 = va.gVertex.curv2;

                    // for the central vertex
                    // The central vertex is not on the border
                    if constexpr(std::is_same_v< InteractionType, MembraneBendingHelfrich >) {
                        impl.forces(
                            force + va.cachedCoordIndex,
                            area, va.gVertex.dAstar / 3,
                            curv, va.gVertex.dCurv,
                            kBending, eqCurv
                        );

                        // Forces from curvature mismatch model for all surface proteins.
                        if(curvatureMismatchEnableCurvatureGeneration) {
                            for(int proteinIndex = 0; proteinIndex < cmParams.numProteins; ++proteinIndex) {
                                const auto& vertexProteinParams = cmParams(proteinIndex, vertexIndex);
                                const auto& proteinParams = ps->proteinCurvatureMismatchParams[proteinIndex];

                                impl.forcesConstArea(
                                    force + va.cachedCoordIndex,
                                    proteinParams.area,
                                    curv, va.gVertex.dCurv,
                                    proteinParams.kBending, proteinParams.eqCurv,
                                    vertexProteinParams.mol
                                );
                            }
                        }

                    }
                    else {
                        impl.forces(
                            force + va.cachedCoordIndex,
                            area, va.gVertex.dAstar / 3,
                            curv2, va.gVertex.dCurv2,
                            kBending
                        );

                        // Quadratic bending does not support protein curvature mismatch.
                    }

                    // for neighboring vertices
                    for (int n = 0; n < va.cachedDegree; ++n) {
                        // "hei_o" is the half edge index from the center vertex to the neighbor vertex.
                        const MT::HalfEdgeIndex hei_o { cvt[mesh.metaAttribute().cachedVertexOffsetLeavingHE(i) + n] };

                        // "vn" is the neighbor vertex index in the mesh.
                        // "vn_i" is the index of x-coordinate of vertex vn in the degree-of-freedom-array.
                        const auto vn = mesh.target(hei_o);
                        const auto vn_i = cvt[mesh.metaAttribute().cachedVertexOffsetNeighborCoord(i) + n];

                        const auto &dArea = mesh.attribute(hei_o).gHalfEdge.dNeighborAstar / 3;

                        if constexpr(std::is_same_v< InteractionType, MembraneBendingHelfrich >) {
                            const auto& dCurv = mesh.attribute(hei_o).gHalfEdge.dNeighborCurv;
                            impl.forces(
                                force + vn_i,
                                area, dArea, curv, dCurv, kBending, eqCurv
                            );

                            // Forces from curvature mismatch model for all surface proteins.
                            if(curvatureMismatchEnableCurvatureGeneration) {
                                for(int proteinIndex = 0; proteinIndex < cmParams.numProteins; ++proteinIndex) {
                                    const auto& vertexProteinParams = cmParams(proteinIndex, vertexIndex);
                                    const auto& proteinParams = ps->proteinCurvatureMismatchParams[proteinIndex];

                                    impl.forcesConstArea(
                                        force + vn_i,
                                        proteinParams.area, curv, dCurv, proteinParams.kBending, proteinParams.eqCurv,
                                        vertexProteinParams.mol
                                    );
                                }
                            }

                        }
                        else {
                            const auto& dCurv2 = mesh.attribute(hei_o).gHalfEdge.dNeighborCurv2;
                            impl.forces(
                                force + vn_i,
                                area, dArea, curv2, dCurv2, kBending
                            );

                            // Quadratic bending does not support protein curvature mismatch.
                        }
                    }
                }

                ++vertexIndex;
            }

        } // End for membrane

    }
};

} // namespace medyan

#endif

#ifndef MEDYAN_Structure_SurfaceMesh_MembraneMeshModifier_hpp
#define MEDYAN_Structure_SurfaceMesh_MembraneMeshModifier_hpp

// NOTE
// This file uses both the membrane geometry, as well as the membrane chemistry

#include "Chemistry/AdsorptionDesorption.hpp"
#include "Structure/SubSystem.h"
#include "Structure/SurfaceMesh/FuncMembraneChem.hpp"
#include "Structure/SurfaceMesh/FuncMembraneGeo.hpp"
#include "Structure/SurfaceMesh/FuncMembraneMech.hpp"

namespace medyan {

//-------------------------------------------------------------------------
// Mesh modifiers
//-------------------------------------------------------------------------

// The HalfEdgeMesh already provides basic mesh operations which handles
// the mesh element connections.
// In these operations, some attributes of the mesh need to be handled
// specifically. For example, in material coordinate system, the
// equilibrium area of neighboring triangles after a edge flip need to be
// re-distributed.
// However, these operations cannot be simply achieved by injecting
// behavior (such as passing functions) to the mesh operators of the mesh
// connection class (HalfEdgeMesh), because only with full attribute info
// can one decide what to do/access/modify before/after the mesh
// reconnection.

template< typename ContextFunc >
inline auto insertVertexOnEdge(
    SubSystem&                                   sys,
    MembraneMeshAttribute::MeshType&             mesh,
    MembraneMeshAttribute::MeshType::EdgeIndex   ei,
    const MembraneMeshAttribute::CoordinateType& newPos
) {
    using MT = MembraneMeshAttribute::MeshType;

    const auto ohei       = mesh.halfEdge(ei);
    const auto ohei_o     = mesh.opposite(ohei);

    const auto opt0       = mesh.polygonType(ohei);
    const bool ist0       = opt0 == MT::PolygonType::triangle;
    const auto opt2       = mesh.polygonType(ohei_o);
    const bool ist2       = opt2 == MT::PolygonType::triangle;

    const auto vi0        = mesh.target(ohei);
    const auto vi2        = mesh.target(ohei_o);

    // Do the vertex insertion
    //---------------------------------
    const auto change = MT::VertexInsertionOnEdge{}(mesh, ei, ElementAttributeModifier<SubSystem, ContextFunc>(&sys));

    const auto redisEqArea = [](
        double& eqArea1, double a1,
        double& eqArea2, double a2,
        double totalEqArea
    ) {
        // Postcondition: sum of eq areas is equal to the total eq area.
        eqArea1 = totalEqArea * a1 / (a1 + a2);
        eqArea2 = totalEqArea * a2 / (a1 + a2);
    };

    // Set attributes
    //---------------------------------
    mesh.attribute(change.viNew).vertex(sys).coord = newPos;

    // Redistribute properties
    //---------------------------------
    double& v0EqArea = mesh.attribute(vi0).vertex(sys).mVertex.eqArea;
    double& v2EqArea = mesh.attribute(vi2).vertex(sys).mVertex.eqArea;
    const double v0OldEqArea = v0EqArea;
    const double v2OldEqArea = v2EqArea;
    if(mesh.metaAttribute().isMechParamsSet) {
        if(ist0) {
            const double a1 = area(sys, mesh, change.tiNew[0]);
            const double a2 = area(sys, mesh, change.tiNew[1]);

            auto& eqArea1 = mesh.attribute(change.tiNew[0]).triangle(sys).mTriangle.eqArea;
            auto& eqArea2 = mesh.attribute(change.tiNew[1]).triangle(sys).mTriangle.eqArea;
            if(mesh.metaAttribute().hasLipidReservoir) {
                eqArea1 = a1;
                eqArea2 = a2;
            } else {
                redisEqArea(eqArea1, a1, eqArea2, a2, eqArea1);
            }
        }
        if(ist2) {
            const double a1 = area(sys, mesh, change.tiNew[2]);
            const double a2 = area(sys, mesh, change.tiNew[3]);

            auto& eqArea1 = mesh.attribute(change.tiNew[2]).triangle(sys).mTriangle.eqArea;
            auto& eqArea2 = mesh.attribute(change.tiNew[3]).triangle(sys).mTriangle.eqArea;
            if(mesh.metaAttribute().hasLipidReservoir) {
                eqArea1 = a1;
                eqArea2 = a2;
            } else {
                redisEqArea(eqArea1, a1, eqArea2, a2, eqArea1);
            }
        }

        // Reset vertex equilibrium area.
        setVertexEqArea(sys, mesh, vi0);
        setVertexEqArea(sys, mesh, vi2);
        setVertexEqArea(sys, mesh, change.viNew);
    }

    if(mesh.metaAttribute().isChemParamsSet) {
        // Relink species and reactions
        setSpeciesForVertex(sys, mesh, change.viNew, mesh.metaAttribute().chemInfo);
        // Relocate some species previously in v0 and v2, if the mechanical parameters are set. Otherwise, leave the new vertex with zero species.
        if(mesh.metaAttribute().isMechParamsSet) {
            auto& spe0 = mesh.attribute(vi0).vertex(sys).cVertex.species;
            auto& spe2 = mesh.attribute(vi2).vertex(sys).cVertex.species;
            auto& spen = mesh.attribute(change.viNew).vertex(sys).cVertex.species;
            for(int i = 0; i < spe0.size(); ++i) {
                auto& rs0 = spe0.findSpeciesByIndex(i)->getRSpecies();
                auto& rs2 = spe2.findSpeciesByIndex(i)->getRSpecies();
                auto& rsn = spen.findSpeciesByIndex(i)->getRSpecies();
                const Size nold0 = rs0.getTrueN();
                const Size nold2 = rs2.getTrueN();
                const Size nnew0 = static_cast<Size>(nold0 * std::min<double>(1.0, v0EqArea / v0OldEqArea));
                const Size nnew2 = static_cast<Size>(nold2 * std::min<double>(1.0, v2EqArea / v2OldEqArea));
                rs0.setN(nnew0);
                rs2.setN(nnew2);
                rsn.setN(nold0 + nold2 - nnew0 - nnew2);
            }
        }
        setInternalReactionsForVertex(sys, mesh, change.viNew, mesh.metaAttribute().chemInfo);
        mesh.forEachHalfEdgeTargetingVertex(change.viNew, [&](MT::HalfEdgeIndex nhei) {
            setDiffusionForHalfEdge(sys, mesh, nhei, mesh.metaAttribute().chemInfo);
            setDiffusionForHalfEdge(sys, mesh, mesh.opposite(nhei), mesh.metaAttribute().chemInfo);
        });
        imitateAdsorptionDesorptionReactions(mesh.attribute(vi0).vertex(sys).cVertex, mesh.attribute(change.viNew).vertex(sys).cVertex);
    }

    return change;
}

// Collapse along half edge hei. Source is removed, and target is kept.
template< typename ContextFunc >
inline auto collapseHalfEdge(
    SubSystem&                                             sys,
    MembraneMeshAttribute::MeshType&                       mesh,
    MembraneMeshAttribute::MeshType::HalfEdgeIndex         hei,
    std::optional< MembraneMeshAttribute::CoordinateType > newPos
) {
    using namespace std;
    using MT = MembraneMeshAttribute::MeshType;

    const auto ei = mesh.edge(hei);
    const auto hei_o = mesh.opposite(hei);
    const auto vi0 = mesh.target(hei);
    const auto vi1 = mesh.target(hei_o);

    // Record properties
    //---------------------------------

    // Record old vertex equilibrium areas around v1.
    static std::vector<double>          vAround1OldEqArea;
    static std::vector<MT::VertexIndex> vAround1Index;
    if(mesh.metaAttribute().isMechParamsSet) {
        const auto degree = mesh.getVertices()[vi1.index].degree;
        vAround1OldEqArea.resize(degree);
        vAround1Index.resize(degree);
        int i = 0;
        mesh.forEachHalfEdgeTargetingVertex(vi1, [&](MT::HalfEdgeIndex nhei) {
            vAround1Index[i] = mesh.target(mesh.opposite(nhei));
            vAround1OldEqArea[i] = mesh.attribute(vAround1Index[i]).vertex(sys).mVertex.eqArea;
            ++i;
        });
    }

    // Accumulate around two vertices, and subtract triangles on the middle edge
    double oldTotalEqArea = 0.0;
    const auto addOldTotalEqArea = [&](MT::HalfEdgeIndex hei) {
        if(mesh.isInTriangle(hei)) {
            oldTotalEqArea += mesh.attribute(mesh.triangle(hei)).triangle(sys).mTriangle.eqArea;
        }
    };
    const auto subtractOldTotalEqArea = [&](MT::HalfEdgeIndex hei) {
        if(mesh.isInTriangle(hei)) {
            oldTotalEqArea -= mesh.attribute(mesh.triangle(hei)).triangle(sys).mTriangle.eqArea;
        }
    };
    if(mesh.metaAttribute().isMechParamsSet) {
        mesh.forEachHalfEdgeTargetingVertex(vi0, addOldTotalEqArea);
        mesh.forEachHalfEdgeTargetingVertex(vi1, addOldTotalEqArea);
        mesh.forEachHalfEdgeInEdge(ei, subtractOldTotalEqArea);
    }

    // Move species vectors out of vertices.
    SpeciesPtrContainerVector spe0, spe1;
    static std::vector<Size> nold1;
    if(mesh.metaAttribute().isChemParamsSet) {
        spe0 = move(mesh.attribute(vi0).vertex(sys).cVertex.species);
        spe1 = move(mesh.attribute(vi1).vertex(sys).cVertex.species);
        nold1.resize(spe1.size());
        for(int si = 0; si < spe1.size(); ++si) {
            auto& rs1 = spe1.findSpeciesByIndex(si)->getRSpecies();
            nold1[si] = rs1.getTrueN();
        }
        // spe0 will be moved to the change.viTo vertex.
        // spe1 will be abandoned and destroyed.
        // However, spe1 cannot be destroyed right now, because there are still reactions pointing to it.
    }

    // Do the edge collapse
    //---------------------------------
    const auto change = MT::HalfEdgeCollapse{}(mesh, hei, ElementAttributeModifier<SubSystem, ContextFunc>(&sys));

    // Set attributes
    //---------------------------------
    if(newPos.has_value()) mesh.attribute(change.viTo).vertex(sys).coord = *newPos;

    // Redistribute properties
    //---------------------------------
    if(mesh.metaAttribute().isMechParamsSet) {
        double newTotalArea = 0.0;
        mesh.forEachHalfEdgeTargetingVertex(change.viTo, [&](MT::HalfEdgeIndex nhei) {
            if(mesh.isInTriangle(nhei)) {
                newTotalArea += area(sys, mesh, mesh.triangle(nhei));
            }
        });
        mesh.forEachHalfEdgeTargetingVertex(change.viTo, [&](MT::HalfEdgeIndex nhei) {
            if(mesh.isInTriangle(nhei)) {
                auto& eqArea = mesh.attribute(mesh.triangle(nhei)).triangle(sys).mTriangle.eqArea;
                eqArea = mesh.metaAttribute().hasLipidReservoir
                    ? area(sys, mesh, mesh.triangle(nhei))
                    : oldTotalEqArea * area(sys, mesh, mesh.triangle(nhei)) / newTotalArea;
            }
        });

        // Reset vertex equilibrium area.
        mesh.forEachHalfEdgeTargetingVertex(change.viTo, [&](MT::HalfEdgeIndex nhei) {
            setVertexEqArea(sys, mesh, mesh.target(mesh.opposite(nhei)));
        });
        setVertexEqArea(sys, mesh, change.viTo);
    }

    if(mesh.metaAttribute().isChemParamsSet) {
        // Relink species and reactions.
        //------------------------------
        auto& speTo = mesh.attribute(change.viTo).vertex(sys).cVertex.species;
        speTo = move(spe0);
        // Redistribute species from spe1.
        const auto transferFrac = [](SpeciesPtrContainerVector& speTo, Index si, double frac) -> Size {
            auto& rsTo = speTo.findSpeciesByIndex(si)->getRSpecies();
            const Size ntransfer = static_cast<Size>(frac * nold1[si]);
            rsTo.setN(rsTo.getTrueN() + ntransfer);
            return ntransfer;
        };
        if(mesh.metaAttribute().isMechParamsSet) {
            static std::vector<double> vAround1ExcessArea;
            vAround1ExcessArea.resize(vAround1Index.size());
            double totalExcessEqArea = 0;
            for(int i = 0; i < vAround1Index.size(); ++i) {
                vAround1ExcessArea[i] = std::max<double>(0, mesh.attribute(vAround1Index[i]).vertex(sys).mVertex.eqArea - vAround1OldEqArea[i]);
                totalExcessEqArea += vAround1ExcessArea[i];
            }
            for(int si = 0; si < spe1.size(); ++si) {
                Size ntransfer = 0;
                if(totalExcessEqArea > 0) {
                    for(int i = 0; i < vAround1Index.size(); ++i) {
                        // Transfer some species to the vertex, based on non-negative excess area.
                        ntransfer += transferFrac(mesh.attribute(vAround1Index[i]).vertex(sys).cVertex.species, si, vAround1ExcessArea[i] / totalExcessEqArea);
                    }
                }
                // Transfer the rest to the change.viTo vertex.
                if(ntransfer > nold1[si]) {
                    log::error("Not enough species to transfer. {}/{} already transferred.", ntransfer, nold1[si]);
                    throw std::runtime_error("Not enough species to transfer.");
                }
                auto& rsTo = speTo.findSpeciesByIndex(si)->getRSpecies();
                rsTo.setN(rsTo.getTrueN() + (nold1[si] - ntransfer));
            }
        }
        else {
            // Mechanical parameters are not set.
            // Transfer all species to the change.viTo vertex.
            for(int si = 0; si < spe1.size(); ++si) {
                auto& rsTo = speTo.findSpeciesByIndex(si)->getRSpecies();
                rsTo.setN(rsTo.getTrueN() + nold1[si]);
            }
        }
        setInternalReactionsForVertex(sys, mesh, change.viTo, mesh.metaAttribute().chemInfo);
        mesh.forEachHalfEdgeTargetingVertex(change.viTo, [&](MT::HalfEdgeIndex nhei) {
            setDiffusionForHalfEdge(sys, mesh, nhei, mesh.metaAttribute().chemInfo);
            setDiffusionForHalfEdge(sys, mesh, mesh.opposite(nhei), mesh.metaAttribute().chemInfo);
        });
        // Adsorption/desorption reactions can be destroyed safely.

        //-----------------------------
        // Now spe0 is empty.
        // Now spe1 can be safely destroyed.
    }

    return change;
}

inline void flipEdge(
    SubSystem&                                     sys,
    MembraneMeshAttribute::MeshType&               mesh,
    MembraneMeshAttribute::MeshType::EdgeIndex     ei
) {
    using MT = MembraneMeshAttribute::MeshType;

    // Precondition:
    //   - ei is not on the border

    // Record old attributes
    //---------------------------------
    double oldTotalEqArea = 0.0;
    if(mesh.metaAttribute().isMechParamsSet) {
        mesh.forEachHalfEdgeInEdge(ei, [&](MT::HalfEdgeIndex hei) {
            oldTotalEqArea += mesh.attribute(mesh.triangle(hei)).triangle(sys).mTriangle.eqArea;
        });
    }

    // Do edge flip
    //---------------------------------
    MT::EdgeFlip{}(mesh, ei);

    // Redistribute attributes
    //---------------------------------
    if(mesh.metaAttribute().isMechParamsSet) {
        double newTotalArea = 0.0;
        mesh.forEachHalfEdgeInEdge(ei, [&](MT::HalfEdgeIndex hei) {
            newTotalArea += area(sys, mesh, mesh.triangle(hei));
        });
        mesh.forEachHalfEdgeInEdge(ei, [&](MT::HalfEdgeIndex hei) {
            auto& eqArea = mesh.attribute(mesh.triangle(hei)).triangle(sys).mTriangle.eqArea;
            eqArea = mesh.metaAttribute().hasLipidReservoir
                ? area(sys, mesh, mesh.triangle(hei))
                : oldTotalEqArea * area(sys, mesh, mesh.triangle(hei)) / newTotalArea;
        });

        // Reset vertex equilibrium area.
        mesh.forEachHalfEdgeInEdge(ei, [&](MT::HalfEdgeIndex hei) {
            setVertexEqArea(sys, mesh, mesh.target(hei));
            setVertexEqArea(sys, mesh, mesh.target(mesh.next(hei)));
        });
    }

    if(mesh.metaAttribute().isChemParamsSet) {
        // Relink reactions
        mesh.forEachHalfEdgeInEdge(ei, [&](MT::HalfEdgeIndex hei) {
            setDiffusionForHalfEdge(sys, mesh, hei, mesh.metaAttribute().chemInfo);
        });
        // No update required for adsorption/desorption reactions.
    }
}

template< typename VecType, std::enable_if_t< IsVecLike<VecType>::value >* = nullptr >
inline void resetVertexCoordinate(
    SubSystem&                      sys,
    Membrane::MeshType&             mesh,
    Membrane::MeshType::VertexIndex vi,
    const VecType&                  newCoord
) {
    using HI = Membrane::MeshType::HalfEdgeIndex;

    // Record old properties
    double oldTotalEqArea = 0;
    if(mesh.metaAttribute().isMechParamsSet) {
        mesh.forEachHalfEdgeTargetingVertex(vi, [&](HI hei) {
            if(mesh.isInTriangle(hei)) {
                const auto ti = mesh.triangle(hei);
                oldTotalEqArea += mesh.attribute(ti).triangle(sys).mTriangle.eqArea;
            }
        });
    }

    // Set vertex new coordinate
    mesh.attribute(vi).getCoordinate(sys) = newCoord;

    // Set new properties
    if(mesh.metaAttribute().isMechParamsSet) {
        double newTotalArea = 0.0;
        mesh.forEachHalfEdgeTargetingVertex(vi, [&](HI hei) {
            if(mesh.isInTriangle(hei)) {
                newTotalArea += medyan::area(sys, mesh, mesh.triangle(hei));
            }
        });
        mesh.forEachHalfEdgeTargetingVertex(vi, [&](HI hei) {
            if(mesh.isInTriangle(hei)) {
                auto& eqArea = mesh.attribute(mesh.triangle(hei)).triangle(sys).mTriangle.eqArea;
                eqArea = mesh.metaAttribute().hasLipidReservoir
                    ? medyan::area(sys, mesh, mesh.triangle(hei))
                    : oldTotalEqArea * medyan::area(sys, mesh, mesh.triangle(hei)) / newTotalArea;
            }
        });

        // Reset affected vertex equilibrium area.
        setVertexEqArea(sys, mesh, vi);
        mesh.forEachHalfEdgeTargetingVertex(vi, [&](HI hei) {
            setVertexEqArea(sys, mesh, mesh.target(mesh.opposite(hei)));
        });
    }
}


} // namespace medyan

#endif

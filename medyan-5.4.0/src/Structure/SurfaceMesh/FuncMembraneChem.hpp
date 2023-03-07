#ifndef MEDYAN_Structure_SurfaceMesh_MembraneMeshChemistry_hpp
#define MEDYAN_Structure_SurfaceMesh_MembraneMeshChemistry_hpp

#include <algorithm>
#include <vector>

#include "Rand.h"
#include "Structure/SubSystem.h"
#include "Structure/SurfaceMesh/Membrane.hpp"

namespace medyan {

inline void setSpeciesForVertex(
    SubSystem&                                   sys,
    MembraneMeshAttribute::MeshType&             mesh,
    MembraneMeshAttribute::MeshType::VertexIndex vi,
    const MembraneMeshChemistryInfo&             info
) {
    auto& vs = mesh.attribute(vi).vertex(sys).cVertex.species;

    vs.clear();
    for(const auto& name : info.diffusingSpeciesNames) {
        vs.addSpecies< Species >(name, 0, max_ulim, SpeciesType::DIFFUSING, RSpeciesType::REG);
    }
}

// Precondition: no reaction is associated with any of the species.
inline void clearSpeciesForVertex(
    SubSystem&                                   sys,
    MembraneMeshAttribute::MeshType&             mesh,
    MembraneMeshAttribute::MeshType::VertexIndex vi
) {
    auto& vs = mesh.attribute(vi).vertex(sys).cVertex.species;

    vs.clear();
}

// Helper function: get certain species in CVertex from the species indices
inline auto indicesToSpecies(
    CVertex&                  cv,
    const std::vector< int >& iri
) {
    vector< Species* > vs(iri.size());
    for(int i = 0; i < iri.size(); ++i) {
        vs[i] = cv.species.findSpeciesByIndex(iri[i]);
    }
    return vs;
}

// Precondition: Species must have been set in the vertex
inline void setInternalReactionsForVertex(
    SubSystem&                                   sys,
    MembraneMeshAttribute::MeshType&             mesh,
    MembraneMeshAttribute::MeshType::VertexIndex vi,
    const MembraneMeshChemistryInfo&             info
) {
    auto& cv = mesh.attribute(vi).vertex(sys).cVertex;

    cv.reactions.clear();
    for(const auto& ri : info.internalReactions) {
        cv.reactions.push_back(make_unique< ReactionDy >(
            indicesToSpecies(cv, ri.reactantSpeciesIndices),
            indicesToSpecies(cv, ri.productSpeciesIndices),
            ReactionType::REGULAR,
            ri.rateConstant,
            // volume is not set here
            1.0,
            // rate = rateConstant * area^(1 - numReactants)
            1 - static_cast<int>(ri.reactantSpeciesIndices.size())
        ));
    }
}

inline void clearInternalReactionsForVertex(
    SubSystem&                                   sys,
    MembraneMeshAttribute::MeshType&             mesh,
    MembraneMeshAttribute::MeshType::VertexIndex vi
) {
    auto& cv = mesh.attribute(vi).vertex(sys).cVertex;

    cv.reactions.clear();
}

// Preconditions:
//   - Species must have been set in the vertex
//   - Reactions must have been set in the vertex
//   - 1-ring area around the vertex is updated
inline void setInternalReactionRatesForVertex(
    SubSystem&                                   sys,
    MembraneMeshAttribute::MeshType&             mesh,
    MembraneMeshAttribute::MeshType::VertexIndex vi
) {
    using MT = MembraneMeshAttribute::MeshType;

    const auto& va = mesh.attribute(vi);

    for(auto& r : va.vertex(sys).cVertex.reactions) {
        // Using "volume frac" as volume
        r->setVolumeFrac(va.gVertex.astar / 3);
    }
}

// Update adsorption reaction rates due to change of accessible area.
inline void updateAdsorptionReactionRateAtVertex(SubSystem& sys, StableVectorIndex<Vertex> vsi) {
    auto& cv = sys.vertices[vsi].cVertex;
    for(auto& par : cv.adsorptionReactions) {
        par->setRateMulFactor(std::max<FP>(0, cv.accessibleArea), ReactionBase::RateMulFactorType::mechanochemical);
        par->updatePropensity();
    }
}

// Callback for surface diffusion events.
// Will update accessible areas for affected vertices.
struct MembraneDiffusionReactionCallback {
    // Area usage of bound protein.
    FP useArea = 0;
    // Pointer to system.
    SubSystem* ps = nullptr;
    // Vertex indices.
    StableVectorIndex<Vertex> viFrom {};
    StableVectorIndex<Vertex> viTo   {};

    void operator()(ReactionBase* pr) const {
        auto& aaFrom = ps->vertices[viFrom].cVertex.accessibleArea;
        auto& aaTo   = ps->vertices[viTo].cVertex.accessibleArea;
        aaFrom += useArea;
        aaTo   -= useArea;

        // Update reaction rate for the adsorption reactions.
        updateAdsorptionReactionRateAtVertex(*ps, viFrom);
        updateAdsorptionReactionRateAtVertex(*ps, viTo);
    }
};

// Precondition: Species must have been set in the vertices
inline void setDiffusionForHalfEdge(
    SubSystem&                                     sys,
    MembraneMeshAttribute::MeshType&               mesh,
    MembraneMeshAttribute::MeshType::HalfEdgeIndex hei,
    const MembraneMeshChemistryInfo&               info
) {
    const auto viTo   = mesh.target(hei);
    const auto viFrom = mesh.target(mesh.opposite(hei));
    auto& objvTo   = mesh.attribute(viTo  ).vertex(sys);
    auto& objvFrom = mesh.attribute(viFrom).vertex(sys);
    auto& ch       = mesh.attribute(hei).halfEdge->cHalfEdge;
    auto& cvTo     = objvTo  .cVertex;
    auto& cvFrom   = objvFrom.cVertex;

    ch.diffusionReactions.clear();
    for(const auto& di : info.diffusion) {
        auto pr = std::make_unique< ReactionDy >(
            indicesToSpecies(cvFrom, { di.speciesIndex }),
            indicesToSpecies(cvTo,   { di.speciesIndex }),
            ReactionType::membraneDiffusion,
            di.diffusionCoeff,
            // volume is not set here
            1.0,
            // rate = diffCoeff * area^-1 * shapeFactor
            // where shapeFactor = 0.5 * (cot α + cot β)
            -1
        );
        // Add callback to update accessible areas.
        pr->connect(MembraneDiffusionReactionCallback {
            di.area,
            &sys,
            objvFrom.sysIndex,
            objvTo.sysIndex,
        });

        ch.diffusionReactions.push_back({
            std::move(pr),
            di.speciesIndex,
        });
    }
}

inline void clearDiffusionForHalfEdge(
    SubSystem&                                     sys,
    MembraneMeshAttribute::MeshType&               mesh,
    MembraneMeshAttribute::MeshType::HalfEdgeIndex hei
) {
    auto& ch = mesh.attribute(hei).halfEdge->cHalfEdge;

    ch.diffusionReactions.clear();
}


// Preconditions:
// - Up-to-date curvature mismatch energy computed by bending force field.
//
// Note:
// - If curvatureMismatchProteinIndex is -1, the curvature mismatch energy will be 0.
inline floatingpoint findVertexSpeciesFreeEnergy(
    SubSystem& sys,
    int        curvatureMismatchProteinIndex,
    int        vertexLoopIndex // indexed by membrane -> vertex sequence.
) {
    floatingpoint en = 0;
    if(curvatureMismatchProteinIndex != -1) {
        en += sys.allVertexCurvatureMismatchParams(curvatureMismatchProteinIndex, vertexLoopIndex).energy;
    }
    return en;
}


// Preconditions:
// - Species must have been set in the vertices
// - Reactions must have been set in the half edge
// - 1-ring areas around the vertices are updated
// - cot θ are set in neighbor triangles.
// - Vectorization (vertex loopIndex) is up-to-date.
//
// The diffusion is from the source to the target of the halfedge.
inline void setDiffusionRatesForHalfEdge(
    SubSystem&                                     sys,
    MembraneMeshAttribute::MeshType&               mesh,
    MembraneMeshAttribute::MeshType::HalfEdgeIndex hei
) {
    using MT = MembraneMeshAttribute::MeshType;

    const auto hei_o   = mesh.opposite(hei);
    const auto hei_n   = mesh.next(hei);
    const auto hei_on  = mesh.next(hei_o);

    const auto viFrom  = mesh.target(hei_o);
    const auto viTo    = mesh.target(hei);
    const auto viNext  = mesh.target(hei_n);
    const auto viPrev  = mesh.target(hei_on);

    const bool heiSide = mesh.isInTriangle(hei);
    const bool heiOSide = mesh.isInTriangle(hei_o);
    // The cot(theta) values.
    //   [ opposing hei, opposing hei_o, at v_from (hei side), at v_from (hei_o side) ]
    double ctHei = 0, ctHeiO = 0, ctVFromHei = 0, ctVFromHeiO = 0;
    if(heiSide) {
        ctHei = mesh.attribute(hei_n).gHalfEdge.cotTheta;
        ctVFromHei = mesh.attribute(mesh.prev(hei)).gHalfEdge.cotTheta;
    }
    if(heiOSide) {
        ctHeiO = mesh.attribute(hei_on).gHalfEdge.cotTheta;
        ctVFromHeiO = mesh.attribute(hei_o).gHalfEdge.cotTheta;
    }
    const double sumCotTheta = ctHei + ctHeiO;

    for(auto& r : mesh.attribute(hei).halfEdge->cHalfEdge.diffusionReactions) {

        // The energy component measures the gradient of free energy (in kT) at this point.
        double energyComponent = 0.0;
        {
            const auto cmIndex = sys.indexMapMembraneDiffusingSpeciesToCurvatureMismatchSetup[r.membraneDiffusingSpeciesIndex];
            const double energyFrom = findVertexSpeciesFreeEnergy(sys, cmIndex, mesh.attribute(viFrom).vertex(sys).loopIndex);
            const double energyTo   = findVertexSpeciesFreeEnergy(sys, cmIndex, mesh.attribute(viTo  ).vertex(sys).loopIndex);

            if(heiSide) {
                const double energyNext = findVertexSpeciesFreeEnergy(sys, cmIndex, mesh.attribute(viNext).vertex(sys).loopIndex);
                energyComponent += (energyTo - energyFrom) * ctHei + (energyTo - energyNext) * ctVFromHei;
            }
            if(heiOSide) {
                const double energyPrev = findVertexSpeciesFreeEnergy(sys, cmIndex, mesh.attribute(viPrev).vertex(sys).loopIndex);
                energyComponent += (energyTo - energyFrom) * ctHeiO + (energyTo - energyPrev) * ctVFromHeiO;
            }

            energyComponent /= 6;
        }

        // Using "volume frac" as volume
        r.reaction->setVolumeFrac(mesh.attribute(viFrom).gVertex.astar / 3);
        r.reaction->setRateMulFactor(
            std::max(0.5 * sumCotTheta - energyComponent / kT, 0.0),
            ReactionBase::RateMulFactorType::diffusionShape
        );
    }
}

inline void setSpeciesAndReactions(
    SubSystem&                        sys,
    MembraneMeshAttribute::MeshType&  mesh,
    const MembraneMeshChemistryInfo&  info
) {
    using namespace std;
    using MT = MembraneMeshAttribute::MeshType;

    // Set species
    //---------------------------------
    for(MT::VertexIndex vi {0}; vi < mesh.numVertices(); ++vi) {
        setSpeciesForVertex(sys, mesh, vi, info);
    }

    // Set reactions
    //---------------------------------

    // General reactions
    for(MT::VertexIndex vi {0}; vi < mesh.numVertices(); ++vi) {
        setInternalReactionsForVertex(sys, mesh, vi, info);
    }

    // Diffusion reactions
    for(MT::HalfEdgeIndex hei {0}; hei < mesh.numHalfEdges(); ++hei) {
        setDiffusionForHalfEdge(sys, mesh, hei, info);
    }
}

inline void clearSpeciesAndReactions(
    SubSystem&                       sys,
    MembraneMeshAttribute::MeshType& mesh
) {
    using namespace std;
    using MT = MembraneMeshAttribute::MeshType;

    // Remove reactions.
    //---------------------------------

    // General reactions.
    for(MT::VertexIndex vi {0}; vi < mesh.numVertices(); ++vi) {
        clearInternalReactionsForVertex(sys, mesh, vi);
    }

    // Diffusion reactions.
    for(MT::HalfEdgeIndex hei {0}; hei < mesh.numHalfEdges(); ++hei) {
        clearDiffusionForHalfEdge(sys, mesh, hei);
    }

    // Remove species.
    //---------------------------------
    for(MT::VertexIndex vi {0}; vi < mesh.numVertices(); ++vi) {
        clearSpeciesForVertex(sys, mesh, vi);
    }
}

// Preconditions:
//   - Species must have been set in all vertices
//   - Reactions must have been set in all half edges
//   - 1-ring areas around all vertices are updated
//   - cot θ are set in all triangles
inline void setReactionRates(
    SubSystem&                        sys,
    MembraneMeshAttribute::MeshType&  mesh
) {
    using namespace std;
    using MT = MembraneMeshAttribute::MeshType;

    // General reactions
    for(MT::VertexIndex vi {0}; vi < mesh.numVertices(); ++vi) {
        setInternalReactionRatesForVertex(sys, mesh, vi);
    }

    // Diffusion reactions
    for(MT::HalfEdgeIndex hei {0}; hei < mesh.numHalfEdges(); ++hei) {
        setDiffusionRatesForHalfEdge(sys, mesh, hei);
    }
}

// Auxiliary function to do unified operations on all reactions in mesh
//
// F has the following signature: void(ReactionDy &)
template< typename F >
inline void forEachReactionInMesh(
    SubSystem&                       sys,
    MembraneMeshAttribute::MeshType& mesh,
    F &&                             f
) {
    using MT = MembraneMeshAttribute::MeshType;

    // General / adsorption / desorption reactions
    for(MT::VertexIndex vi {0}; vi < mesh.numVertices(); ++vi) {
        auto& cv = mesh.attribute(vi).vertex(sys).cVertex;
        for(auto& pr : cv.reactions) {
            f(*pr);
        }
        for(auto& pr : cv.adsorptionReactions) {
            f(*pr);
        }
        for(auto& pr : cv.desorptionReactions) {
            f(*pr);
        }
    }

    // Diffusion reactions
    for(MT::HalfEdgeIndex hei {0}; hei < mesh.numHalfEdges(); ++hei) {
        for(auto& pr : mesh.attribute(hei).halfEdge->cHalfEdge.diffusionReactions) {
            f(*pr.reaction);
        }
    }
}


// Set chemistry parameters for membrane mesh.
inline void setChemistry(SubSystem& sys, MembraneMeshAttribute::MeshType& mesh, MembraneMeshChemistryInfo memChemInfo) {
    mesh.metaAttribute().chemInfo = std::move(memChemInfo);
    medyan::setSpeciesAndReactions(sys, mesh, mesh.metaAttribute().chemInfo);

    mesh.metaAttribute().isChemParamsSet = true;
}

// Clear chemistry parameters for membrane mesh.
inline void clearChemistry(SubSystem& sys, MembraneMeshAttribute::MeshType& mesh) {
    mesh.metaAttribute().isChemParamsSet = false;

    medyan::clearSpeciesAndReactions(sys, mesh);
}


// Initialize species copy numbers for membrane mesh.
//
// Preconditions:
// - The geometry of the mesh is up-to-date.
inline void initMeshSpecies(SubSystem& sys, MembraneMeshAttribute::MeshType& mesh, const MembraneInit& memInit) {
    using namespace std;
    using VI = MembraneMeshAttribute::MeshType::VertexIndex;

    // The mesh must have at least 1 vertex.
    if(mesh.numVertices() == 0) {
        log::error("The mesh does not have any vertices, so membrane diffusing species cannot be initialized.");
        throw runtime_error("Mesh does not contain vertices.");
    }

    // Initialize distribution object.
    std::vector<double> allVertexAreas(mesh.numVertices());
    for(VI vi {0}; vi < mesh.numVertices(); ++vi) {
        allVertexAreas[vi.index] = mesh.attribute(vi).gVertex.astar / 3;
    }
    std::discrete_distribution<Index> dist(allVertexAreas.begin(), allVertexAreas.end());

    // Initialize species copy numbers.
    for(auto& eachSpeciesInit : memInit.speciesInitVec) {
        int mol = SpeciesNamesDB::stringToInt(eachSpeciesInit.name);
        int speciesIndexInVertex = mesh.attribute(VI{0}).vertex(sys).cVertex.species.findSpeciesIndexByMolecule(mol);
        for(Index i = 0; i < eachSpeciesInit.copyNumber; ++i) {
            // Choose a random vertex based on the distribution, and add its associated species copy number by 1.
            VI vi { dist(Rand::eng) };
            auto& rspecies = mesh.attribute(vi).vertex(sys).cVertex.species.findSpeciesByIndex(speciesIndexInVertex)->getRSpecies();
            rspecies.setN(rspecies.getTrueN() + 1);
        }
    }
}

} // namespace medyan

#endif

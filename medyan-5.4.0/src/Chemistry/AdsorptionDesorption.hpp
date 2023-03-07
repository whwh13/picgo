#ifndef MEDYAN_Chemistry_AdsorptionDesorption_hpp
#define MEDYAN_Chemistry_AdsorptionDesorption_hpp

#include "Chemistry/ReactionDy.hpp"
#include "Chemistry/Species.h"
#include "Structure/CompartmentGrid.h"
#include "Structure/SubSystem.h"
#include "Structure/SurfaceMesh/FuncMembraneChem.hpp"
#include "Structure/SurfaceMesh/Vertex.hpp"
#include "SysParams.h"

namespace medyan {

// The callback when adsorption/desorption reaction happens.
// This callback updates the accessible area around the vertex.
struct AdsorptionDesorptionReactionCallback {
    // Area usage of the bound protein.
    FP useArea = 0.0;
    // Is this reaction adsorption?
    bool isAdsorption = false;
    // Pointer to system.
    SubSystem* ps = nullptr;
    // Stable index of the vertex.
    StableVectorIndex<Vertex> vIndex {};

    // The callback function.
    void operator()(ReactionBase* r) const {
        auto& accessibleArea = ps->vertices[vIndex].cVertex.accessibleArea;
        if(isAdsorption) {
            accessibleArea -= useArea;
        } else {
            accessibleArea += useArea;
        }

        // Update reaction rate for adsorption.
        updateAdsorptionReactionRateAtVertex(*ps, vIndex);
    }
};

// Initializes adsorption and desorption reactions for a mesh structure.
//
// Preconditions:
//   - Species must have been set in every vertex and in every compartment.
//   - Adsorption/desorption reactions must have not been set yet.
//   - Vertex positions are updated (such that they are registered correctly in compartments.
//
// Notes:
//   - All reactions are stored in vertices.
//   - Reaction rate multipliers are not set here. Reaction callbacks are not attached here.
inline void setAdsorptionDesorptionReactions(
    SubSystem&                                             sys,
    const ChemistryData::ReactionAdsorptionDesorptionInfo& info
) {
    using namespace std;

    auto& compartmentGrid = *sys.getCompartmentGrid();
    const int speciesMol3D = SpeciesNamesDB::stringToInt(info.speciesName3D);
    const int speciesMol2D = SpeciesNamesDB::stringToInt(info.speciesName2D);

    if(sys.vertices.size() == 0) return;

    // Cache index of 2D species in each vertex, given they are the same for all vertices.
    const int speciesIndex2D = sys.vertices.begin()->cVertex.species.findSpeciesIndexByMolecule(speciesMol2D);

    for(auto& pc : compartmentGrid.compartmentList) {
        Species* ps3 = pc->getSpeciesContainer().findSpeciesByMolecule(speciesMol3D);
        if(ps3 == nullptr) {
            ps3 = compartmentGrid.getSpeciesBulk().findSpeciesByMolecule(speciesMol3D);
        }
        if(ps3 == nullptr) {
            log::error("Cannot find 3D species {} in compartment or bulk.", info.speciesName3D);
            throw std::runtime_error("Invalid diffusing species.");
        }

        // Create reaction for all vertices.
        for(auto vi : pc->getVertices()) {
            auto& v = sys.vertices[vi];
            // Adsorption reaction.
            v.cVertex.adsorptionReactions.push_back(
                make_unique< ReactionDy >(
                    vector<Species*> { ps3 },
                    vector<Species*> { v.cVertex.species.findSpeciesByIndex(speciesIndex2D) },
                    ReactionType::membraneAdsorption,
                    // rate = onRate * volume^-1 * area
                    info.onRate,
                    // 3D volume is not set here.
                    1.0,
                    // Inversely depends on 3D volume.
                    -1
                )
            );

            // Desorption reaction.
            v.cVertex.desorptionReactions.push_back(
                make_unique< ReactionDy >(
                    vector<Species*> { v.cVertex.species.findSpeciesByIndex(speciesIndex2D) },
                    vector<Species*> { ps3 },
                    ReactionType::membraneDesorption,
                    // rate = offRate
                    info.offRate
                )
            );
        }
    }
}
// Same as before, but sets multiple reactions.
inline void setAdsorptionDesorptionReactions(
    SubSystem&                                             sys,
    const std::vector<ChemistryData::ReactionAdsorptionDesorptionInfo>& info
) {
    for(auto& eachInfo : info) {
        setAdsorptionDesorptionReactions(sys, eachInfo);
    }
}


// Clear all adsorption and desorption reactions.
inline void clearAdsorptionDesorptionReactions(SubSystem& sys) {
    for(auto& v : sys.vertices) {
        v.cVertex.adsorptionReactions.clear();
        v.cVertex.desorptionReactions.clear();
    }
}


// Auxiliary function to copy adsorption and desorption reactions to a new vertex.
//
// Notes:
// - Reaction callback is not copied.
inline void imitateAdsorptionDesorptionReactions(
    CVertex& cvFrom,
    CVertex& cvTo
) {
    using namespace std;

    cvTo.adsorptionReactions.clear();
    for(auto& pr : cvFrom.adsorptionReactions) {
        // Find the corresponding 2D species in the new vertex.
        const int speciesMol2D = pr->getRepRSpecies()[1].prs->getSpecies().getMolecule();
        auto pSpecies2D = cvTo.species.findSpeciesByMolecule(speciesMol2D);
        if(pSpecies2D == nullptr) {
            log::error("Species of molecule {} cannot be found in CVertex.", speciesMol2D);
            throw runtime_error("Invalid species molecule.");
        }
        cvTo.adsorptionReactions.push_back(
            make_unique< ReactionDy >(
                vector<Species*> { &(pr->getRepRSpecies()[0].prs->getSpecies()) },
                vector<Species*> { pSpecies2D },
                ReactionType::membraneAdsorption,
                pr->getBareRate(),
                pr->getVolumeFrac(),
                pr->getRateVolumeDepExp()
            )
        );
    }

    cvTo.desorptionReactions.clear();
    for(auto& pr : cvFrom.desorptionReactions) {
        // Find the corresponding 2D species in the new vertex.
        const int speciesMol2D = pr->getRepRSpecies()[0].prs->getSpecies().getMolecule();
        auto pSpecies2D = cvTo.species.findSpeciesByMolecule(speciesMol2D);
        if(pSpecies2D == nullptr) {
            log::error("Species of molecule {} cannot be found in CVertex.", speciesMol2D);
            throw runtime_error("Invalid species molecule.");
        }
        cvTo.desorptionReactions.push_back(
            make_unique< ReactionDy >(
                vector<Species*> { pSpecies2D },
                vector<Species*> { &(pr->getRepRSpecies()[1].prs->getSpecies()) },
                ReactionType::membraneDesorption,
                pr->getBareRate()
            )
        );
    }
}


// Auxiliary function to target adsorption/desorption 3D species in another compartment.
//
// Note:
// - All affected reactions in this function must be passivated, because
//     - This function does not handle propensity update.
//     - This function does not touch dependents in upstream reactions. They are only updated when these reactions are activated.
//     - This function does not create callbacks for new reactions.
// - This function handles species registration update, by creating a new reaction instead of modifying the current one.
inline void retarget3DSpeciesToCompartment(CVertex& cv, Compartment& comp) {
    using namespace std;

    for(auto& pr : cv.adsorptionReactions) {
        // Find 3D species, and only move if the species is diffusing.
        const auto& species3D = pr->getRepRSpecies()[0].prs->getSpecies();
        if(species3D.getType() == SpeciesType::DIFFUSING) {
            // Get molecule ID of 3D species.
            const int speciesMol3D = species3D.getMolecule();
            // Find the molecule in the new compartment.
            auto pSpecies3D = comp.findSpeciesByMolecule(speciesMol3D);
            // Move out the current reaction.
            auto prTemp = move(pr);
            // Create a new reaction using the original parameters.
            pr = make_unique< ReactionDy >(
                vector<Species*> { pSpecies3D },
                vector<Species*> { &(prTemp->getRepRSpecies()[1].prs->getSpecies()) },
                ReactionType::membraneAdsorption,
                prTemp->getBareRate(),
                // 3D volume is not set here.
                1.0,
                // Inversely depends on 3D volume.
                prTemp->getRateVolumeDepExp()
            );
            // The original reaction is destroyed at the end of the scope.
        }
    }

    for(auto& pr : cv.desorptionReactions) {
        // Find 3D species, and only move if the species is diffusing.
        const auto& species3D = pr->getRepRSpecies()[1].prs->getSpecies();
        if(species3D.getType() == SpeciesType::DIFFUSING) {
            // Get molecule ID of 3D species.
            const int speciesMol3D = species3D.getMolecule();
            // Find the molecule in the new compartment.
            auto pSpecies3D = comp.findSpeciesByMolecule(speciesMol3D);
            // Move out the current reaction.
            auto prTemp = move(pr);
            // Create a new reaction using the original parameters.
            pr = make_unique< ReactionDy >(
                vector<Species*> { &(prTemp->getRepRSpecies()[0].prs->getSpecies()) },
                vector<Species*> { pSpecies3D },
                ReactionType::membraneDesorption,
                prTemp->getBareRate()
            );
            // The original reaction is destroyed at the end of the scope.
        }
    }
}


// Set reaction rates for adsorption/desorption reactions.
inline void setAdsorptionDesorptionReactionRates(
    SubSystem&                sys,
    const ChemistryData&      chemData
) {
    using namespace std;
    using VI = Membrane::MeshType::VertexIndex;

    auto& info = chemData.reactionsAdsorptionDesorption;
    auto& sinfo = chemData.speciesMembraneDiffusing;

    auto& grid = *sys.getCompartmentGrid();

    if(sys.vertices.size() == 0) return;

    int mvIndex = 0;
    for(auto& membrane : sys.membranes) {
        auto& mesh = membrane.getMesh();
        for(VI vi{0}; vi < mesh.numVertices(); ++vi) {
            auto& v = mesh.attribute(vi).vertex(sys);
            auto& comp = *grid.vertexCellList.getHead(v.cellElement);
            const auto varea = mesh.attribute(vi).gVertex.astar / 3;
            auto& accessibleArea = v.cVertex.accessibleArea;

            // Find area occupied by proteins.
            floatingpoint proteinArea = 0;
            for(int si = 0; si < sinfo.size(); ++si) {
                proteinArea += v.cVertex.species.findSpeciesByIndex(si)->getRSpecies().getTrueN() * sinfo[si].area;
            }
            accessibleArea = varea - proteinArea;

            for(int ri = 0; ri < info.size(); ++ri) {
                auto& par = v.cVertex.adsorptionReactions[ri];
                auto& pdr = v.cVertex.desorptionReactions[ri];
                auto freeEnergy = findVertexSpeciesFreeEnergy(
                    sys,
                    sys.indexMapDesorptionReactionToCurvatureMismatchSetup[ri],
                    mvIndex
                );
                auto& species3D = par->getRepRSpecies()[0].prs->getSpecies();
                auto speciesIndex2D = v.cVertex.species.findSpeciesIndexByMolecule(par->getRepRSpecies()[1].prs->getSpecies().getMolecule());
                auto volume3D = species3D.getType() == SpeciesType::DIFFUSING
                    // Diffusing.
                    ? (grid.compartmentVolume * comp.getVolumeFrac())
                    // Bulk.
                    : (grid.compartmentVolume * grid.getCompartments().size());

                par->setRateMulFactor(std::max<FP>(0, accessibleArea), ReactionBase::RateMulFactorType::mechanochemical);
                par->setVolumeFrac(volume3D);
                par->clearSignaling();
                par->connect(AdsorptionDesorptionReactionCallback {
                    sinfo[speciesIndex2D].area,
                    true,
                    &sys,
                    v.sysIndex,
                });

                pdr->setRateMulFactor(
                    exp(freeEnergy / kT),
                    ReactionBase::RateMulFactorType::mechanochemical
                );
                pdr->clearSignaling();
                pdr->connect(AdsorptionDesorptionReactionCallback {
                    sinfo[speciesIndex2D].area,
                    false,
                    &sys,
                    v.sysIndex,
                });
            }

            ++mvIndex;
        }
    }
}

} // namespace medyan

#endif

#ifndef MEDYAN_Chemistry_EmissionAbsorption_hpp
#define MEDYAN_Chemistry_EmissionAbsorption_hpp

#include "Chemistry/ReactionDy.hpp"
#include "Chemistry/Species.h"
#include "SysParams.h"

namespace medyan {

// Note:
// - If the species is used in any other reaction not contained in this struct, the those reactions should be deleted before the struct is deleted.
// - The 2 reactions must correspond to: this_species <-> diffusing_or_bulk_species.
struct EmissionAbsorptionContainer {
    using DyRateType = ReactionEmissionAbsorptionSetup::DyRateType;
    std::unique_ptr< Species > pspecies;
    std::unique_ptr< ReactionDy > prEmi;
    std::unique_ptr< ReactionDy > prAbs;
    DyRateType dyRateType = DyRateType::none;
};

template< typename Context >
auto setEmiAbs(Context& sys, Index compartmentIndex, const std::vector<ReactionEmissionAbsorptionSetupInit>& vecEmiAbsInit, const ChemistryData& chemData) {
    using namespace std;

    std::vector<EmissionAbsorptionContainer> res;

    auto& grid = *sys.getCompartmentGrid();
    auto& comp = grid.getCompartment(compartmentIndex);
    auto& chemSim = *sys.pChemSim;

    res.reserve(vecEmiAbsInit.size());
    for(auto& ri : vecEmiAbsInit) {
        auto& r = ri.setup;
        auto speciesGeneralIndex = chemData.findIndexSpeciesGeneral(r.speciesName1).value();
        auto rType =
            chemData.speciesGeneral[speciesGeneralIndex].rspeciesType == "CONST"
            ? RSpeciesType::CONST : RSpeciesType::REG;
        const bool useDiffusing = r.speciesType2 == "diffusing";
        auto ps3d =   useDiffusing ? comp.findSpeciesByName(r.speciesName2) : grid.findSpeciesBulkByName(r.speciesName2);
        auto volume = useDiffusing ? grid.compartmentVolume : grid.gridVolume;

        EmissionAbsorptionContainer ea;
        // Setup species.
        ea.pspecies = make_unique<Species>(r.speciesName1, ri.initNumSpecies1, max_ulim, SpeciesType::general, rType);
        // Setup reactions.
        ea.prEmi = make_unique<ReactionDy> (
            vector<Species*> { ea.pspecies.get() },
            vector<Species*> { ps3d },
            ReactionType::emission,
            r.emiRateConst
        );
        ea.prAbs = make_unique<ReactionDy> (
            vector<Species*> { ps3d },
            vector<Species*> { ea.pspecies.get() },
            ReactionType::absorption,
            r.absRateConst,
            volume * comp.getVolumeFrac(),
            -1
        );
        // Setup dynamic rates.
        ea.dyRateType = r.dyRateType;

        // Activate the reactions.
        chemSim.addReaction(ea.prEmi.get());
        ea.prEmi->activateReaction();
        chemSim.addReaction(ea.prAbs.get());
        ea.prAbs->activateReaction();

        res.push_back(move(ea));
    }

    return res;
}

// Update the diffusing species in emi-abs reactions.
template< typename Context >
void updateCompartmentEmiAbs(Context& sys, EmissionAbsorptionContainer& ea, Index newci) {
    auto& grid = *sys.getCompartmentGrid();
    auto& chemSim = *sys.pChemSim;

    auto ps3d = &ea.prEmi->getRepRSpecies()[1].prs->getSpecies();
    if(ps3d->getType() == SpeciesType::DIFFUSING) {
        auto mol3d = ps3d->getMolecule();
        // Find in the new compartment the same molecule.
        auto& newc = grid.getCompartment(newci);
        auto newps3d = newc.findSpeciesByMolecule(mol3d);
        // Update the reactions.
        {
            ea.prEmi->passivateReaction();
            chemSim.removeReaction(ea.prEmi.get());

            auto oldEmi = std::move(ea.prEmi);
            ea.prEmi = make_unique<ReactionDy> (
                vector<Species*> { &oldEmi->getRepRSpecies()[0].prs->getSpecies() },
                vector<Species*> { newps3d },
                oldEmi->getReactionType(),
                oldEmi->getBareRate(),
                oldEmi->getVolumeFrac(),
                oldEmi->getRateVolumeDepExp()
            );

            chemSim.addReaction(ea.prEmi.get());
            ea.prEmi->activateReaction();
        }
        {
            ea.prAbs->passivateReaction();
            chemSim.removeReaction(ea.prAbs.get());

            auto oldAbs = std::move(ea.prAbs);
            ea.prAbs = make_unique<ReactionDy> (
                vector<Species*> { newps3d },
                vector<Species*> { &oldAbs->getRepRSpecies()[1].prs->getSpecies() },
                oldAbs->getReactionType(),
                oldAbs->getBareRate(),
                grid.compartmentVolume * newc.getVolumeFrac(),
                oldAbs->getRateVolumeDepExp()
            );

            chemSim.addReaction(ea.prAbs.get());
            ea.prAbs->activateReaction();
        }
    }
}

} // namespace medyan

#endif

#ifndef MEDYAN_Structure_CompartmentChem_hpp
#define MEDYAN_Structure_CompartmentChem_hpp

#include "Chemistry/ChemSim.h"
#include "Structure/CompartmentGrid.h"

namespace medyan {

//-----------------------------------------------------------------------------
// Compartment reaction management.
//-----------------------------------------------------------------------------

// Auxiliary function to compute base diffusion reaction rate from diffusion coefficient.
inline floatingpoint getDiffusionReactionRate(floatingpoint diffusionCoefficient, const CompartmentGrid& grid, int dirNeighbor) {
    auto& length = grid.compartmentLengths.at(dirNeighbor / 2);
    return diffusionCoefficient / (length * length);
}

// Generate diffusion reactions between this compartment and another.
//
// Parameters:
// - rxns:        The reaction base vector to append to.
// - grid:        The compartment grid object.
// - cindex1:     The index of the 1st compartment. Must be active.
// - dirNeighbor: Neighbor direction from 1st compartment to the 2nd compartment. Must point to a valid neighbor.
// - forwardOnly: If true, only c1->c2 diffusion will be added. Otherwise, c1<->c2 will both be added.
inline void generateDiffusionReactions(std::vector<ReactionBase*>& rxns, CompartmentGrid& grid, int cindex1, int dirNeighbor, bool forwardOnly) {
    const auto dirOppo = Compartment::getOppositeNeighborDirection(dirNeighbor);
    auto& c1 = grid.getCompartment(cindex1);
    const auto cindex2 = c1.getNeighborIndices()[dirNeighbor];
    auto& c2 = grid.getCompartment(cindex2);

    for(auto &sp1 : c1.getSpeciesContainer().species()) {
        int mol1 = sp1->getMolecule();
        const auto diffRate = getDiffusionReactionRate(c1.getDiffusionCoefficient(mol1), grid, dirNeighbor);
        if(diffRate<0)  continue;

        if(c2.isActivated()) {
            // Scale the diffusion rate according to the contacting areas.
            double scaleFactor = 0.5 * (c1.getPartialArea()[dirNeighbor] + c2.getPartialArea()[dirOppo]) / grid.compartmentAreas[dirNeighbor / 2];

            const float volumeFrac = c1.getVolumeFrac();
            // cout << "VolumeFraction = " << volumeFrac << endl;
            Species *sp2 = c2.getSpeciesContainer().findSpeciesByMolecule(mol1);

            auto r = std::make_unique<DiffusionReaction>(std::initializer_list<Species*>{sp1.get(), sp2}, diffRate, false, volumeFrac);
            r->setRateMulFactor(scaleFactor, ReactionBase::diffusionShape);
            rxns.push_back(r.get());
            c1.addDiffusionReaction(std::move(r));

            if(!forwardOnly) {
                // Generate inward diffusion reaction
                auto r = std::make_unique<DiffusionReaction>(std::initializer_list<Species*>{sp2, sp1.get()}, diffRate, false, c2.getVolumeFrac());
                r->setRateMulFactor(scaleFactor, ReactionBase::diffusionShape);
                rxns.push_back(r.get());
                c2.addDiffusionReaction(std::move(r));
            }
        }
    }
}

// Generate all diffusion reactions for this compartment and its neighbors.
// Returns a vector of reactionbases that was just added.
inline std::vector<ReactionBase*> generateAllDiffusionReactions(CompartmentGrid& grid, int cindex, bool outwardOnly) {
    auto& c = grid.getCompartment(cindex);
    
    std::vector<ReactionBase*> rxns;

    if(c.isActivated()) {
        for(int dir = 0; dir < Compartment::numNeighborDirections; ++dir) {
            auto cnindex = c.getNeighborIndices()[dir];
            if(cnindex != -1) {
                generateDiffusionReactions(rxns, grid, cindex, dir, outwardOnly);
            }
        }
    }
    return rxns;
}

// Remove diffusion reactions from the neighboring compartment to the specified compartment.
//
// Parameters:
// - grid:        The compartment grid object.
// - cindex:      The index of the 1st compartment.
// - dirNeighbor: Neighbor direction from 1st compartment to the 2nd compartment. Must point to a valid neighbor.
// - chem:        The chemical simulation engine.
inline void removeDiffusionReactions(CompartmentGrid& grid, int cindex, int dirNeighbor, ChemSim& chem) {
    auto& c = grid.getCompartment(cindex);
    auto cnindex = c.getNeighborIndices()[dirNeighbor];
    auto& cn = grid.getCompartment(cnindex);

    //look for neighbor's diffusion reactions
    vector<ReactionBase*> to_remove;

    for(auto &r : cn.getDiffusionReactionContainer().reactions()) {

        auto rs = r.get()->rspecies()[1];
        if(rs->getSpecies().getParent() == &c) {

            r->passivateReaction();

            chem.removeReaction(r.get());

            to_remove.push_back(r.get());
        }

    }

    //remove them
    for(auto &r : to_remove)
        cn.getDiffusionReactionContainer().removeReaction(r);

}

// Remove all diffusion reactions between this compartment and its neighbors.
inline void removeAllDiffusionReactions(CompartmentGrid& grid, int cindex, ChemSim& chem) {
    auto& c = grid.getCompartment(cindex);

    //remove all diffusion reactions that this has ownership of
    for(auto &r : c.getDiffusionReactionContainer().reactions()) {
        r->passivateReaction();
        chem.removeReaction(r.get());
    }

    c.getDiffusionReactionContainer().clear();

    // Remove neighboring diffusing reactions with this compartment.
    for(int dir = 0; dir < Compartment::numNeighborDirections; ++dir) {
        if(c.getNeighborIndices()[dir] != -1) {
            removeDiffusionReactions(grid, cindex, dir, chem);
        }
    }
}

// Transfer all species copy numbers from this compartment to neighboring active compartments.
// If no neighboring active compartments are present, throw an error.
inline void transferSpecies(CompartmentGrid& grid, int cindex, int i) {
    //i axis
    //-1 all directions
    //0 X
    //1 Y
    //2 Z
    //3 all directions
    //get active neighbors
    auto& c = grid.getCompartment(cindex);
    vector<Compartment*> activeNeighbors;

    for(auto cnindex : c.getNeighborIndices()) if(cnindex != -1) {
        auto& cn = grid.getCompartment(cnindex);
        auto ncoord = cn.centerCoord;

        if(cn.isActivated()){
            if(i < 0 || i == 3)
                activeNeighbors.push_back(&cn);
            else if(areEqual(distance(ncoord, c.centerCoord), abs(c.centerCoord[i]-ncoord[i])))
                activeNeighbors.push_back(&cn);
        }
    }

    assert(activeNeighbors.size() != 0
           && "Cannot transfer species to another compartment... no neighbors are active");
    if(i >= 0 && i<3 && activeNeighbors.size()>1){
        cout<<"Error transferring species along an axis. More than 1 neighbor. Exiting. "<< endl;
        exit(EXIT_FAILURE);
    }

    //go through species
    Species* sp_neighbor;
    vector<Species*> sp_neighbors;

    for(auto &sp : c.getSpeciesContainer().species()) {

        int copyNumber = sp->getN();
        auto nit = activeNeighbors.begin();

        if(sp->getType() == SpeciesType::DIFFUSING) {
            while(copyNumber > 0) {
                sp->down();

                //choose a random active neighbor
                auto neighbor = *nit;

                sp_neighbor = neighbor->findSpeciesByName(sp->getName());

                //add to list if not already
                auto spit = find(sp_neighbors.begin(),
                                 sp_neighbors.end(),
                                 sp_neighbor);

                if(spit == sp_neighbors.end())
                    sp_neighbors.push_back(sp_neighbor);

                //increase copy number

                sp_neighbor->up();

                //reset if we've looped through
                if(++nit == activeNeighbors.end())
                    nit = activeNeighbors.begin();
                copyNumber--;

            }
        }

        //activate all reactions changed
        for(auto spn : sp_neighbors)
            spn->updateReactantPropensities();
        for(auto &sp : c.getSpeciesContainer().species())
            sp->updateReactantPropensities();
    }
}

// Share all species copy numbers from this compartment to neighboring active compartments.
// If no neighboring active compartments are present, throw an error.
inline void shareSpecies(CompartmentGrid& grid, int cindex, int i) {
    //i axis
    //-1 all directions
    //0 X
    //1 Y
    //2 Z
    //3 all directions
    //get active neighbors
    auto& c = grid.getCompartment(cindex);
    vector<Compartment*> activeNeighbors;

    for(auto cnindex : c.getNeighborIndices()) if(cnindex != -1) {
        auto& cn = grid.getCompartment(cnindex);
        auto ncoord = cn.centerCoord;
        if(cn.isActivated()){
            if(i < 0 || i == 3)
                activeNeighbors.push_back(&cn);
            else if(areEqual(distance(ncoord, c.centerCoord), abs(c.centerCoord[i]-ncoord[i])))
                activeNeighbors.push_back(&cn);
        }
    }

    assert(activeNeighbors.size() != 0
           && "Cannot share species to another compartment... no neighbors are active");
    if(i >= 0 && i<3 && activeNeighbors.size()>1){
        cout<<"Error sharing species along an axis. More than 1 neighbor. Exiting."<< endl;
        exit(EXIT_FAILURE);
    }
    //go through species
    Species* sp_neighbor;
    vector<Species*> sp_neighbors;

    for(auto &sp : c.getSpeciesContainer().species()) {
        auto nit = activeNeighbors.begin();
        auto neighbor = *nit;
        sp_neighbor = neighbor->findSpeciesByName(sp->getName());
        int copyNumber = sp_neighbor->getN();
        int lowerlimit = (int) sp_neighbor->getN()/2;
        if(sp->getType() == SpeciesType::DIFFUSING) {
            while(copyNumber > lowerlimit) {
                sp_neighbor->down();

                //add to list if not already
                auto spit = find(sp_neighbors.begin(),
                                 sp_neighbors.end(),
                                 sp_neighbor);

                if(spit == sp_neighbors.end())
                    sp_neighbors.push_back(sp_neighbor);

                //increase copy number
                sp->up();
                //reset if we've looped through
                if(++nit == activeNeighbors.end())
                    nit = activeNeighbors.begin();
                neighbor = *nit;
                sp_neighbor = neighbor->findSpeciesByName(sp->getName());
                copyNumber--;

            }
        }

        //activate all reactions changed
        for(auto spn : sp_neighbors)
            spn->updateReactantPropensities();
        for(auto &sp : c.getSpeciesContainer().species())
            sp->updateReactantPropensities();

    }
}



//-----------------------------------------------------------------------------
// Compartment status update.
//-----------------------------------------------------------------------------

// Activate a compartment.
inline void activate(CompartmentGrid& grid, int cindex, ChemSim& chem, Compartment::ActivateReason reason = Compartment::ActivateReason::Whole) {
    /**************************************************************************
    The diffusion-reactions with the already activated neighbors would be added
    for both directions.
    **************************************************************************/

    auto& c = grid.getCompartment(cindex);
    assert(!c.isActivated() && "Compartment is already activated.");

    //set marker
    c.setAsActive();
    //add all diffusion reactions
    auto rxns = generateAllDiffusionReactions(grid, cindex, false);

    for(auto &r : rxns) {
        chem.addReaction(r);
        r->activateReaction(); // Conditionally activate the new diffusion reactions
    }

    if(reason == Compartment::ActivateReason::Whole) {
        shareSpecies(grid, cindex, SysParams::Boundaries().transfershareaxis);
    }

}

// Update activation of a compartment.
inline void updateActivation(CompartmentGrid& grid, int cindex, ChemSim& chem, Compartment::ActivateReason reason = Compartment::ActivateReason::Whole) {
    auto& c = grid.getCompartment(cindex);
    const auto volumeFrac = c.getVolumeFrac();

    if(c.isActivated()) {
        // Update the reaction rates for diffusions in both directions.
        for(int dir = 0; dir < Compartment::numNeighborDirections; ++dir) {
            auto cnindex = c.getNeighborIndices()[dir];
            if(cnindex == -1) {
                continue;
            }

            auto& cn = grid.getCompartment(cnindex);
            if(!cn.isActivated()) {
                continue;
            }

            // For any activated neighbor.
            for(auto &sp_this : c.getSpeciesContainer().species()) {
                int molecule = sp_this->getMolecule();
                const auto baseDiffRate = getDiffusionReactionRate(c.getDiffusionCoefficient(molecule), grid, dir);
                if(baseDiffRate<0)  continue;

                Species *sp_neighbor = cn.getSpeciesContainer().findSpeciesByMolecule(molecule);

                auto dirOppo = Compartment::getOppositeNeighborDirection(dir);
                // Scale the diffusion rate according to the contacting areas
                double scaleFactor =
                    0.5 * (c.getPartialArea()[dir] + cn.getPartialArea()[dirOppo]) /
                    grid.compartmentAreas[dir / 2];

                // Update outward reaction rate
                for(auto& r: c.getDiffusionReactionContainer().reactions())
                    if(sp_this.get() == &r->rspecies()[0]->getSpecies() && sp_neighbor == &r->rspecies()[1]->getSpecies()) {
                        r->setVolumeFrac(volumeFrac);
                        r->recalcRateVolumeFactor();
                        r->setRateMulFactor(scaleFactor, ReactionBase::diffusionShape);
                    }
                // We also update inward reaction rate here to ensure that neighbors are always on the same page.
                // Update inward reaction rate
                for(auto& r: cn.getDiffusionReactionContainer().reactions())
                    if(sp_this.get() == &r->rspecies()[1]->getSpecies() && sp_neighbor == &r->rspecies()[0]->getSpecies()) {
                        r->setVolumeFrac(cn.getVolumeFrac());
                        r->recalcRateVolumeFactor();
                        r->setRateMulFactor(scaleFactor, ReactionBase::diffusionShape);
                    }
            }

        }
    } else {
        activate(grid, cindex, chem, reason);
    }

    // Update the internal reaction rates
    for(auto& r: c.getInternalReactionContainer().reactions()) {
        r->setVolumeFrac(volumeFrac);
        r->recalcRateVolumeFactor();
    }
}


// Deactivate a compartment.
// Has the following sid effects:
// 0) Initially checks that all cylinders are removed
//    from this compartment. A compartment cannot be deactivated
//    unless this condition is already true.
// 1) Transfers copy numbers of all diffusing species from this
//    compartment and moves them to a neighboring active compartment.
//    If there are no neighboring active compartments, an error will result.
// 2) Removes all diffusion reactions involving diffusing species
//    in this compartment.
inline void deactivate(CompartmentGrid& grid, int cindex, ChemSim& chem, bool init) {
    auto& c = grid.getCompartment(cindex);

    // Assert no cylinders in this compartment.
    assert((c.getCylinders().size() == 0)
           && "Compartment cannot be deactivated when containing active cylinders.");

    //set marker
    c.setAsInactive();

    if(!init) transferSpecies(grid, cindex, SysParams::Boundaries().transfershareaxis);
    removeAllDiffusionReactions(grid, cindex, chem);
}

} // namespace medyan

#endif


//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_CController_h
#define MEDYAN_CController_h

#include <memory>

#include "common.h"
#include "Compartment.h"
#include "Chemistry/ChemSim.h"
#include "Structure/CompartmentChem.hpp"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Structure/SurfaceMesh/FuncMembraneChem.hpp"
#include "Structure/SubSystem.h"


namespace medyan {
// Forward declarations.
class ChemManager;

/// Used to intialize, control, and run the chemical components of a simulation

/*!
 *  ChemController is a class used by the SubSystem to instantiate, control, and run
 *  the chemical dynamics of a simulation. It has functions to initialize a chemical 
 *  system, which, based on a choice of the reaction-diffusion algorithm as well as the 
 *  type of manager which controls the reactions in the simulation, as well as run the 
 *  chemical components of the simulation.
 *
 *  The controller initializes all chemical objects used, including ChemSim
 *  and ChemManager to the correct implementations, given that they are implemented.
 */
class CController {
   
private:
    SubSystem* _subSystem; ///< A pointer to the SubSystem
    
    //@{
    /// Holding pointers to control all chemistry.
    ChemManager* _chemManager;
    //@}
    
public:

    /// Constructor which sets subsystem pointer
    CController(SubSystem* s) : _subSystem(s) {}
    
    ///Activate a compartment. Wrapper function for Compartment::activate().
    void activate(CompartmentGrid& grid, int cindex, Compartment::ActivateReason reason = Compartment::ActivateReason::Whole) { medyan::activate(grid, cindex, *_subSystem->pChemSim, reason); }
    /// Update activation of a compartment. Wrapper function for Compartment::updateActivation()
    void updateActivation(CompartmentGrid& grid, int cindex, Compartment::ActivateReason reason = Compartment::ActivateReason::Whole) { medyan::updateActivation(grid, cindex, *_subSystem->pChemSim, reason); }
    ///Deactivate a compartment. Wrapper function for Compartment::deactivate().
    void deactivate(CompartmentGrid& grid, int cindex, bool init=false) { medyan::deactivate(grid, cindex, *_subSystem->pChemSim, init); }
    
    /// Initialize the ChemSim algorithm as well as the ChemManager
    ///@param chemAlgorithm - a string defining the chemical algorithm to be used
    ///@param chemInitializer - a string defining the chemical manager used
    void initialize(string& chemAlgorithm, ChemistryData& chem, DissipationTracker* dt, medyan::SimulConfig& sc);

    void initializerestart(floatingpoint restartime, floatingpoint _minimizationTime, ChemistryData& chemData);

    // Things to be done before chemistry simulation
    void beforeRun() const {
        // Link all membrane mesh reactions and activate them
        for(auto& m : _subSystem->membranes) {
            medyan::forEachReactionInMesh(*_subSystem, m.getMesh(), [this](ReactionDy& r) {
                _subSystem->pChemSim->addReaction(&r);
                r.activateReaction();
            });
        }
    }
    // Things to be done after chemistry simulation
    void afterRun() const {
        // Unlink all membrane mesh reactions
        for(auto& m : _subSystem->membranes) {
            medyan::forEachReactionInMesh(*_subSystem, m.getMesh(), [this](ReactionDy& r) {
                r.passivateReaction();
                _subSystem->pChemSim->removeReaction(&r);
            });
        }
    }

    ///Run chemistry for a given amount of time
    bool run(FP time, ChemistryData& chemData, const SimulConfig& conf);
    
    ///Run chemistry for a given number of reaction steps
    bool runSteps(int steps, ChemistryData& chemData, const SimulConfig& conf);
    
    ///Remove set of reactions at runtime, specified by input
    void removeReactions();

    bool crosschecktau();
    
    ChemSim* getCS() const { return _subSystem->pChemSim.get(); }

    static ofstream _crosscheckdumpFilechem;
};

} // namespace medyan

#endif

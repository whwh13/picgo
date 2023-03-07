
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

#ifndef MEDYAN_Controller_h
#define MEDYAN_Controller_h

#include <memory> // unique_ptr

#include "common.h"
#include "Output.h"
#include "MedyanConfig.hpp"
#include "Restart.h"
#include "Chemistry/DissipationTracker.h"
#include "Controller/CController.h"
#include "Controller/DRController.h"
#include "Controller/GController.h"
#include "Controller/MController.h"
#include "Controller/SpecialProtocols.hpp"
#include "Structure/SubSystem.h"
#include "Structure/Trackable.h"
#include "Structure/SurfaceMesh/AdaptiveMesh.hpp"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Structure/SurfaceMesh/MembraneRegion.hpp"


namespace medyan {

/// Used to initialize, manage, and run an entire simulation.

/*!
 *  The Controller is initialized in the main program, and initializes the SubSystem
 *  given an initial input directory. After initialization of all member controllers,
 *  the Controller class can run a simulation given the already read input parameters, 
 *  by iteratively switching between mechanical equilibration and stochastic chemical 
 *  steps. The Controller is also responsible for updating positions, reaction rates,
 *  and neighbors lists through its corresponding sub-controllers and the Subsystem.
 */
class Controller {

private:
    
    SubSystem _subSystem; ///< A pointer to the subsystem that this controls

    MController _mController;   ///< Mechanical controller used
    CController _cController;   ///< Chemical controller used
    GController _gController;   ///< Geometry controller used
    DRController _drController; ///< Dynamic rate controller used
    SpecialProtocolManager specialProtocolManager_;
    
    vector<std::unique_ptr< Output >> _outputs; ///< Vector of specified outputs
    vector<std::unique_ptr< Output >> _outputdump; ///<Vector of outputs that correspond to a datadump
    // Output for snapshots.
    SnapshotOutput output_;
    RockingSnapshot* _rSnapShot;

    floatingpoint _runTime;          ///< Total desired runtime for simulation

    floatingpoint _minimizationTime;  ///< Frequency of mechanical minimization

    std::unique_ptr<adaptive_mesh::MembraneMeshAdapter> _meshAdapter; ///< Used in adaptive remeshing algorithm
    std::unique_ptr< MembraneRegion< Membrane > > _regionInMembrane; // The region that is inside the outermost membrane

    floatingpoint _initialSlowDownTime = 0.0; ///< Time cut off for more frequent minimization

    
    
    //@{
    /// Same parameter set as timestep, but in terms of chemical
    /// reaction steps. Useful for small runs and debugging.
    floatingpoint _runSteps;
    floatingpoint _snapshotSteps;
    
    floatingpoint _minimizationSteps;
    floatingpoint _neighborListSteps;
    vector<Compartment*> activatecompartments;
    multimap<int,Compartment*> fCompmap;
    multimap<int,Compartment*> bCompmap;

    Restart* _restart;
    //@}
    floatingpoint bounds[2], bounds_prev[2];
    ///INITIALIZATION HELPER FUNCTIONS
    
    /// Set up an initial configuration of a network
    /// For now, only [Bubbles](@ref Bubble) and [Filaments](@ref Filament)
    /// can be initialized before the simulation begins. Any other elements
    /// should be initialized in this function.
    void setupInitialNetwork(const CommandLineConfig&, SimulConfig&);
    
    /// Setup any special structures needed
    void setupSpecialStructures(medyan::SimulConfig&);
    
    ///RUNTIME HELPER FUNCTIONS
    
    /// Move the boundary based on the timestep
    void moveBoundary(floatingpoint deltaTau);

    // Update compartments activity based on boundary, membrane, etc.
    // Also update partial volumes and reaction rates
    // Used after each mechanical minimization
    void updateActiveCompartments();
    
    ///Activate/deactivate compartments based on the longest filament (along Xaxis).
    void activatedeactivateComp();
    void ControlfrontEndCompobsolete();
    void ControlbackEndCompobsolete();
    void ControlfrontbackEndComp();
    /// Update the positions of all elements in the system
    void updatePositions(double& tp);

    /// Update the reaction rates of all elements in the system
    void updateReactionRates(const SimulConfig&);
    
    /// Update neighbors lists, called in run
    void updateNeighborLists(const CommandLineConfig&, SimulConfig& sc);

    // Initialize special protocols.
    void initializeSpecialProtocols(const SimulConfig&);
    /// Execute any special protocols needed, for example,
    /// making Linker and Filament species static
    void executeSpecialProtocols(const SimulConfig& conf) {
        specialProtocolManager_.execute(_subSystem, conf);
    }

    /// Reset counters on all elements in the system
    void resetCounters();

    /// Helper function to remesh the membranes
    void membraneAdaptiveRemesh();
    
    double tp = SysParams::Chemistry().makeRateDependTime;
    double threforce = SysParams::Chemistry().makeRateDependForce;
    
public:
    void printtimers();
    floatingpoint chemistrytime = 0.0;
    floatingpoint minimizationtime = 0.0;
    floatingpoint nltime = 0.0;
    floatingpoint nl2time = 0.0;
    floatingpoint bmgrvectime = 0.0;
    floatingpoint bmgrtime = 0.0;
    floatingpoint rxnratetime = 0.0;
    floatingpoint updateposition = 0.0;
    floatingpoint outputtime =0.0;
    floatingpoint specialtime = 0.0;
    floatingpoint updatepositioncylinder = 0.0;
    floatingpoint updatepositionmovable=0.0;
    floatingpoint whileloop = 0.0;

    Controller() :
        _mController(&_subSystem),
        _cController(&_subSystem),
        _gController(&_subSystem)
    {
        Trackable::_subSystem = &_subSystem;
    }

    ~Controller() {
        finalize();
    }

    // Initialize the system, given parameters from the command line.
    // Returns the simulation configuration.
    SimulConfig initialize(const CommandLineConfig& cmdConfig);

    // Run the simulation.
    void run(const CommandLineConfig&, SimulConfig&);

    // Finalize the system.
    // Some elements of the system cannot be torn down automatically, and this function takes care of those.
    void finalize() {
        // Inter-element cleanup.
        log::debug("Clearing adsorption desorption reactions...");
        clearAdsorptionDesorptionReactions(_subSystem);

        // Element cleanup.
        log::debug("Tearing down fixed vertex attachments...");
        for(auto it = _subSystem.fixedVertexAttachments.begin(); it != _subSystem.fixedVertexAttachments.end(); ++it) {
            SubSystemFunc{}.removeTrackable<FixedVertexAttachment>(_subSystem, _subSystem.fixedVertexAttachments.indexat(it));
        }
        log::debug("Tearing down membranes...");
        for(auto it = _subSystem.membranes.begin(); it != _subSystem.membranes.end(); ++it) {
            clearChemistry(_subSystem, it->getMesh());
            SubSystemFunc{}.removeTrackable<Membrane>(_subSystem, _subSystem.membranes.indexat(it));
        }

        log::debug("Tearing down bubbles...");
        for(auto it = _subSystem.afms.begin(); it != _subSystem.afms.end(); ++it) {
            it->clearChemistry(_subSystem);
            SubSystemFunc{}.removeTrackable<AFM>(_subSystem, _subSystem.afms.indexat(it));
        }
        for(auto it = _subSystem.mtocs.begin(); it != _subSystem.mtocs.end(); ++it) {
            it->clearChemistry(_subSystem);
            SubSystemFunc{}.removeTrackable<MTOC>(_subSystem, _subSystem.mtocs.indexat(it));
        }
    }
};

} // namespace medyan

#endif

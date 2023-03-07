
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

#ifndef MEDYAN_ChemSim_h
#define MEDYAN_ChemSim_h

#include "common.h"
    
namespace medyan {
//FORWARD DECLARATIONS
class ReactionBase;
class DissipationTracker;


/// Used to manage running a network of chemical reactions.

/*! ChemSim implements a Strategy pattern, allowing custom algorithms for running 
 *  stochastic chemical simulations. It itself has a pointer to a  single static 
 *  implementation of ChemSimImpl. After the specific algorithm is chosen and ChemSim 
 *  is instantiated, ChemSim can be used to manage simulations, through such methods 
 *  as run(steps) etc.
 */
class ChemSim {
public:
    DissipationTracker* dt = nullptr;

    virtual ~ChemSim() = default;
    
    /// After all initial reactions have been added via addReaction(...) method,
    /// invoke initialize() prior to invoking run()
    virtual void initialize() = 0;

    //necessary to call during restart to reset global time to necessary time.
    virtual void initializerestart(floatingpoint time) = 0;

    /// Add Reaction *r to the chemical network which needs to be simulated
    virtual void addReaction(ReactionBase *r) = 0;
    
    /// Remove Reaction *r from the simulated chemical network 
    virtual void removeReaction(ReactionBase *r) = 0;
    
    /// Run the chemical dynamics for a set amount of time
    virtual bool run(floatingpoint time) = 0;
    
    /// Run the chemical dynamics for a set amount of reaction steps
    virtual bool runSteps(int steps) = 0;
    
    /// Mainly used for debugging: print chemical reactions in the network at
    /// this moment
    virtual void printReactions() const = 0;

    /// Cross checks all reactions in the network for firing time.
    virtual bool crosschecktau() const = 0;
    
    DissipationTracker* getDT() const { return dt; }
};

} // namespace medyan

#endif

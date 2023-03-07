
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

#ifndef MEDYAN_ChemGillespieImpl_h
#define MEDYAN_ChemGillespieImpl_h

#include <vector>
#include <random>

#include "common.h"

#include "Reaction.h"
#include "ChemRNode.h"
#include "Chemistry/ChemSim.h"

namespace medyan {

//FORWARD DECLARATIONS
class RNodeGillespie;
class ChemGillespieImpl;
    
/// Used by ChemGillespieImpl to implement the cached version of the Gillespie algorithm.

/*! RNodeGillespie manages a single chemical reaction within the Gillespie algorithm. 
 *  When the propensity drops to zero, the RNodeGillespie can execute the 
 *  passivateReaction() method. Alternatively, passivated RNodeGillespie can be
 *  activated via activateReaction(). The main part of the Gillespie algoritm is 
 *  implemented in the makeStep() method.
 */
class RNodeGillespie : public RNode {
public:
    /// Ctor:
    /// @param *r is the Reaction object corresponding to this RNodeGillespie
    /// @param &chem_Gillespie is a refernce to ChemGillespieImpl object, which does
    /// the overall management of the Gillespie scheme (e.g. random distribution
    /// generators, etc.)
    RNodeGillespie(ReactionBase *r, ChemGillespieImpl &chem_Gillespie);
    
    /// Copying is not allowed
    RNodeGillespie(const RNodeGillespie& rhs) = delete;
    
    /// Assignment is not allowed
    RNodeGillespie& operator=(RNodeGillespie &rhs) = delete;
    
    /// Dtor: The RNode pointer of the tracked Reaction object is set to nullptr
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~RNodeGillespie() noexcept;
            
    /// Returns a pointer to the Reaction which corresponds to this RNodeGillespie.
    ReactionBase* getReaction() const {return _react;};
    
    /// Return the currently stored propensity, "a", for this Reaction.
    /// @note The propensity is not recomputed in this method, so it potentially
    /// can be out of sync.
    floatingpoint getPropensity() const {return _a;}
    
    /// Set the propensity, "a", for this Reaction.
    void setA(floatingpoint a) {_a=a;}
    
    /// Return the propensity, "a", associated with the penultimate step
    /// of this Reaction.
    /// @note The propensity is not recomputed in this method, so it potentially
    /// can be out of sync.
    floatingpoint getPenultStepPropensity() const {return _a_prev;}

    /// Set the propensity, "a", associated with the penultimate step of this Reaction.
    void setPenultA(floatingpoint a_prev) {_a_prev=a_prev;}
    
    /// (Re)Compute and return the propensity associated with this Reaction.
    /// Remembers the penultimate step propensity as well
    floatingpoint reComputePropensity() {
        _a_prev=_a;
        _a=_react->computePropensity();
        return _a;
    }
    
    /// Set the the penultimate step propensity to zero and compute
    /// the current propensity.
    void reset() {
        _a_prev = 0;
        _a = _react->computePropensity();
    }
    
    /// This method calls the corresponding Reaction::makeStep() method of the underyling Reaction object
    void makeStep() {_react->makeStep();}
    
    /// Forwards the call to the similarly named method of ChemGillespieImpl
    virtual void activateReaction();
    
    /// Forwards the call to the similarly named method of ChemGillespieImpl
    virtual void passivateReaction();
    
    /// Forwards the call to the similarly named method of ChemGillespieImpl
    bool isPassivated() const {return _react->isPassivated();}
    
    /// Print information about this RNodeGillespie such as "a", "a_penult" and the
    /// Reaction which this RNodeGillespie tracks.
    void printSelf() const;
    
    /// Print the RNode objects which are dependents of this RNode (via the tracked Reaction object dependencies)
    void printDependents() const;
private:
    ChemGillespieImpl &_chem_Gillespie; ///< A reference to the ChemGillespieImpl
                                        ///< which containts the heap, random number
                                        ///< generators, etc.
    ReactionBase *_react; ///< The pointer to the associated Reaction object. The
                          ///<corresponding memory is not managed by RNodeGillespie.
    floatingpoint _a; ///< The propensity associated with the Reaction. It may be outdated
               ///< and may need to be recomputed if needed.
    floatingpoint _a_prev; ///< The propensity associated with the penultimate
                    ///< step of this Reaction.
};



/// Implements a slightly optimized version of the Gillespie algorithm.

/*! ChemGillespieImpl manages the Gillespie algorithm at the level of the network of 
 *  reactions. Reaction objects can be added and removed from the ChemGillespieImpl 
 *  instance. The propensities of all Reactions are cached, and they are recomputed 
 *  only when the copy number of associated reactant species gets changed (due to the
 *  currently occurring Reaction).
 *  @note The algorithm used by this class relies on tracking dependent
 *  Reactions, so TRACK_DEPENDENTS must be defined. The algorithm can work
 *  both when TRACK_ZERO_COPY_N and TRACK_UPPER_COPY_N are either defined or
 *  undefined. In the former case, Reaction objects with zero propensities
 *  stop being treated as dependents until their propensities become nonzero again.
 */
class ChemGillespieImpl : public ChemSim {
public:
    /// Ctor: Seeds the random number generator, sets global time to 0.0
    ///and the number of reactions to 0
    ChemGillespieImpl() :
        _uniform_distr(), _a_total(0),_n_reacts(0) { resetTime(); }
    
    /// Copying is not allowed
    ChemGillespieImpl(const ChemGillespieImpl &rhs) = delete;
    
    /// Assignment is not allowed
    ChemGillespieImpl& operator=(ChemGillespieImpl &rhs) = delete;
    
    ///Dtor: The reaction network is cleared. The RNodeGillespie objects will be
    /// destructed, but Reaction objects will stay intact.
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~ChemGillespieImpl() noexcept;
    
    /// Return the number of reactions in the network.
    size_t getSize() const {return _map_rnodes.size();}
    
    /// Return the current global time (as defined in the Gillespie algorithm)
    floatingpoint getTime() const {return _t;}
    
    /// Sets the global time to 0.0
    void resetTime() {_t=0.0; syncGlobalTime(); }
    
    /// Sets global time variable to ChemGillespieImpl's global time
    void syncGlobalTime() {global_time=_t; }
            
    /// Add ReactionBase *r to the network
    virtual void addReaction(ReactionBase *r);
    
    /// Remove ReactionBase *r from the network
    virtual void removeReaction(ReactionBase *r);

    //sets global time to restart time when called.
    virtual void initializerestart(floatingpoint restarttime);
    
    /// Unconditionally compute the total propensity associated with the network.
    floatingpoint computeTotalA();
    
    /// Returns a random time tau, drawn from the exponential distribution,
    /// with the propensity given by a.
    floatingpoint generateTau(floatingpoint a);
    
    /// Returns a random number between 0 and 1, drawn from the uniform distribution
   floatingpoint generateUniform();
    
    /// This function iterates over all RNodeGillespie objects in the network,
    /// activating all Reaction objects and calling reset().
    /// The total propentsity for the network is computed.
    /// @note This method needs to be called before calling run(...).
    /// @note If somewhere in the middle of simulaiton initialize() is called, it will
    /// be analogous to starting the simulation from scratch, except with the Species
    /// copy numbers given at that moment in time. The global time is reset to zero
    /// again.
    virtual void initialize();
    
    /// This method runs the Gillespie algorithm for the given amount of time.
    /// @return true if successful.
    virtual bool run(floatingpoint time) {
        
        floatingpoint endTime = _t + time;
        
        while(_t < endTime) {
            bool success = makeStep(endTime);
            if(!success)
                return false;
        }
        return true;
    }
    
    /// This method runs the Gillespie algorithm for the given amount of reaction steps.
    /// @return true if successful.
    virtual bool runSteps(int steps) {
        
        for(int i = 0; i < steps; i++) {
            
            bool success = makeStep();
            if(!success)
                return false;
        }
        return true;
    }
    
    /// This method is used to track the change in the total propensity of the network
    /// as the previously passivated ReactionBase *r has become activated
    void activateReaction(ReactionBase *r);
    
    /// This method is used to track the change in the total propensity of the network
    /// as the ReactionBase *r has become passivated
    void passivateReaction(ReactionBase *r);
    
    /// Prints all RNodes in the reaction network
    virtual void printReactions() const;

    /// Cross checks all reactions in the network for firing time.
    virtual bool crosschecktau() const {
        log::warn("Cannot check for tau in reactions in ChemGillespieImpl.h");
        return true;
    };
    
private:

    /// This subroutine, along with with passivateReaction() and activateReaction()
    /// implements a cached version of the Gillespie algorithm.
    /// First Tau the time to the next reaction event is sampled. 
    /// If it is after endTime, then time is set to endTime and Returns true. 
    /// Otherwise a uniform random number is generated to select which
    /// Reaction has occurred. Instead of computing the total propensity of the network
    /// from scratch, the cached value is being modified as Reaction events occur.
    /// Returns true if successful, false if endTime==inf and there are no reactions.
    bool makeStep(floatingpoint endTime = std::numeric_limits<floatingpoint>::infinity());

    //sets glocal time to specified value. To be used only during restart.
    void setTime(floatingpoint timepoint){ _t=timepoint; syncGlobalTime();}
private:
    #ifdef DEBUGCONSTANTSEED
    map<ReactionBase*, unique_ptr<RNodeGillespie>> _map_rnodes;
    #else
    unordered_map<ReactionBase*, unique_ptr<RNodeGillespie>> _map_rnodes;
	#endif
    ///< The database of RNodeGillespie
                                                                          ///< objects, representing the reaction network
    exponential_distribution<floatingpoint> _exp_distr; ///< Adaptor for the exponential distribution
    uniform_real_distribution<floatingpoint> _uniform_distr;
    floatingpoint _t; ///< global time
    floatingpoint _a_total;
    size_t _n_reacts; ///< number of reactions in the network
};

} // namespace medyan

#endif


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

/// \section nrm_sec Next Reaction Method

/// The algorithm implemented here is based on the following reference:
/// ** Michael A. Gibson, and Jehoshua Bruck J. Phys. Chem. A, 2000,
/// 104 (9), 1876-1889 **

#ifndef MEDYAN_ChemNRMImpl_h
#define MEDYAN_ChemNRMImpl_h

#include <vector>
#include <random>

#include <boost/heap/pairing_heap.hpp>

#include "common.h"

#include "Reaction.h"
#include "Chemistry/ChemSim.h"
#include "ChemRNode.h"

#ifdef BOOST_MEM_POOL
#include <boost/pool/pool.hpp>
#include <boost/pool/pool_alloc.hpp>
#endif

#define BOOST_POOL_MEM_PQNODE
#define BOOST_POOL_MEM_RNODENRM
#define BOOST_POOL_MEM_HEAP_ELEMENT

namespace medyan {

//FORWARD DECLARATIONS
class PQNode;
class RNodeNRM;
class ChemNRMImpl;

#if defined BOOST_MEM_POOL && defined BOOST_POOL_MEM_HEAP_ELEMENT
typedef boost::fast_pool_allocator<PQNode> fast_allocator_t;
typedef boost::heap::pairing_heap<PQNode,
        boost::heap::allocator<fast_allocator_t>> boost_heap;
typedef boost::heap::pairing_heap<PQNode,
        boost::heap::allocator<fast_allocator_t>>::handle_type handle_t;
#else
typedef boost::heap::pairing_heap<PQNode> boost_heap;
typedef boost::heap::pairing_heap<PQNode>::handle_type handle_t;
#endif
    
/// Priority Queue Node. It is stored as an element of a heap, such as
/// boost::heap::pairing_heap<PQNode>. There will be an associated heap handle which
/// can be used to dynamically access PQNode's fields
/// (e.g. boost::heap::pairing_heap<PQNode>::handle_type)

/*! 
 *  This is a simple structure which holds a pointer to RNodeNRM and stores the last 
 *  computed tau assoicated with this reaction. On some occasions tau needs to be 
 *  recomputed, for example when another reaction took place which affects this 
 *  reaction. handle_t of the corresponding heap element is used to get access to tau 
 *  (and RNodeNRM).
*/    
class PQNode {
public:
    /// Ctor
    /// @param *rnode is a pointer to the RNodeNRM instance which this PQNode tracks
    /// (no ownership - make sure it is not null)
    /// @note tau is set to infinity in the constructor
    PQNode(RNodeNRM *rnode) : _rn(rnode), _tau (numeric_limits<floatingpoint>::infinity()) {}

    /// Dtor: we only reassign _rn to null, the actual pointer needs to be deleted
    /// somewhere else to avoid memory leaks
    ~PQNode(){_rn=nullptr;}
    
    /// The PQ heap requires a comparision operator for ordering elements within the
    /// heap.
    /// @note In a bit of hack, the less operator is actually defined via real greater
    /// comparison of tau's of two corresponding PQNode objects. We do this so the top
    /// PQ node has the smallest tau and not largest.
    inline bool operator<(PQNode const &rhs) const{
        return _tau > rhs._tau;
    }
    
#ifdef BOOST_MEM_POOL
#ifdef BOOST_POOL_MEM_PQNODE
    /// Advanced memory management
    void* operator new(size_t size);
    
    void operator delete(void* ptr) noexcept;
#endif
#endif
private: 
    RNodeNRM *_rn; ///< Pointer to the reaction node (RNodeNRM) which this PQNode
                   ///< represents (or tracks)
    floatingpoint _tau;   ///< tau for this reaction for the Gibson-Bruck NRM algoritm
private:
    // Think of PQNode as a simple C-like structure (i.e. no methods),
    // but that is private to ChemNRMImpl and RNodeNRM
    friend class ChemNRMImpl;
    friend class RNodeNRM;
};

    
/// Reaction Node for the Next Reaction Method.

/*! 
 *  RNodeNRM manages a single chemical reaction within the NRM algorithm. It has a 
 *  pointer to the PQ element containing the Reaction via a handle_t object (and hence 
 *  can modify both the corresponding PQNode, such as PQNode's tau or the underlying 
 *  Reaction instance). RNodeNRM can recompute tau if needed and has auxiliary method
 *  for computing reaction's propensity. When the propensity drops to zero, the RNodeNRM 
 *  can execute the passivateReaction() method. Alternatively, passivated RNodeNRM can 
 *  be activated via activateReaction(). The main part of the NRM algoritm is
 *  implemented in the makeStep() method. 
 */
class RNodeNRM : public RNode {
public:
    /// Ctor: 
    /// @param *r is the Reaction object corresponding to this RNodeNRM
    /// @param &chem_nrm is a refernce to ChemNRMImpl object, which does the overall
    /// management of the NRM scheme (e.g. it gives acces to the heap itself, random
    /// distribution generators, etc.)
    RNodeNRM(ReactionBase *r, ChemNRMImpl &chem_nrm);
    
    /// Copying is not allowed
    RNodeNRM(const RNodeNRM& rhs) = delete;

    /// Assignment is not allowed
    RNodeNRM& operator=(RNodeNRM &rhs) = delete;
    
    /// Dtor: 1) Erases the corresponding PQNode element in the heap via the handle;
    /// 2) The RNode pointer of the tracked Reaction object is set to nullptr
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~RNodeNRM() noexcept;
    
    /// This methods recomputes the reaction propensity based on current coefficients
    /// and the rate, and then obtains a corresponding new tau drawn from an exponential
    /// distribution. @note This methods modifies tau which is stored in the associated
    /// PQNode object. However, it does not update the heap - this needs to be done
    /// seprately.
    void generateNewRandTau();
    
    
    /// This method reuses a previous tau to generate a new tau based on a change
    /// in reaction propensity.@note This methods modifies tau which is stored in the associated
    /// PQNode object. However, it does not update the heap - this needs to be done seperately.
    void generateNewTau();
    
    /// Returns a pointer to the Reaction which corresponds to this RNodeNRM.
    inline ReactionBase* getReaction() const {return _react;};
    
    /// The heap is updated only with respect to the specific PQNode element which
    /// presumably was modified (e.g. via generateNewRandTau()).
    void updateHeap();
    
    /// Returns the tau from the the associated PQNode element
    inline floatingpoint getTau() const {return (*_handle)._tau;}
    
    /// Sets the tau of the the associated PQNode element. Does not update the heap.
    inline void setTau(floatingpoint tau) {(*_handle)._tau=tau;}
    
    /// Returns a handle to the associated PQNode element, which can be used to access
    /// tau, for example.
    inline handle_t& getHandle() {return _handle;}
    
    /// Return the currently stored propensity, "a", for this Reaction.
    /// @note The propensity is not recomputed in this method, so it potentially can be
    /// out of sync.
    inline floatingpoint getPropensity() const {return _a;}
    
    /// (Re)Compute and return the propensity associated with this Reaction.
    inline floatingpoint reComputePropensity() {_a=_react->computePropensity(); return _a;}

    /// Compute and return the product of reactant copy numbers. This method can be used
    /// to quickly check whether this reaction needs to be passivated, if the returned
    /// result is zero.
    inline int getProductOfReactants () {return _react->getProductOfReactants ();};
    
    /// This method calls the corresponding Reaction::makeStep() method of the
    /// underyling Reaction object
    inline void makeStep() {_react->makeStep();}
    
    /// When this method is called, a new tau is computed and the corresponding PQNode
    /// element is updated in the heap.
    /// @note This method does not activate the Reaction itself, but instead only deals
    /// with the activation process related to the corresponding PQNode element.
    virtual void activateReaction();
    
    /// When this method is called, reaction's tau is set to infinity, the propensity is
    /// set to 0, and the corresponding PQNode element is updated in the heap.
    /// @note This method does not passivate the Reaction itself, but instead only deals
    /// with the activation process related to the corresponding PQNode element.
    virtual void passivateReaction();
    
    /// Return true if the Reaction is currently passivated
    inline bool isPassivated() const {return _react->isPassivated();}
    
    /// Print information about this RNodeNRM such as tau, a and the Reaction which this
    /// RNodeNRM tracks.
    void printSelf() const;
    
    /// Print the RNode objects which are dependents of this RNode (via the tracked
    /// Reaction object dependencies)
    void printDependents() const;
    
#ifdef BOOST_MEM_POOL
#ifdef BOOST_POOL_MEM_RNODENRM
    /// Advanced memory management
    void* operator new(size_t size);
    
    void operator delete(void* ptr) noexcept;
#endif
#endif
private:
    ChemNRMImpl &_chem_nrm; ///< A reference to the ChemNRMImpl which containts the
                            ///<heap, random number generators, etc.
    handle_t _handle; ///< The handle to the associated PQNode element in the heap.
    ReactionBase *_react; ///< The pointer to the associated Reaction object. The
                          ///< corresponding memory is not managed by RNodeNRM.
    floatingpoint _a; ///< The propensity associated with the Reaction.
               ///< It may be outdated and may need to be recomputed if needed.
};


    
/// The chemical NRM implementation. 

/*! 
 *  ChemNRMImpl manages the NRM algorithm at the level of the network of reactions. In 
 *  particular, this class contains the NRM heap and the exponential random number 
 *  generator. Reaction objects can be added and removed from the ChemNRMImpl instance.
 */
class ChemNRMImpl : public ChemSim {
public:
    /// Ctor: Seeds the random number generator, sets global time to 0.0 and the number
    /// of reactions to 0
    ChemNRMImpl() { resetTime(); }
    
    /// Copying is not allowed
    ChemNRMImpl(const ChemNRMImpl &rhs) = delete;
    
    /// Assignment is not allowed
    ChemNRMImpl& operator=(ChemNRMImpl &rhs) = delete;
    
    ///Dtor: The reaction network is cleared. The RNodeNRM objects will be destructed,
    /// but Reaction objects will stay intact. When the RNodeNRM objects are destructed,
    /// they in turn destruct the corresponding PQNode element, setting the RNode
    /// pointer of the Reaction object to null. At the end, the heap itself will go out
    /// of scope.
    virtual ~ChemNRMImpl();
    
    /// Return the number of reactions in the network.
    inline auto getSize() const { return _map_rnodes.size(); }
    
    /// Return the current global time (as defined in the NRM algorithm)
    inline floatingpoint getTime() const {return _t;}
    
    /// Sets global time to 0.0
    inline void resetTime() {_t=0.0; syncGlobalTime(); }
    
    /// Sets global time variable to ChemNRMImpl's global time
    inline void syncGlobalTime() {global_time=_t;}
    
    /// Return a pointer to the heap
    inline boost_heap* getHeap() {return &_heap;}
    
    /// Add Reaction *r to the network
    virtual void addReaction(ReactionBase *r);
    
    /// Remove Reaction *r from the network
    virtual void removeReaction(ReactionBase *r);

    //sets global time to restart time when called.
    virtual void initializerestart(floatingpoint restarttime);
    
    /// A pure function (without sideeffects), which returns a random time tau, drawn
    /// from the exponential distribution, with the propensity given by a.
    floatingpoint generateTau(floatingpoint a);
        
    /// This function iterates over all RNodeNRM objects in the network, generating new
    /// tau-s for each case and subsequently updating the heap. It needs to be called
    /// before calling run(...).
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
            if(!success) {
                return false;
            }
        }
        return true;
    }
    
    /// Prints all RNodes in the reaction network
    virtual void printReactions() const;

    /// Cross checks all reactions in the network for firing time.
    virtual bool crosschecktau() const;
    
private:
    /// This is a somewhat complex subroutine which implements the main part of the
    /// Gibson-Bruck NRM algoritm. See the implementation for details. After this method
    /// returns, roughly the following will have happened:
    /// 1) If the lowest reaction time is after endTime,
    /// the time is set to endTime and returns true. Otherwise,
    /// 2) The Reaction corresponding to the lowest tau RNodeNRM is
    /// executed and the corresponding Species copy numbers are changed
    /// 3) A new tau is computed from this Reaction and the corresponding PQNode element
    /// is updated in the heap
    /// 4) The other affected Reaction objects are found, their taus are recomputed and
    /// corresponding PQNode elements are updated in the heap.
    /// 5) For the Reaction and associated Species signals are emitted, if these objects
    /// broadcast signals upon change.
    /// Returns true if successful, and false if the heap is exchausted and there no
    /// more reactions to fire and endTime==inf.
    bool makeStep(floatingpoint endTime = std::numeric_limits<floatingpoint>::infinity());

    //sets glocal time to specified value. To be used only during restart.
    void setTime(floatingpoint timepoint){ _t=timepoint; syncGlobalTime();}
private:
    unordered_map<ReactionBase*, unique_ptr<RNodeNRM>> _map_rnodes; ///< The database of RNodeNRM objects,
                                                                    ///< representing the reaction network
    boost_heap _heap; ///< A priority queue for the NRM algorithm,
                      ///< containing PQNode elements
    exponential_distribution<floatingpoint> _exp_distr; ///< Adaptor for the exponential distribution
    floatingpoint _t; ///< global time
};

} // namespace medyan
    
#endif


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

#ifndef MEDYAN_RSpecies_h
#define MEDYAN_RSpecies_h

#include <functional>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <cassert>
#include <stdexcept>

#include "common.h"
#include "SysParams.h"

namespace medyan {
///Enumeration for RSpecies types
#undef CONST // Fuck <windows.h>
enum class RSpeciesType {
    REG, AVG, CONST
};

//FORWARD DECLARATIONS
class Species;
class RSpecies;
class ReactionBase;

/// vr stands for vector of Reactions
typedef vector<ReactionBase*>::iterator vr_iterator; 
typedef vector<ReactionBase*>::const_iterator vr_const_iterator; 

/// vsp stands for vector of RSpecies
typedef vector<RSpecies*>::iterator vrsp_iterator; 
typedef vector<RSpecies*>::const_iterator vrsp_const_iterator; 

/// Represents the reactive aspect of chemical molecules. It tracks their copy
/// number and can be used in [Reactions](@ref Reaction).

/*!  This abstract class represents the reactivity of chemical species. The name RSpecies stems
 *   from reacting species. RSpecies tracks the copy number of molecules and the 
 *   [Reactions](@ref Reaction) in which it is involed (@see Reaction).
 *   @note Each intantiation of RSpecies is unique, and hence, cannot be neither copied
 *   nor moved (C++11). This has significant implications - e.g., RSpecies cannot be 
 *   used in vector<RSpecies>. Instead, one should use either vector<RSpecies*> if not 
 *   owning the RSpecies pointers, or vector<unique_ptr<RSpecies>> if owning the 
 *   RSpecies pointers. A special allocator can be written such that dynamically 
 *   allocated RSpecies (through new), are arranged contigiously in memory.
 */
class RSpecies {
    friend class Species;
    friend class CMonomer;
    /// Reactions calls addAsReactant(), removeAsReactant() - which other classes should not call
protected: //Variables
    vector<ReactionBase *> _as_reactants = {}; ///< a vector of [Reactions]
                                               ///< (@ref Reaction) where this RSpecies
                                               ///< is a Reactant
    vector<ReactionBase *> _as_products = {};  ///< a vector of [Reactions]
                                               ///< (@ref Reaction) where this RSpecies
                                               ///< is a Product
    Species& _species; ///< reference to the **parent** Species object
    species_copy_t _n; ///< Current copy number of this RSpecies
#ifdef TRACK_UPPER_COPY_N
    species_copy_t _ulim; ///< Upper limit for the copy number, afterwards all
                          ///< reactions leading to further accum. are turned off
#endif

    // The callbacks associated with change of n of this RSpecies (usually when a single step of this Reaction occurs).
    std::vector< std::function< void(RSpecies*, int) >> callbacks_;


    RSpeciesType _type; ///< The RSpecies type
    
//CONSTRUCTORS
    
    /// Constructor
    /// @param parent - the Species object to which this RSpecies belongs
    /// @param n - copy number
    RSpecies (Species &parent, species_copy_t n, species_copy_t ulim)
    : _species(parent), _n(n) {
#ifdef TRACK_UPPER_COPY_N
        _ulim = ulim;
#endif
    }
    /// deleted copy constructor - each RSpecies is uniquely created by the parent
    /// Species
    RSpecies(const RSpecies &r) = delete;
    
    /// deleted move constructor - each RSpecies is uniquely created by the parent
    /// Species
    RSpecies(RSpecies &&r) = delete;
    
    /// deleted assignment operator - each RSpecies is uniquely created by the parent
    /// Species
    RSpecies& operator=(RSpecies&) = delete;
    
public:
           
    /// Sets the copy number for this RSpecies. 
    /// @param n should be a non-negative number, but no checking is done in run time
    /// @note The operation does not emit any signals about the copy number change.
    void setN(species_copy_t n) {
        _n=n;
    }

#ifdef TRACK_UPPER_COPY_N        
    /// Sets the upper limit for the copy number of this RSpecies.
    /// @param ulim should be a non-negative number, but no checking is done in run time
    void setUpperLimitForN(species_copy_t ulim) {_ulim=ulim;}
    
    /// Return the upper limit for the copy number of this RSpecies
    species_copy_t getUpperLimitForN() const {return _ulim;}
#endif
    
    // \internal This methods is called by the Reaction class during construction
    // of the Reaction where this RSpecies is involved as a Reactant
    void addAsReactant(ReactionBase *r){_as_reactants.push_back(r);}
    
    // \internal This methods is called by the Reaction class during construction
    // of the Reaction where this RSpecies is involved as a Product    
    void addAsProduct(ReactionBase *r){_as_products.push_back(r);}
    
    // \internal This method is called by the Reaction class during destruction
    // of the Reaction where this RSpecies is involved as a Reactant
    void removeAsReactant(const ReactionBase *r) {
        
        if(!_as_reactants.empty()) {
            auto rxit = find(_as_reactants.begin(),_as_reactants.end(),r);
            if(rxit!=_as_reactants.end()){
                _as_reactants.erase(rxit);
            }
            else {
                
                cout << "Did not find the as_reactant" << endl;
                
            }
        }
    }
    // \internal This method is called by the Reaction class during destruction
    // of the Reaction where this RSpecies is involved as a Product
    void removeAsProduct(const ReactionBase *r) {
        
        if (!_as_products.empty()) {
            auto rxit = find(_as_products.begin(),_as_products.end(),r);
            if(rxit!=_as_products.end()){
                _as_products.erase(rxit);
            }
            else {
                cout << "Did not find the as_product" << endl;
            }
        }
        
    }
    /// \internal Attempts to activate previously passivated [Reactions](@ref Reaction)
    /// where this RSpecies is involved as a Reactant. Usually, the Reaction was first
    /// passivated, for example if the RSpecies copy number of one of the reactants
    /// dropeed to zero. This attempt may not succeed if there are still other reactants
    /// in the same Reaction with zero copy count.
    void activateAssocReactantReactions();
    
    /// \internal Attempts to activate previously passivated [Reactions](@ref Reaction)
    /// where this RSpecies is involved as a Product. Usually, the Reaction was first
    /// passivated, for example if the RSpecies copy number of one of the reactants went
    /// up to max_ulim.
    void activateAssocProductReactions();
    
    /// \internal Passivates all [Reactions](@ref Reaction) where this RSpecies is among
    /// the reactants.
    void passivateAssocReactantReactions();
    
    /// \internal Passivates all [Reactions](@ref Reaction) where this RSpecies is among
    /// the products.
    void passivateAssocProductReactions();
           
    /// Clear all callbacks associated with this RSpecies; all associated slots will be
    /// destroyed
    void clearSignaling();
    
public:
    /// It is required that all [Reactions](@ref Reaction) associated with this
    /// RSpecies are destructed before this RSpecies is destructed. Most of the time,
    /// this will occur naturally. If not, an assertion will ungracefully terminate the
    /// program.
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~RSpecies() noexcept;
    
    /// Broadcasts signal indicating that the copy number of this RSpecies has changed
    /// This method should usually called by the code which runs the chemical dynamics
    /// (i.e. Gillespie-like algorithm)
    void emitSignal(int delta) {
        for(auto& callback : callbacks_) {
            callback(this, delta);
        }
    }
    
    /// return parent Species as a reference
    inline Species& getSpecies() {return _species;}
    
    /// return parent Species as a const reference
    inline const Species& getSpecies() const {return _species;}
    
    /// Return the full name of this Species in a string format (e.g. "Arp2/3{Bulk}")
    string getFullName() const;
            
    /// Return vector<ReactionBase *>, which contains pointers to all [Reactions](@ref
    /// Reaction) where this RSpecies is involved as a Reactant
    inline vector<ReactionBase *>& reactantReactions(){return _as_reactants;}
    
    /// Return vector<ReactionBase *>, which contains pointers to all [Reactions](@ref
    /// Reaction) where this RSpecies is involved as a Product
    inline vector<ReactionBase *>& productReactions(){return _as_products;}
    
    /// Return vector<ReactionBase *>::iterator, which points to the beginning of all 
    /// [Reactions](@ref Reaction) where this RSpecies is involved as a Reactant
    inline vr_iterator beginReactantReactions() {return _as_reactants.begin();}
    
    /// Return vector<ReactionBase *>::iterator, which points to the beginning of all 
    /// [Reactions](@ref Reaction) where this RSpecies is involved as a Product
    inline vr_iterator beginProductReactions() {return _as_products.begin();}
    
    /// Return vector<ReactionBase *>::iterator, which points to the end of all 
    /// [Reactions](@ref Reaction) where this RSpecies is involved as a Reactant    
    inline vr_iterator endReactantReactions() {return _as_reactants.end();}
    
    /// Return vector<ReactionBase *>::iterator, which points to the end of all 
    /// [Reactions](@ref Reaction) where this RSpecies is involved as a Product
    inline vr_iterator endProductReactions() {return _as_products.end();}
    
    /// Increases the copy number of this RSpecies
    virtual void up() = 0;
    /// Decreases the copy number of this RSpecies
    virtual void down() = 0;
    
    /// Return the effective copy number of this RSpecies
    virtual float getN() const = 0;
    /// Return the true copy number of this RSpecies
    inline species_copy_t getTrueN() const {return _n;}
    
};

/// A RSpecies implementation that changes copy number regularly.
/*!
 * The RSpeciesReg class represents a typical RSpecies. Upon calling up() or down(),
 * the copy number is incremented normally. getN() and getTrueN() both return the 
 * true copy number tracked in the base class.
 *
 * The RSpeciesReg can be used for any type of Species in the chemical system.
 */
class RSpeciesReg : public RSpecies {
    
public:
    /// Constructors
    RSpeciesReg (Species &parent, species_copy_t n=0, species_copy_t ulim=max_ulim)
        : RSpecies(parent, n, ulim) {
    
        _type = RSpeciesType::REG;
    }
    /// deleted copy constructor - each RSpeciesRegular is uniquely created by the parent
    /// Species
    RSpeciesReg(const RSpeciesReg &r) = delete;
    
    /// deleted move constructor - each RSpeciesRegular is uniquely created by the parent
    /// Species
    RSpeciesReg(RSpeciesReg &&r) = delete;
    
    /// deleted assignment operator - each RSpeciesRegular is uniquely created by the
    /// parent Species
    RSpeciesReg& operator=(RSpeciesReg&) = delete;
    
#ifdef BOOST_MEM_POOL
    /// Advanced memory management
    void* operator new(size_t size);
    
    void operator delete(void* ptr) noexcept;
#endif
    
    /// If the copy number changes from 0 to 1, calls a
    /// "callback"-like method to activate previously passivated [Reactions](@ref
    /// Reaction), where this RSpecies is a Reactant.
    /// Also emits a signal with the change in copy number if attached.
    virtual void up() {
        
        _n+=1;
#ifdef TRACK_ZERO_COPY_N
        if(_n==1)
            activateAssocReactantReactions();
#endif
#ifdef TRACK_UPPER_COPY_N
        if(_n==_ulim)
            passivateAssocProductReactions();
#endif
        
        emitSignal(+1);
    }
    
    /// If the copy number changes from 1 to 0, calls a
    /// "callback"-like method to passivate previously activated [Reactions](@ref
    /// Reaction), where this RSpecies is a Reactant.
    /// Also emits a signal with the change in copy number if attached.
    virtual void down() {
        
#ifdef TRACK_UPPER_COPY_N
        species_copy_t prev_n = _n;
#endif
        _n-=1;
//        if(_n > 100000){
//            std::cout << "At time = " << tau() << endl;
//            std::cout << "Wrong CopyNum = " << _n << " for species: " << getFullName() << endl;
//        }
#ifdef TRACK_ZERO_COPY_N
        if(_n == 0)
            passivateAssocReactantReactions();
#endif
#ifdef TRACK_UPPER_COPY_N
        if(prev_n == _ulim)
            activateAssocProductReactions();
#endif
        
        emitSignal(-1);
    }
    /// Return the true copy number
    virtual float getN() const {return (float)_n;}
};


/// A constant RSpecies whose copy number does not change.
/*!
 *  The RSpeciesConst class represents a RSpecies that does not change copy number upon 
 *  reacting. This is used for species whose copy numbers are constant throughout the 
 *  duration of the chemical simulation.
 *
 *  The RSpeciesAvg can only be used for SpeciesBulk or SpeciesDiffusing.
 */
class RSpeciesConst : public RSpecies {
    
public:
    /// Constructors
    RSpeciesConst (Species &parent, species_copy_t n=0, species_copy_t ulim=max_ulim)
        : RSpecies(parent, n, ulim){
        
        _type = RSpeciesType::CONST;
    }
    /// deleted copy constructor - each RSpeciesConst is uniquely created by the parent
    /// Species
    RSpeciesConst(const RSpeciesConst &r) = delete;
    
    /// deleted move constructor - each RSpeciesConst is uniquely created by the parent
    /// Species
    RSpeciesConst(RSpeciesConst &&r) = delete;
    
    /// deleted assignment operator - each RSpeciesConst is uniquely created by the
    /// parent Species
    RSpeciesConst& operator=(RSpeciesConst&) = delete;
    
#ifdef BOOST_MEM_POOL
    /// Advanced memory management
    void* operator new(size_t size);
    
    void operator delete(void* ptr) noexcept;
#endif
    
    //@{
    /// In constant species, do nothing. Copy numbers do not change.
    /// Emits a signal with the zero change in copy number if attached.
    virtual void up() {
        emitSignal(0);
    }
    virtual void down() {
        emitSignal(0);
    }
    //@}
    /// Return the true copy number
    virtual float getN() const {return (float)_n;}
};

/// An average RSpecies that tracks the average copy number over a number of events.
/*!
 *  The RSpeciesAvg class represents a RSpecies that, while updating the true copy
 *  number of the species, also tracks the average copy number over a specified number
 *  of events. up() and down() update the true copy number while adding to a running
 *  average of copy number events. When the number of events to average (named _numEvents)
 *  is reached, the RSpecies computes the new average value and updates accordingly.
 *
 *  The RSpeciesAvg can only be used for SpeciesBulk or SpeciesDiffusing.
 *
 *  @note - Still buggy for lower concentration averaging. Use with caution. Since this
 *          only produces speedups for larger diffusive concentrations, this choice would
 *          be sub-optimal anyway.
 *
 */
class RSpeciesAvg : public RSpecies {
    
friend class Species;
    
private:
    int _numEvents;      ///< The number of events to generate an average copy number
    int _eventCount = 0; ///< Tracking the number of events since last average
    
    floatingpoint _localTau;  ///< The time since the last event occured
    floatingpoint _firstTau;  ///< Time since an average update occured
    
    floatingpoint _average;   ///< The current average value
    
    bool _newAvg = false; ///< If we just calculated a new average. Used by Reaction to mark the
                          ///< updating of its dependencies accordingly.
    
    /// A map of copy number and time (representing the time the species has been at this
    /// copy number). These values will be collected until a new average is needed, and the
    /// map will be used to compute an average and then reset.
    unordered_map<species_copy_t, floatingpoint> _copyNumbers;
    
    ///Compute the time average of the currently tracked, previous copy numbers
    void computeAverageN() {
        
        floatingpoint totalTau = tau() - _firstTau;
        floatingpoint average = 0;
        
        //average values
        for(auto it : _copyNumbers)
            average += it.first * it.second;
        
        average /= totalTau;
        
        //clear map and return
        _copyNumbers.clear();
        
        _average = average;
    }
    
public:
    /// Constructors
    /// @param numEvents - the number of copy number change events to use before
    /// computing a new average. When this number is reached, a new average is created
    /// using the _copyNumbers map (time averaged) and replaces the previous average.
    RSpeciesAvg (Species &parent, species_copy_t n=0, species_copy_t ulim=max_ulim)

        : RSpecies(parent, n, ulim) {
        
        _type = RSpeciesType::AVG;
            
        //set first average to n, first time update
        _average = n;
        _localTau = _firstTau = tau();
    }
    /// deleted copy constructor - each RSpeciesAvg is uniquely created by the parent
    /// Species
    RSpeciesAvg(const RSpeciesAvg &r) = delete;
    
    /// deleted move constructor - each RSpeciesAvg is uniquely created by the parent
    /// Species
    RSpeciesAvg(RSpeciesAvg &&r) = delete;
    
    /// deleted assignment operator - each RSpeciesAvg is uniquely created by the
    /// parent Species
    RSpeciesAvg& operator=(RSpeciesAvg&) = delete;
    
#ifdef BOOST_MEM_POOL
    /// Advanced memory management
    void* operator new(size_t size);
    
    void operator delete(void* ptr) noexcept;
#endif
    
    /// Whether we just calculated a new average
    bool newAverage() {return _newAvg;}
    
    /// Set the number of events to use averaging
    void setNumEvents(int numEvents) {_numEvents = numEvents;}
    
    /// Increase the true copy number.
    /// Add the old copy number to the running copy number map.
    /// If a new average is needed, compute it and reset accordingly.
    /// If the copy number changes from 1 to 0, calls a
    /// "callback"-like method to passivate previously activated [Reactions](@ref
    /// Reaction), where this RSpecies is a Reactant.
    /// Also emits a signal with the change in copy number if attached.
    virtual void up() {
        
        // if initialization, set avg to be true
        // copy number and return
        floatingpoint deltaTau = tau() - _localTau;
        
        if(areEqual(deltaTau, 0.0)) {
            
            _average = ++_n;
            return;
        }

        //add old n to map
        _copyNumbers[_n] += deltaTau;
        
        //compute a new avg if we need it
        if(++_eventCount > _numEvents) {
            computeAverageN();
            
            _eventCount = 0;
            _firstTau = tau();
            
            _newAvg = true;
        }
        else _newAvg = false;
        
        //increase copy number, reset local tau
        _n++; _localTau = tau();
        
#ifdef TRACK_ZERO_COPY_N
        if(_n==1)
            activateAssocReactantReactions();
#endif
#ifdef TRACK_UPPER_COPY_N
        if(_n==_ulim)
            passivateAssocProductReactions();
#endif
        
        emitSignal(+1);
    }
    
    /// Decrease the true copy number.
    /// Add the old copy number to the running copy number map.
    /// If a new average is needed, compute it and reset accordingly.
    /// If the copy number changes from 1 to 0, calls a
    /// "callback"-like method to passivate previously activated [Reactions](@ref
    /// Reaction), where this RSpecies is a Reactant.
    /// Also emits a signal with the change in copy number if attached.
    virtual void down() {
        
#ifdef TRACK_UPPER_COPY_N
        species_copy_t prev_n = _n;
#endif
        floatingpoint deltaTau = tau() - _localTau;
        
        //add old n to map
        _copyNumbers[_n] += deltaTau;
        
        //compute a new avg if we need it
        if(++_eventCount > _numEvents) {
            computeAverageN();
            
            //clear map and return
            _copyNumbers.clear();
            
            _eventCount = 0;
            _firstTau = tau();
            
            _newAvg = true;
        }
        else _newAvg = false;
        
        //decrease copy number, reset local tau
        _n--; _localTau = tau();
        
#ifdef TRACK_ZERO_COPY_N
        if(_n == 0)
            passivateAssocReactantReactions();
#endif
#ifdef TRACK_UPPER_COPY_N
        if(prev_n == _ulim)
            activateAssocProductReactions();
#endif
        
        emitSignal(-1);
    }
    
    /// Return the current average
    virtual float getN() const {return _average;}
};

/// Print self into an iostream
ostream& operator<<(ostream& os, const RSpecies& s);


/// A factory class to create RSpecies objects
class RSpeciesFactory {
    
public:
    ///Create an RSpecies object
    static std::unique_ptr< RSpecies > createRSpecies(
        Species &      parent,
        species_copy_t n = 0,
        species_copy_t ulim = max_ulim,
        RSpeciesType   t = RSpeciesType::REG
    ) {
        //create the appropriate rspecies
        //average
        if(t == RSpeciesType::AVG) return std::make_unique< RSpeciesAvg >(parent, n, ulim);
        //constant
        else if(t == RSpeciesType::CONST) return std::make_unique< RSpeciesConst >(parent, n, ulim);
        //regular
        else if(t == RSpeciesType::REG) return std::make_unique< RSpeciesReg >(parent, n, ulim);
        
        else {
            LOG(ERROR) << "Species type is not recognized";
            throw std::runtime_error("Species type is incorrect");
        }
    }
};

} // namespace medyan

#endif

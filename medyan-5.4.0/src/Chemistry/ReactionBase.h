
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

#ifndef MEDYAN_ReactionBase_h
#define MEDYAN_ReactionBase_h

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <set>
#include <unordered_set>
#include <utility>

#include "common.h"
#include "Chemistry/Species.h"

namespace medyan {
//FORWARD DECLARATIONS
class CBound;
class RNode;
class Composite;
class ReactionBase;
class SpeciesPtrContainerVector;

///Enumeration for type of reaction
enum class ReactionType {
    unspecified,
    REGULAR, DIFFUSION,
    emission, absorption,
    POLYMERIZATIONPLUSEND, POLYMERIZATIONMINUSEND,
    DEPOLYMERIZATIONPLUSEND, DEPOLYMERIZATIONMINUSEND,
    LINKERBINDING, MOTORBINDING, LINKERUNBINDING, MOTORUNBINDING,
    MOTORWALKINGFORWARD, MOTORWALKINGBACKWARD,
    AGING, FILAMENTCREATION, FILAMENTDESTRUCTION,
    BRANCHING, BRANCHUNBINDING, SEVERING,
    membraneDiffusion,
    membraneAdsorption, membraneDesorption,
};

constexpr const char* text(ReactionType t) {
    switch(t) {
        case ReactionType::REGULAR:                  return "regular";
        case ReactionType::DIFFUSION:                return "diffusion";
        case ReactionType::emission:                 return "emission";
        case ReactionType::absorption:               return "absorption";
        case ReactionType::POLYMERIZATIONPLUSEND:    return "polymerizationPlusEnd";
        case ReactionType::POLYMERIZATIONMINUSEND:   return "polymerizationMinusEnd";
        case ReactionType::DEPOLYMERIZATIONPLUSEND:  return "depolymerizationPlusEnd";
        case ReactionType::DEPOLYMERIZATIONMINUSEND: return "depolymerizationMinusEnd";
        case ReactionType::LINKERBINDING:            return "linkerBinding";
        case ReactionType::MOTORBINDING:             return "motorBinding";
        case ReactionType::LINKERUNBINDING:          return "linkerUnbinding";
        case ReactionType::MOTORUNBINDING:           return "motorUnbinding";
        case ReactionType::MOTORWALKINGFORWARD:      return "motorWalkingForward";
        case ReactionType::MOTORWALKINGBACKWARD:     return "motorWalkingBackward";
        case ReactionType::AGING:                    return "aging";
        case ReactionType::FILAMENTCREATION:         return "filamentCreation";
        case ReactionType::FILAMENTDESTRUCTION:      return "filamentDestruction";
        case ReactionType::BRANCHING:                return "branching";
        case ReactionType::BRANCHUNBINDING:          return "branchUnbinding";
        case ReactionType::SEVERING:                 return "severing";
        case ReactionType::membraneDiffusion:        return "membraneDiffusion";
        case ReactionType::membraneAdsorption:       return "membraneAdsorption";
        case ReactionType::membraneDesorption:       return "membraneDesorption";
        default:                                     return "unspecified";
    }
}


/// Represents an abstract interface for simple chemical reactions of the form
/// A + B -> C.

/*! ReactionBase provides an interface for managing a chemical reaction. It is an 
 *  abstract interface, so it cannot be directly instantiated, but other concrete 
 *  classes may be derived from it. ReactionBase may have a composite object as a 
 *  parent. A signaling interface may be used to make callbacks when some event, such 
 *  as a single reaction step, has been executed.
 *  The ReactionBase indicates a forward process only. For processes in both directions, 
 *  e.g. A <-> B, two ReactionBases objects need to be defined, corresponding to A->B 
 *  and B->A.
 *  A ReactionBase tracks other ReactionBase objects that are affected if this 
 *  ReactionBase is executed. A ReactionBase may be set up such that it "signals" when a 
 *  ReactionBase event happens, in which case the corresponding callbacks are called.
 */
class ReactionBase {
public:
    // The type of callback invoked when this reaction is fired.
    using CallbackType = std::function< void(ReactionBase*) >;

protected:
	#ifdef DEBUGCONSTANTSEED
	using dependentdatatype = unordered_set<ReactionBase*, HashbyId<ReactionBase*>,
            customEqualId<ReactionBase*>>;
	#else
    using dependentdatatype = unordered_set<ReactionBase*>; ///< Pointers to
    // ReactionBase objects that depend
                                               ///< on this ReactionBase being executed
	#endif
    dependentdatatype _dependents;
    RNode* _rnode; ///< A pointer to an RNode object which is used
                   ///< to implement a Gillespie-like algorithm (e.g. NRM)
    
    Composite *_parent; ///< A pointer to a Composite object to which
                        ///< this Reaction belongs
    
    FP _rate;      ///< the rate for this ReactionBase
    FP _rate_bare; ///< the bare rate for this ReactionBase (original rate)

    static size_t  _Idcounter;


    
    std::vector<CallbackType> callbacks_;

#if defined TRACK_ZERO_COPY_N || defined TRACK_UPPER_COPY_N
    bool _passivated; ///< Indicates whether the ReactionBase is currently passivated
#endif
    ReactionType _reactionType = ReactionType::unspecified; ///< Reaction type enumeration
    
    bool _isProtoCompartment = false;///< Reaction is in proto compartment
    ///< (Do not copy as a dependent, not in ChemSim)
    
    CBound* _cBound = nullptr; ///< CBound that is attached to this reaction
    
    
    float _gnum = 0.0;
    
    string _hrcdid = "DNT";
    
    float _linkRateForward = 0.0;
    
    float _linkRateBackward = 0.0;

    FP _volumeFrac; ///< Used in compartments to store volume fraction of the compartment
    int _rateVolumeDepExp; ///< Exponent of rate dependency on volume
    ///< Dependence on bulk properties are NOT considered currently

    size_t _Id = 0;

public:
    // Multiplicative factors used to update rate of a reaction.
    enum RateMulFactorType {
        // The volume factor is automatically adjusted by the volumeFrac parameter, and should not be set manually.
        VOLUMEFACTOR,
        // Factors that can be set manually.
        diffusionShape, mechanochemical, MOTORWALKCONSTRAINTFACTOR,
        RESTARTPHASESWITCH, MANUALRATECHANGEFACTOR1,
        // This must be the last item, indicating the total number of items of this enum.
        RATEMULFACTSIZE
    };
    std::array<FP, RATEMULFACTSIZE> _ratemulfactors;

    void setRateMulFactor(FP factor, RateMulFactorType type){
        //Call to this function should always be followed with call to updatePropensity

        if(factor == _ratemulfactors[type]) return;

        if(_ratemulfactors[type] == 0.0 || std::isinf(_ratemulfactors[type])) {
            _rate = _rate_bare;
            _ratemulfactors[type] = factor;
            for(unsigned i = 0; i < RateMulFactorType::RATEMULFACTSIZE; i++)
                _rate *= _ratemulfactors[i];
        } else {
            _rate = _rate * factor / _ratemulfactors[type];
            _ratemulfactors[type] = factor;
        }
    }
    auto getRateMulFactor(RateMulFactorType type) const {
        return _ratemulfactors[underlying(type)];
    }

    /// The main constructor:
    /// @param rate - the rate constant (full volume) for this ReactionBase
    ReactionBase (FP rate, bool isProtoCompartment, FP volumeFrac=1.0, int rateVolumeDepExp=0);
    
    /// No copying (including all derived classes)
    ReactionBase (const ReactionBase &rb) = delete;
    
    /// no assignment (including all derived classes)
    ReactionBase& operator=(ReactionBase &rb) = delete;
    
    /// Destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~ReactionBase() noexcept { }
    
    /// Copy this reaction using SpeciesPtrContainerVector &spcv as a source of
    /// analogous Species.
    /// @return the cloned ReactionBase pointer.
    /// @note the receving object should take care of the memory management of the
    /// returned pointer
    ReactionBase* clone(const SpeciesPtrContainerVector &spcv) {
        return cloneImpl(spcv);
    }
    
    /// (Private) implementation of the clone() method to be elaborated in derived classes
    virtual ReactionBase* cloneImpl(const SpeciesPtrContainerVector &spcv) = 0;
    
    /// Returns a pointer to the first element of the array<RSpecies*>. This pointer can
    /// be used to iterate over RSpecies* if necessary (also by relying on getM() and
    /// size() to determine the iteration limits). The corresponding array<RSpecies*> is
    /// defined by the derived classes.
    virtual RSpecies** rspecies() = 0;
    
    //aravind, June 30, 2016.
    vector<string> getreactantspecies(){
        vector<string> returnvector;
        for(int i=0;i<2;i++){
            returnvector.push_back((*(rspecies()+i))->getSpecies().getName());
        }
        return returnvector;
    }

    
    vector<string> getReactantSpecies() {
        vector<string> returnvector;
        for(auto i=0U;i<getM();++i){
            string name = (*(rspecies()+i))->getSpecies().getName();
            string namecut = name.substr(0,name.find("-",0));
            returnvector.push_back(namecut);
        }
        return returnvector;
    }
    
    vector<string> getProductSpecies() {
        vector<string> returnvector;
        for(auto i=getM();i<size();++i){
        string name = (*(rspecies()+i))->getSpecies().getName();
        string namecut = name.substr(0,name.find("-",0));
        returnvector.push_back(namecut);
        }
        return returnvector;
    }
    
    vector<species_copy_t> getReactantCopyNumbers()  {
        vector<species_copy_t> returnvector;
        for(auto i=0U;i<getM();i++)
        {returnvector.push_back((*(rspecies()+i))->getN());}
        return returnvector;
    }
    
    vector<species_copy_t> getProductCopyNumbers()  {
        vector<species_copy_t> returnvector;
        for(auto i=getM();i<size();i++)
        {returnvector.push_back((*(rspecies()+i))->getN());}
        return returnvector;
    }

    
    
    ///Set reaction type
    void setReactionType(ReactionType rxnType) {_reactionType = rxnType;}
    
    ///Get reaction type
    ReactionType getReactionType() {return _reactionType;}
    
    void setGNumber(floatingpoint gnum) {_gnum = gnum;};

	floatingpoint getGNumber() {return _gnum;};
    
    void setHRCDID(string hrcdid) {_hrcdid = hrcdid;};
    
    string getHRCDID() {return _hrcdid;};
    
    ///Set CBound
    void setCBound(CBound* cBound) {_cBound = cBound;}
    ///Get CBound
    CBound* getCBound() {return _cBound;}

    

    // Sets the scaled rate based on volume dependence.
    void recalcRateVolumeFactor() {
        // This can automatically set the "_rate" as scaled value of "rate"

        // Some possibilities of the exponent are implemented specifically to decrease the use of "pow"
        switch(_rateVolumeDepExp) {
            case 0:
                setRateMulFactor(1.0f, VOLUMEFACTOR); break;
            case -1:
	            setRateMulFactor(1.0f / _volumeFrac, VOLUMEFACTOR); break;
            default:
                if(_volumeFrac == 1.0f) setRateMulFactor(1.0f, VOLUMEFACTOR);
                else setRateMulFactor(std::pow(_volumeFrac, _rateVolumeDepExp), VOLUMEFACTOR);
                break;
        }
    }

    auto getRateVolumeDepExp() const { return _rateVolumeDepExp; }

    /// Getter and setter for compartment volume fraction
    auto getVolumeFrac()const { return _volumeFrac; }
    void setVolumeFrac(FP volumeFrac) {
        _volumeFrac = volumeFrac;
        recalcRateVolumeFactor();
    }
    
    /// Sets the RNode pointer associated with this ReactionBase to rhs. Usually is
    /// called only by the Gillespie-like algorithms.
    void setRnode(RNode *rhs) {_rnode=rhs;}
    /// Get the RNode pointer
    RNode* getRNode() {return _rnode;}
    
    /// Returns the rate associated with this ReactionBase.
    auto getRate() const {return _rate;}
    
    /// Returns the bare rate associated with this ReactionBase
    auto getBareRate() const {return _rate_bare;}
    
    ///aravind June 24, 2016
    void setBareRate(FP a) {
        _rate_bare = a;
        _rate = _rate_bare;
        for(unsigned i = 0; i < RateMulFactorType::RATEMULFACTSIZE; i++)
            _rate *= _ratemulfactors[i];
    }
    /// Returns a pointer to the RNode associated with this ReactionBase.
    RNode* getRnode() const {return _rnode;}
    
    /// Returns the number of reactant RSpecies
    unsigned short getM() const {return getMImpl();}
    
    /// (Private) implementation of the getM() method to be elaborated in derived classes.
    virtual unsigned short getMImpl() const = 0;
    
    /// Returns the number of product RSpecies
    unsigned short getN() const {return getNImpl();}
    
    /// (Private) implementation of the getN() method to be elaborated in derived classes.
    virtual unsigned short getNImpl() const = 0;
    
    /// Returns the total number of RSpecies
    unsigned short size() const {return sizeImpl();}
    
    /// (Private) implementation of the size() method to be elaborated in derived classes.
    virtual unsigned short sizeImpl() const = 0;
    
    /// Return the parent Composite object pointer, to which this Reaction belongs to.
    /// If not present, return a nullptr.
    Composite* getParent() {return _parent;}
    
    /// Set the parent Composite object pointer to which this Reaction belongs to.
    void setParent (Composite *other) {_parent=other;}
    
    /// Returns true if this Reaction has a parent Composite object to which it
    /// belongs to.
    bool hasParent() const {return _parent!=nullptr? true : false;}
    
    /// Get the root parent (i.e. follow the pointers of parentage until the root node
    /// in the Composition hieararchy)
    Composite* getRoot();
    
    /// Computes the product of the copy number of all reactant RSpecies.
    /// Can be used to quickly determine whether this ReactionBase should be allowed to
    /// activate - if one of the reactants has a copy number equal to zero, then zero is
    /// returned, indicating that this ReactionBase should not be (yet) activated.
    floatingpoint getProductOfReactants () const {return getProductOfReactantsImpl();}
    
    /// (Private) implementation of the getProductOfReactants() method to be elaborated
    /// in derived classes.
    virtual floatingpoint getProductOfReactantsImpl() const = 0;
    
    /// Computes the product of the copy number of all product RSpecies minus maximum
    /// allowed copy number. Can be used to quickly determine whether this ReactionBase
    /// should be allowed to activate - if one of the products has a copy number equal
    /// to the maximum allowed, then zero is returned, indicating that this ReactionBase
    /// should not be (yet) activated.
    floatingpoint getProductOfProducts () const {return getProductOfProductsImpl();}
    
    /// (Private) implementation of the getProductOfProducts() method to be elaborated
    /// in derived classes.
    virtual floatingpoint getProductOfProductsImpl() const = 0;
    
    /// Return true if the ReactionBase is currently passivated
#if defined TRACK_ZERO_COPY_N || defined  TRACK_UPPER_COPY_N
    bool isPassivated() const {return _passivated;}
#else
    bool isPassivated() const {return false;}
#endif
    
    /// Returns true of this Reaction contains Species *s either as a reactant or a
    /// product
    bool containsSpecies(Species *s) const {return containsSpeciesImpl(s);}
    
    /// (Private) implementation of the containsSpecies() method to be elaborated in
    /// derived classes.
    virtual bool containsSpeciesImpl(Species *s) const = 0;
    
    // Clear all callbacks attached to this reaction.
    void clearSignaling ();
    
    /// Connect the callback, react_callback to a signal corresponding to
    /// ReactionBase *r.
    /// @param CallbackType callback - a function
    /// object to be called (a slot)
    void connect(CallbackType callback);
    
    /// Broadcasts signal indicating that the ReactionBase event has taken place
    /// This method is only called by the code which runs the chemical dynamics (i.e.
    /// Gillespie-like algorithm)
    void emitSignal() {
        for(auto& callback : callbacks_) {
            callback(this);
        }
    }
    
    /// Return a const reference to the set of dependent ReactionBases
    /// @note One can obtain two different lists of affected ReactionBases:
    /// 1) via getAffectedReactionBases(), where the copy numbers do influence the
    /// dependencies, and 2) via dependents(), where dependencies stop being counted
    /// if the copy numbers of reactant species drop to 0.
    const dependentdatatype& dependents() const {return _dependents;}
    
    /// Returns true if two ReactionBase objects are equal.
    /// Two ReactionBase objects are equal if each of their reactants and products
    /// are equal
    friend bool operator==(const ReactionBase& a, const ReactionBase& b)
    {
        if(typeid(a) != typeid(b))
            return false;
        return a.is_equal(b);
        // Invoke virtual is_equal via derived subclass of a
        // (which should be the same as b)
    }
    
    /// (Private) implementation of the operator==(...) method to be elaborated
    /// in derived classes.
    virtual bool is_equal(const ReactionBase& b) const = 0;
    
    /// Return true if two ReactionBase are not equal.
    /// @see operator ==(const ReactionBase& a, const ReactionBase& b) above
    friend bool operator !=(const ReactionBase& a, const ReactionBase& b){
        return !(a==b);
    }
    
    /// Fire the ReactionBase - make a single step, where reactant RSpecies copy
    /// numbers are decreased by one, and the product RSpecies copy numbers are
    /// increased by one.
    /// @note This method does not send a ReactionBase event Signal. The latter is
    /// usually done from within a Gillespie-like algorithm.
    void makeStep() {makeStepImpl();}
    
    /// (Private) implementation of the makeStep() method to be elaborated in
    /// derived classes.
    virtual void makeStepImpl() = 0;
    
    /// Compute the ReactionBase propensity that is needed by a Gillespie like
    /// algorithm, rate*reactant_1.getN()*reactant_2.getN()...
    floatingpoint computePropensity() const {return computePropensityImpl();}
    
    /// (Private) implementation of the computePropensity() method to be elaborated
    /// in derived classes.
    virtual floatingpoint computePropensityImpl() const = 0;
    
    /// Usually is applied to ReactionBase objects with propensity of 0 (e.g. when one
    /// of the copy numbers of reactants has dropped to 0. This method call notifies all
    /// other ReactionBase objects that may affect this ReactionBase to stop tracking
    /// this ReactionBase. Eventually, activateReaction() may be called to restart
    /// tracking, if the propensity stops being 0.
    void passivateReaction() {passivateReactionImpl();}
    
    /// (Private) implementation of the passivateReaction() method to be elaborated in
    /// derived classes.
    virtual void passivateReactionImpl() = 0;
    
    /// Requests that ReactionBase objects that may affect this Reaction to start
    /// tracking it, which can be used to follow ReactionBase objects whose propensities
    /// change upong firing of some ReactionBase. This request is acted upon
    /// unconditionally.
    void activateReactionUnconditional() {activateReactionUnconditionalImpl();}
    
    virtual void activateReactionUnconditionalImpl() = 0;
    
    /// Requests that Reaction objects that may affect this Reaction to start tracking
    /// it, which can be used to follow Reaction objects whose propensities change upon
    /// firing of some Reaction. This request will be ignored if the Reaction's
    /// propensity is still zero.
    void activateReaction();
    
    /// Performs a simple updating of the propensity of this ReactionBase. Does not change
    /// dependents, only updates the RNode if initialized.
    void updatePropensity() {updatePropensityImpl();}
    
    virtual void updatePropensityImpl() = 0;
    
    /// Print the ReactionBases that are affacted by this ReactionBase being fired
    void printDependents() ;
    
    /// Reset the list of ReactionBase objects that are affected when this
    /// ReactionBase is fired
    /// @note This method is "expensive" because it computes from scratch the
    /// dependencies. Importantly, the copy numbers of molecules do not influence the
    /// result of this function. \sa dependents()
    virtual void setDependentReactions() = 0;
    
    /// Request that the ReactionBase *r adds this ReactionBase to its list of
    /// dependents which it affects.
    void registerNewDependent(ReactionBase *r);
    
    /// Request that the ReactionBase *r removes this ReactionBase from its list of
    /// dependents which it affects.
    /// This is usually requested when the ReactionBase propensity drops to zero (i.e.
    /// via passivateReactionBase()).
    void unregisterDependent(ReactionBase *r);
    
    virtual void printToStream(ostream& os) const = 0;
    
    /// Print self into an iostream
    friend ostream& operator<<(ostream& os, const ReactionBase& rr)
    {
    	rr.printToStream(os);
        return os;
    }
    
    ///Whether the dependencies should be updated
    virtual bool updateDependencies() = 0;

    size_t getId() const { return _Id;}
};

} // namespace medyan

#endif

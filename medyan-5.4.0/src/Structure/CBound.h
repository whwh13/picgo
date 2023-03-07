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

#ifndef MEDYAN_CBound_h
#define MEDYAN_CBound_h

#include "common.h"

#include "Species.h"
#include "ReactionBase.h"

namespace medyan {
//FORWARD DECLARATIONS
class Compartment;
class SubSystem;
class CCylinder;

/// Represents a chemical object that is bound to a Filament.
/*!
 *  The CBound class is an abstract representation of a chemically bound object to a 
 *  Filament (Could be a Linker, MotorGhost, BranchingPoint, etc). Each CBound object 
 *  has a pointer to the corresponding [SpeciesBounds] (@ref SpeciesBound) on a Filament.
 *  Different implementations of CBound class will have different functions to bind,
 *  move, etc. See documentation of subclass for more details on function.]
 *
 *  The CBound class also controls the element's unbinding reaction, including its
 *  creation and destruction in the chemical system.
 */
class CBound {
    
protected:
    // The filament types this CBound binds to.
    int filType1_ = 0;
    int filType2_ = 0;
    
    SpeciesBound* _firstSpecies = nullptr; ///< Corresponding first species on Filament
    SpeciesBound* _secondSpecies = nullptr;///< Corresponding second species on Filament
    
    Compartment* _compartment; ///< Compartment this CBound is in
    
    CCylinder* _cc1 = nullptr; ///< Pointer to first CCylinder
    CCylinder* _cc2 = nullptr; ///< Pointer to second CCylinder
    
    short _position1; ///< position of first head
    short _position2 = 0; ///< position of second head
    
    //@{
    ///Reaction rates
    float _onRate = 0.0;
    float _offRate = 0.0;
    //@}
    
    ReactionBase* _offRxn; ///< The off reaction for this bound object
    
public:
    /// Constructor, just sets species
    CBound(int filType1, int filType2, Compartment* c, CCylinder* cc1, CCylinder* cc2,
           short position1, short position2)
    
        : filType1_(filType1), filType2_(filType2), _compartment(c), _cc1(cc1), _cc2(cc2),
          _position1(position1), _position2(position2) {}
    
    /// Virtual destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~CBound() noexcept {
        if(_firstSpecies != nullptr) _firstSpecies->removeCBound();
        if(_secondSpecies != nullptr) _secondSpecies->removeCBound();
    }
    
    /// Set first species
    void setFirstSpecies(SpeciesBound* species) {
        ///remove from old first species
        if(_firstSpecies != nullptr) _firstSpecies->removeCBound();
        
        _firstSpecies = species;
        _firstSpecies->setCBound(this);
    }
    /// Get first species
    SpeciesBound* getFirstSpecies() {return _firstSpecies;}
    
    /// Set second species
    void setSecondSpecies(SpeciesBound* species) {
        ///remove from old second species
        if(_secondSpecies != nullptr) _secondSpecies->removeCBound();
        
        _secondSpecies = species;
        _secondSpecies->setCBound(this);
    }
    /// Get second species
    SpeciesBound* getSecondSpecies() {return _secondSpecies;}
    
    //@{
    /// Position getters
    short getFirstPosition() {return _position1;}
    short getSecondPosition() {return _position2;}
    //@}
    
    /// Set first CCylinder
    void setFirstCCylinder(CCylinder* cc) {_cc1 = cc;}
    /// Get first CCylinder
    CCylinder* getFirstCCylinder() {return _cc1;}
    
    /// Set second CCylinder
    void setSecondCCylinder(CCylinder* cc) {_cc2 = cc;}
    /// Get first CCylinder
    CCylinder* getSecondCCylinder() {return _cc2;}
    
    /// Get compartment that this CBound is in
    Compartment* getCompartment() {return _compartment;}
    
    //@{
    /// On rate management
    void setOnRate(float rate) {_onRate = rate;}
    floatingpoint getOnRate(){return _onRate;}
    //@}
    
    //@{
    /// Off rate management
    void setOffRate(float rate) {_offRate = rate;}
    floatingpoint getOffRate(){return _offRate;}
    //@}
    
    /// Set all rates at once
    void setRates(float onRate, float offRate) {
        _onRate = onRate;
        _offRate = offRate;
    }
    
    //@{
    /// Off reaction management
    virtual void createOffReaction(ReactionBase* onRxn, SubSystem* ps) = 0;
    
    void setOffReaction(ReactionBase* offRxn) {
        _offRxn = offRxn;
        _offRxn->setCBound(this);
    }
    ReactionBase* getOffReaction() {return _offRxn;}
    //@}
};

} // namespace medyan

#endif


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

#ifndef MEDYAN_CMotorGhost_h
#define MEDYAN_CMotorGhost_h

#include "common.h"

#include "CBound.h"

#define SPECIESM_BINDING_INDEX 1
#define SPECIESM_UNBINDING_INDEX 2
#define SPECIESM_DIFFUSING_INDEX_OFFRXN 2

namespace medyan {
//FORWARD DECLARATIONS
class MotorGhost;
class SubSystem;
class CCylinder;
class Compartment;

/// A class to represent the chemical component of a MotorGhost.
/*!
 *  The CMotorGhost class contains chemical info of the parent MotorGhost.
 *
 *  Extending CBound, this class tracks its corresponding Species and unbinding Reaction.
 */
class CMotorGhost : public CBound {
    
private:
    MotorGhost* _pMotorGhost; ///< Pointer to parent
    
public:
    ///Constructor
    ///@param pos1 - monomer index on first cylinder
    ///@param pos2 - monomer index on second cylinder
    CMotorGhost(
        int motorSpeciesIndex1, int motorSpeciesIndex2,
        Compartment* c,
        CCylinder* cc1, CCylinder* cc2, int position1, int position2);
    
    ///Destructor, removes off reaction from system
    ~CMotorGhost();
    
    /// Copy constructor, standard
    CMotorGhost(const CMotorGhost& rhs, Compartment* c)
    
        : CBound(rhs.filType1_, rhs.filType2_, c, rhs._cc1, rhs._cc2, rhs._position1, rhs._position2),
          _pMotorGhost(rhs._pMotorGhost) {
        
        //set species
        setFirstSpecies(rhs._firstSpecies);
        setSecondSpecies(rhs._secondSpecies);

        //set reaction
        setOffReaction(rhs._offRxn);
        //set rates
        setOnRate(rhs._onRate);
        setOffRate(rhs._offRate);
              
    }
    
    /// Assignment is not allowed
    CMotorGhost& operator=(CMotorGhost &rhs) = delete;
    
    /// Clone, calls copy constructor
    virtual CMotorGhost* clone(Compartment* c) {
        CMotorGhost* cm = new CMotorGhost(*this, c);
        _offRxn = nullptr; return cm;
    }
    
    /// Set parent
    void setMotorGhost(MotorGhost* MotorGhost) {_pMotorGhost = MotorGhost;}
    /// Get parent
    MotorGhost* getMotorGhost() {return _pMotorGhost;}
        
    
    /// Create the off reaction for this MotorGhost
    virtual void createOffReaction(ReactionBase* onRxn, SubSystem* ps);
    
    /// Move the motor head to a new position chemically
    void moveMotorHead(CCylinder* cc,
                       short oldPosition,
                       short newPosition,
                       short speciesMotorIndex,
                       short boundType,
                       SubSystem* ps);
    
    /// Move the motor head to a new CCylinder chemically
    void moveMotorHead(CCylinder* oldCC,
                       CCylinder* newCC,
                       short oldPosition,
                       short newPosition,
                       short speciesMotorIndex,
                       short boundType,
                       SubSystem* ps);
    
    void printReaction();

    Species* getDiffusingSpecies(){
        RSpecies** rs = _offRxn->rspecies();
        Species* sfb = &(rs[SPECIESM_DIFFUSING_INDEX_OFFRXN]->getSpecies());
        return sfb;
    }
};

} // namespace medyan

#endif

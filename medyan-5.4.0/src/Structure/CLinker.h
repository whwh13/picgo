
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

#ifndef MEDYAN_CLinker_h
#define MEDYAN_CLinker_h

#include "common.h"

#include "CBound.h"

#define SPECIESL_BINDING_INDEX 1
#define SPECIESL_DIFFUSING_INDEX_OFFRXN 2

namespace medyan {
//FORWARD DECLARATIONS
class Linker;
class SubSystem;
class CCylinder;
class Compartment;

/// To represent the chemical component of a Linker.
/*! 
 *  The CLinker class contains chemical info of the parent Linker.
 *
 *  Extending CBound, this class tracks its corresponding Species and unbinding Reaction.
 */
class CLinker : public CBound {
    
private:
    Linker* _pLinker; ///< Pointer to parent

public:
    ///Constructor
    ///@param pos1 - monomer index on first cylinder
    ///@param pos2 - monomer index on second cylinder
    CLinker(
        int linkerSpeciesIndex1,
        int linkerSpeciesIndex2,
        Compartment* c,
        CCylinder* cc1, CCylinder* cc2, int position1, int position2);
    
    ///Destructor, removes off reaction from system
    ~CLinker();
    
    /// Copy constructor, standard
    CLinker(const CLinker& rhs, Compartment* c)
        : CBound(rhs.filType1_, rhs.filType2_, c, rhs._cc1, rhs._cc2, rhs._position1, rhs._position2),
          _pLinker(rhs._pLinker){
        
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
    CLinker& operator=(CLinker &rhs) = delete;
    
    /// Clone, calls copy constructor
    virtual CLinker* clone(Compartment* c) {
        CLinker* cl = new CLinker(*this, c);
        _offRxn = nullptr; return cl;
    }
    
    /// Set parent
    void setLinker(Linker* linker) {_pLinker = linker;}
    /// Get parent 
    Linker* getLinker() {return _pLinker;}
    
    /// Create the off reaction for this Linker
    virtual void createOffReaction(ReactionBase* onRxn, SubSystem* ps);

    Species* getDiffusingSpecies(){
        RSpecies** rs = _offRxn->rspecies();
        Species* sfb = &(rs[SPECIESL_DIFFUSING_INDEX_OFFRXN]->getSpecies());
        return sfb;
    }
    
};

} // namespace medyan

#endif

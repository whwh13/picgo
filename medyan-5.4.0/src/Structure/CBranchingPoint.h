
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

#ifndef MEDYAN_CBranchingPoint_h
#define MEDYAN_CBranchingPoint_h

#include "common.h"

#include "CBound.h"
#include "Compartment.h"

#define SPECIESB_BINDING_INDEX 0
#define SPECIESB_DIFFUSING_INDEX_OFFRXN 2
#define SPECIESA_DIFFUSING_INDEX_ONRXN 1

namespace medyan {
//FORWARD DECLARATIONS
class BranchingPoint;
class SubSystem;

/// A class to represent the chemical component of a BranchingPoint.
/*!
 *  The CBranchingPoint class contains chemical info of the parent BranchingPoint.
 *
 *  Extending CBound, this class tracks its corresponding Species and unbinding Reaction.
 */
class CBranchingPoint : public CBound {
    
private:
    BranchingPoint* _pBranchingPoint; ///< Pointer to parent
    
    short _branchType; ///< Branching point type

    string diffusingactinspeciesname="";
public:
    /// Default constructor and destructor
    /// @param pos - monomer index on first cylinder
    CBranchingPoint(short branchType, Compartment* c,
                    CCylinder* cc1, CCylinder* cc2, int position);
    ///Destructor, removes off reaction from system
    ~CBranchingPoint();
    
    /// Copy constructor, standard
    CBranchingPoint(const CBranchingPoint& rhs, Compartment* c)
        : CBound(rhs.filType1_, rhs.filType2_, c, rhs._cc1, rhs._cc2, rhs._position1, rhs._position2),
          _pBranchingPoint(rhs._pBranchingPoint) {
        
          //set species
          setFirstSpecies(rhs._firstSpecies);
            
          //set reaction
          setOffReaction(rhs._offRxn);
            
          //set rates
          setOnRate(rhs._onRate);
          setOffRate(rhs._offRate);
          setdiffusingactinspeciesname(rhs.getdiffusingactinspeciesname());

    }
    
    /// Assignment is not allowed
    CBranchingPoint& operator=(CBranchingPoint &rhs) = delete;
    
    /// Clone, calls copy constructor
    virtual CBranchingPoint* clone(Compartment* c) {
        
        CBranchingPoint* cb = new CBranchingPoint(*this, c);
        _offRxn = nullptr; return cb;
    }
    
    /// Set parent
    void setBranchingPoint(BranchingPoint* BranchingPoint) {
        _pBranchingPoint = BranchingPoint;
    }
    /// Get parent
    BranchingPoint* getBranchingPoint() {return _pBranchingPoint;}
    
    virtual void createOffReaction(ReactionBase* onRxn, SubSystem* ps);

    Species* getDiffusingBranchSpecies(){
        RSpecies** rs = _offRxn->rspecies();
        Species* sfb = &(rs[SPECIESB_DIFFUSING_INDEX_OFFRXN]->getSpecies());
        return sfb;
    }

    string getdiffusingactinspeciesname() const { return diffusingactinspeciesname;}

    void setdiffusingactinspeciesname(string _diffusingactinspeciesname) {
        diffusingactinspeciesname = _diffusingactinspeciesname;
    }
};

} // namespace medyan

#endif

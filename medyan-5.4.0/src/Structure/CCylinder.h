
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

#ifndef MEDYAN_CCylinder_h
#define MEDYAN_CCylinder_h

#include <map>
#include <unordered_set>

#include "common.h"

#include "CMonomer.h"

#include "SysParams.h"

namespace medyan {
//FORWARD DECLARATIONS
class Compartment;
class Cylinder;
class ChemSim;
class ChemManager;

/// Holds all [CMonomers](@ref CMonomer) and [Reactions](@ref Reaction) associated with it.
/*! 
 *  The CCylinder class has lists of the [CMonomers](@ref CMonomer) that it contains.
 *
 *  A CCylinder also stores the [Reactions] (@ref Reaction) associated with it 
 *  internally, as well as cross-cylinder reactions (still "owned" by this CCylinder,
 *  but has reactants and/or products that contain a species in another CCylinder). 
 *  This information is stored in the reaction map. This map can be used to add and 
 *  delete reactions when a CCylinder changes compartment or is removed from the system.
 *
 *  Lastly, the CCylinder stores a list of other CCylinders that it has cross-cylinder
 *  reactions with, but does not have ownership of these reactions. This list can be 
 *  used to change or delete reactions when needed, but will not be tied to a 
 *  Compartment transfer of this CCylinder.
 */
class CCylinder {
    
friend class CController;
friend class ChemManager;
    
private:
#ifdef OPTIMOUT
    chrono::high_resolution_clock::time_point mins, mine, minsi, minei;
#endif
    #ifdef DEBUGCONSTANTSEED
	using internalreactiondatatype = unordered_set<ReactionBase*, HashbyId<ReactionBase*>,
			customEqualId<ReactionBase*>>;
//	using internalreactiondatatype = set<ReactionBase*, customcompareId<ReactionBase*>>;
	using reactingcylindersdatatype = set<CCylinder*, customcompareId<CCylinder*>>;
	using crosscylinderreactionsdatatype = map<CCylinder*,unordered_set<ReactionBase*, HashbyId<ReactionBase*>,
			customEqualId<ReactionBase*>>>;
	#else
    using internalreactiondatatype = unordered_set<ReactionBase*>;
    using reactingcylindersdatatype = unordered_set<CCylinder*>;//This can be rewritten
    // as a set of CCylinder IDs.
    using crosscylinderreactionsdatatype = map<CCylinder*,unordered_set<ReactionBase*>>;
	#endif

    vector<unique_ptr<CMonomer>> _monomers; ///< List of monomers
    
    ///REACTION CONTAINERS
    internalreactiondatatype _internalReactions;///< Set of internal reactions associated

	reactingcylindersdatatype _reactingCylinders;   ///< Set of ccylinders that this ccylinder has
                                                    ///< reactions with, but not ownership
    crosscylinderreactionsdatatype _crossCylinderReactions;
    ///< Map of cross-cylinder reactions owned
    
    Compartment* _compartment; ///< Compartment this ccylinder is in
    Cylinder* _pCylinder;      ///< Parent cylinder
    
    short _size = 0; ///< Maximum length
    
    inline static medyan::ChemSim* _chemSim = nullptr;   ///< A pointer to the ChemSim, initialized by CController

    bool initialized = false;//CCylinder is completely self-consistent only after
    // initializeCCylinder function is called from ChemManager. This variable will help
    // us track whether a CCylinder has been initialized properly.
    
public:
    /// Default constructor, sets compartment and cylinder
    CCylinder(Compartment* C, Cylinder* c);
    
    /// Copy constructor
    /// @note This constructor will create a new CCylinder with different Species and
    /// [Reactions](@ref Reaction) within the Compartment that is chosen as a parameter
    /// to the constructor. The copied and original CCylinder will not share reactions
    /// or species, but will be copied into a new Compartment.
    CCylinder(const CCylinder& rhs, Compartment* c);
    
    /// Assignment is not allowed
    CCylinder& operator=(CCylinder &rhs) = delete;
    
    /// Default destructor, explicitly removes CMonomers(including their species, rxns)
    /// Removes all reactions associated with this CCylinder, including ones owned by
    /// this as well as other CCylinders.
    ~CCylinder();
    
    /// Clone, calls copy constructor
    CCylinder* clone(Compartment* c) {
        return new CCylinder(*this, c);
    }
    
    /// Get compartment
    Compartment* getCompartment() {return _compartment;}
    
    /// Set parent
    void setCylinder(Cylinder* c) {_pCylinder = c;}
    /// Get parent
    Cylinder* getCylinder() {return _pCylinder;}
    
    /// Add a monomer
    void addCMonomer(CMonomer* monomer) {
        _monomers.emplace_back(unique_ptr<CMonomer>(monomer));
    }
    /// Get monomer at an index
    /// @note no check on index
    CMonomer* getCMonomer(int index) {return _monomers[index].get();}
    
    ///Get list of reactions associated
    const internalreactiondatatype& getInternalReactions() {return _internalReactions;}
    
    ///Get list of reacting cylinders associated
    const reactingcylindersdatatype& getReactingCylinders() {return _reactingCylinders;}
    
    ///Get map of reactions associated
    crosscylinderreactionsdatatype& getCrossCylinderReactions() {
        return _crossCylinderReactions;
    }
    
    //@{
    /// Reaction management function
    void addInternalReaction(ReactionBase* r);
    void removeInternalReaction(ReactionBase* r);
    void removeAllInternalReactions();
    
    void addCrossCylinderReaction(CCylinder* other, ReactionBase* r);
    void removeCrossCylinderReaction(CCylinder* other, ReactionBase* r);
    void removeCrossCylinderReactions(CCylinder* other);
    void removeAllCrossCylinderReactions();
    
    void addReactingCylinder(CCylinder* other);
    void removeReactingCylinder(CCylinder* other);
    void removeAllReactingCylinders();
    
    void passivatefilreactions();
    void activatefilreactions();
    
    void passivatefilcrossreactions();
    void activatefilcrossreactions();
    
    //@}
    
    /// Get all reactions that this CCylinder has ownership of
    vector<ReactionBase*> getAllReactions();
    
    /// Print
    void printCCylinder();
    
    //Check if all chemical species held in this CCylinder
    //are self-consistent. For debugging purposes only.
    bool isConsistent();
    
    /// Get size in number of monomers
    short getSize() {return _size;}
    
    /// Get the type
    short getType();

    int getId();

    bool isinitialized(){return initialized;}
    
};

} // namespace medyan

#endif

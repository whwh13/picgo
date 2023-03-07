
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_BindingManager_h
#define MEDYAN_BindingManager_h

#include <unordered_map>
#include <unordered_set>
#include <random>

#include "common.h"

#include "NeighborListImpl.h"
#include "HybridBindingSearchManager.h"
#include "ReactionBase.h"

#include "SysParams.h"
#include "Rand.h"

#define B_RXN_INDEX 2
#define ML_RXN_INDEX 0

namespace medyan {
//FORWARD DECLARATIONS
class SubSystem;
class ReactionBase;
class CCylinder;
class Compartment;
class Cylinder;


///Enumeration for nucleation zone type. Used by BranchingManager.
enum NucleationZoneType {
    ALL, BOUNDARY, TOPBOUNDARY, SIDEBOUNDARY, RIGHTBOUNDARY, MEMBRANE
};

/// To store and manage binding reactions.

/*!
 *  FilamentBindingManager is used to store a binding Reaction on [Filaments](@ref Filament)
 *  in Compartment. Contains the binding reaction, possible binding sites, and integers
 *  representing the binding species involved in the reaction. Classes that extend this will
 *  implement their own data structures for holding possible reaction sites, etc.
 *
 *  The main function of this class is to call updatePossibleBindings(), which will
 *  update the possible binding sites if the binding reaction is called in this compartment.
 *  Also contains functions to replace ccylinders in the structure, as well as remove a ccylinder
 *  from the structure when it is removed from the subsystem.
 *
 *  When the binding reaction is fired, the manager will choose a random binding based on its current
 *  data structure state, and perform whatever callback is associated with that binding reaction.
 */
class FilamentBindingManager {

friend class ChemManager;

protected:
    ReactionBase* _bindingReaction; ///< The binding reaction for this compartment

    Compartment* _compartment; ///< Compartment this is in

    short _boundInt; ///< Integer index in CMonomer of bound chemical value.
                     ///< @note - THIS ALSO REPRESENTS THE SPECIES TYPE THAT IS MANAGED.

    vector<short> _filamentIDvec; ///< The filament type to operate on

    Species* _bindingSpecies; ///< The binding species that this manager tracks.
                              ///< Resposible for all copy number changes

    short _nlIndex = 0; ///<Index of this manager (for access of neighbor lists)
    short _mIndex = 0;  ///<Index of this manager (for access in other compartments)


    static SubSystem *_subSystem; ///< Ptr to the SubSystem



#if defined(NLSTENCILLIST) || defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
    vector<floatingpoint> bindingSites;
#endif

    ///helper function to update copy number and reactions
    void updateBindingReaction( int oldN, int newN) {

        int diff = newN - oldN;
        //update copy number
        if(diff > 0) {
            while (diff != 0) {
                _bindingSpecies->up();
                diff--;
            }
        }
        else if(diff < 0) {
            while (diff != 0) {
                _bindingSpecies->down();
                diff++;
            }
        }
        else {} //do nothing
        //check if matching
#ifdef NLORIGINAL

        assert((_bindingSpecies->getN() == numBindingSites())
               && "Species representing binding sites does not match \
                   number of binding sites held by the manager.");
#endif
#if defined(NLSTENCILLIST) || defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
                assert((_bindingSpecies->getN() == numBindingSitesstencil())
               && "Species representing binding sites does not match \
                   number of binding sites held by the manager.");
#endif
        _bindingReaction->updatePropensity();
    }
public:
    FilamentBindingManager(ReactionBase* reaction,
                           Compartment* compartment,
                           short boundInt,
                           vector<short> filamentIDvec)

    : _bindingReaction(reaction), _compartment(compartment),
      _boundInt(boundInt),
      _filamentIDvec(filamentIDvec){

    }
    virtual ~FilamentBindingManager() = default;

    //@{
    ///add possible binding reactions that could occur
#ifdef NLORIGINAL
    virtual void addPossibleBindings(CCylinder* cc, short bindingSite) = 0;

    virtual void addPossibleBindings(CCylinder* cc) = 0;
    //@}

    //@{
    /// Remove all bindings including this cylinder
    virtual void removePossibleBindings(CCylinder* cc, short bindingSite) = 0;
    virtual void removePossibleBindings(CCylinder* cc) = 0;
    //@}
    ///update all possible binding reactions that could occur
    virtual void updateAllPossibleBindings() = 0;
     /// Get current number of binding sites
    virtual int numBindingSites() = 0;

    /// ARAVIND ADDED FEB 17 2016. append possible bindings.
    virtual void appendpossibleBindings(short boundInt, short bountCCylinder* ccyl1,
    		CCylinder* ccyl2,
    		short site1, short site2)=0;

    virtual void clearpossibleBindings()=0;

    virtual void printbindingsites() = 0;
#endif
    ///Get the bound species integer index
    short getBoundInt() {return _boundInt;}

    ///Set the index of this manager, for access to NeighborList
    void setNLIndex(int index) {_nlIndex = index;}

    ///Set the index of this manager, for access to other managers
    void setMIndex(int index) {_mIndex = index;}
    ///Check consistency and correctness of binding sites. Used for debugging.
    virtual bool isConsistent() = 0;

    ///get the filament that the species binds to aravind June 30, 2016.
    auto& getfilamentTypes() const { return _filamentIDvec; }

    ///aravind, June 30,2016.
    vector<string> getrxnspecies(){return _bindingReaction->getreactantspecies();}

#if defined(NLSTENCILLIST) || defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
    virtual void addPossibleBindingsstencil(CCylinder* cc) = 0;
    virtual void addPossibleBindingsstencil(CCylinder* cc, short bindingSite) = 0;
    ///update all possible binding reactions that could occur using stencil NL
    virtual void updateAllPossibleBindingsstencil() = 0;
    virtual void appendpossibleBindingsstencil(short boundInt, CCylinder* ccyl1, CCylinder* ccyl2, short site1,
                                               short site2)=0;
    virtual void clearpossibleBindingsstencil()=0;
    virtual int numBindingSitesstencil() = 0;
    virtual void removePossibleBindingsstencil(CCylinder* cc) = 0;
    virtual void removePossibleBindingsstencil(CCylinder* cc, short bindingSite) = 0;
    virtual void crosscheck(){};
    virtual void replacepairPossibleBindingsstencil
            (unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>> _pbs) =0;
    virtual void printbindingsitesstencil() = 0;
    virtual void setHNLID(short id, short idvec[2]) = 0;
#endif
    void updateBindingReaction(int newN) {
        int oldN = _bindingSpecies->getN();
        int diff = newN - oldN;
        //update copy number
        if(diff > 0) {
            while (diff != 0) {
                _bindingSpecies->up();
                diff--;
            }
        }
        else if(diff < 0) {
            while (diff != 0) {
                _bindingSpecies->down();
                diff++;
            }
        }
        else {} //do nothing
        //check if matching
        assert((_bindingSpecies->getN() == newN)
               && "Species representing binding sites does not match \
                   number of binding sites held by the manager.");
        _bindingReaction->updatePropensity();
    }

};


/// Manager for Filament and BranchingPoint creation
class BranchingManager : public FilamentBindingManager {

friend class ChemManager;

private:
    ///Nucleation zone type, to define where nucleation should occur
    NucleationZoneType _nucleationZone;

    ///If using a nucleation zone, nucleating distance from the boundary
    floatingpoint _nucleationDistance = 0.0;

    ///possible bindings at current state
//    #elif defined(NLORIGINAL)
    unordered_set<tuple<CCylinder*, short>> _possibleBindings;
//    #elif defined(NLSTENCILLIST)
    ///possible bindings at current state in Stencil list.
    unordered_set<tuple<CCylinder*, short>> _possibleBindingsstencil;
    vector<tuple<tuple<CCylinder*, short>, tuple<CCylinder*, short>>> _branchrestarttuple; //Used only during restart conditions.

public:
    BranchingManager(ReactionBase* reaction,
                     Compartment* compartment,
                     short boundInt,
                     vector<short> _filamentIDvec,
                     NucleationZoneType zone = NucleationZoneType::ALL,
                     floatingpoint nucleationDistance = numeric_limits<floatingpoint>::infinity());
    ~BranchingManager() {}

    //@{
#ifdef NLORIGINAL
    ///add possible binding reactions that could occur
    virtual void addPossibleBindings(CCylinder* cc, short bindingSite);

    virtual void addPossibleBindings(CCylinder* cc);
    //@}

    //@{
    /// Remove all bindings including this cylinder
    virtual void removePossibleBindings(CCylinder* cc, short bindingSite);
    virtual void removePossibleBindings(CCylinder* cc);
    //@}

    ///update all possible binding reactions that could occur
    virtual void updateAllPossibleBindings();

    virtual int numBindingSites() {

        return _possibleBindings.size();
    }

    /// Choose a random binding site based on current state
    tuple<CCylinder*, short> chooseBindingSite() {
        assert((_possibleBindings.size() != 0)
               && "Major bug: Branching manager should not have zero binding \
                  sites when called to choose a binding site.");
        int randomIndex = Rand::randInteger(0, _possibleBindings.size() - 1);
        auto it = _possibleBindings.begin();

        advance(it, randomIndex);

        return *it;
    }
    /// ARAVIND ADDED FEB 17 2016. append possible bindings to be used for restart
    virtual void appendpossibleBindings(short boundInt, CCylinder* ccyl1, CCylinder* ccyl2, short site1,
                                        short site2);

    virtual void clearpossibleBindings() {
        floatingpoint oldN=numBindingSites();
        _possibleBindings.clear();
        updateBindingReaction(oldN,0);
    }
    virtual unordered_set<tuple<CCylinder*, short>> getpossibleBindings(){
        return _possibleBindings;
    }

    virtual void printbindingsites();
#endif
    virtual bool isConsistent();
    vector<tuple<tuple<CCylinder*, short>, tuple<CCylinder*, short>>> getbtuple() {
        return _branchrestarttuple;
    }
#if defined(NLSTENCILLIST) || defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
    virtual void addPossibleBindingsstencil(CCylinder* cc);
    virtual void addPossibleBindingsstencil(CCylinder* cc, short bindingSite);
    ///update all possible binding reactions that could occur using stencil NL
    virtual void updateAllPossibleBindingsstencil();
    virtual void appendpossibleBindingsstencil(short boundInt, CCylinder* ccyl1, CCylinder* ccyl2, short site1,
                                               short site2);
    virtual void clearpossibleBindingsstencil() {
        floatingpoint oldN=numBindingSitesstencil();
        _possibleBindingsstencil.clear();
        updateBindingReaction(oldN,0);
    }
    virtual int numBindingSitesstencil() {

        return _possibleBindingsstencil.size();
    }
    virtual void removePossibleBindingsstencil(CCylinder* cc, short bindingSite);
    virtual unordered_set<tuple<CCylinder*, short>> getpossibleBindingsstencil(){
        return _possibleBindingsstencil;
    }
    tuple<CCylinder*, short> chooseBindingSitesstencil() {

        assert((_possibleBindingsstencil.size() != 0)
               && "Major bug: Branching manager should not have zero binding \
                  sites when called to choose a binding site.");

        int randomIndex = Rand::randInteger(0, _possibleBindingsstencil.size() - 1);
        auto it = _possibleBindingsstencil.begin();

        advance(it, randomIndex);

        return *it;
    }
    virtual void removePossibleBindingsstencil(CCylinder* cc);
    virtual void crosscheck();
    virtual void replacepairPossibleBindingsstencil
            (unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>>){};
    virtual void setHNLID(short id, short idvec[2]){};

    virtual void printbindingsitesstencil();
#endif

#ifdef CUDAACCL_NL
    floatingpoint *gpu_distance;
    int *gpu_zone;
    int *gpu_numpairs = NULL;
    void assigncudavars();
    void freecudavars();
    floatingpoint* getdistancesCUDA(){return gpu_distance;}
    int* getzoneCUDA();
    int* getnumpairsCUDA();
    int* getpossiblebindingssizeCUDA(){ return gpu_numpairs;}
#endif
};

/// Manager for Linker binding.
/*!
 *  LinkerBindingManager controls linker binding in a compartment.
 *  Manages a multimap of possible binding sites.
 */
class LinkerBindingManager : public FilamentBindingManager {

friend class ChemManager;

private:
    float _rMin; ///< Minimum reaction range
    float _rMax; ///< Maximum reaction range
	float _rMinsq = _rMin * _rMin; ///< Minimum reaction range squared
	float _rMaxsq = _rMax * _rMax; ///< Maximum reaction range squared
    vector<floatingpoint> bindingsites1;
	vector<floatingpoint> bindingsites2;
    short HNLID = 1000;
    short _Hbsmidvec[2]= {1000,1000};
    int dBInt = 1;
    int dBI = SysParams::Chemistry().linkerbindingskip-1;

    //possible bindings at current state. updated according to neighbor list
    unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>> _possibleBindings;

    unordered_map<tuple<CCylinder*, short>, vector<tuple<CCylinder*, short>>> _reversePossibleBindings;


        //possible bindings at current state. updated according to neighbor list stencil
    unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>>
            _possibleBindingsstencil;
    //static neighbor list
    static vector<CylinderCylinderNL*> _neighborLists;

    // The linker species index on each filament type.
    int linkerSpeciesIndices_[2] {-1, -1};

public:
    LinkerBindingManager(
        ReactionBase* reaction,
        Compartment* compartment,
        short linkerType,
        vector<short> _filamentIDvec,
        int linkerSpeciesIndex1,
        int linkerSpeciesIndex2,
        float rMax, float rMin);


    //@{
#ifdef NLORIGINAL
    ///add possible binding reactions that could occur
    virtual void addPossibleBindings(CCylinder* cc, short bindingSite);

    virtual void addPossibleBindings(CCylinder* cc);
    //@}

    //@{
    /// Remove all bindings including this cylinder
    virtual void removePossibleBindings(CCylinder* cc, short bindingSite);
    virtual void removePossibleBindings(CCylinder* cc);
    //@}

    ///update all possible binding reactions that could occur
    virtual void updateAllPossibleBindings();

    virtual int numBindingSites() {

        return _possibleBindings.size();
    }

    /// ARAVIND ADDED FEB 17 2016. append possible bindings.
    virtual void appendpossibleBindings(short boundInt, CCylinder* ccyl1, CCylinder* ccyl2, short site1,
                                        short site2);

    //get possible bindings.
    virtual unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>> getpossibleBindings(){
        return _possibleBindings;
    }

    /// ARAVIND ADDED FEB 18 2016. clear possible bindings.
    virtual void clearpossibleBindings() {
        floatingpoint oldN=numBindingSites();
        _possibleBindings.clear();
        updateBindingReaction(oldN,0);
    }

        /// Choose random binding sites based on current state
    vector<tuple<CCylinder*, short>> chooseBindingSites();

    virtual void printbindingsites();
#endif
    //@{
    /// Getters for distances
    float getRMin() {return _rMin;}
    float getRMax() {return _rMax;}
    //@}

    // Get linker species indices on each filament type.
    auto& getLinkerSpeciesIndices() const { return linkerSpeciesIndices_; }

    virtual bool isConsistent();

#if defined(NLSTENCILLIST) || defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
    virtual void addPossibleBindingsstencil(CCylinder* cc);
    virtual void addPossibleBindingsstencil(CCylinder* cc, short bindingSite);
    ///update all possible binding reactions that could occur using stencil NL
    virtual void updateAllPossibleBindingsstencil();
    virtual void appendpossibleBindingsstencil(short boundInt, CCylinder* ccyl1, CCylinder* ccyl2, short site1,
                                               short site2);
    virtual void clearpossibleBindingsstencil();
    virtual int numBindingSitesstencil();
    virtual unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>>
    getpossibleBindingsstencil(){
        return _possibleBindingsstencil;
    }
    vector<tuple<CCylinder*, short>> chooseBindingSitesstencil();
    virtual void removePossibleBindingsstencil(CCylinder* cc, short bindingSite);
    virtual void removePossibleBindingsstencil(CCylinder* cc);
    virtual void crosscheck();
    virtual void replacepairPossibleBindingsstencil
            (unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>> _pbs){
//        _possibleBindingsstencil.clear();
        _possibleBindingsstencil = _pbs;
        int oldN = _bindingSpecies->getN();
        int newN = numBindingSitesstencil();
        updateBindingReaction(oldN, newN);
    }
    virtual void setHNLID(short id, short idvec[2]){
        HNLID = id;
        _Hbsmidvec[0] = idvec[0];
        _Hbsmidvec[1] = idvec[1];
    };

    virtual void printbindingsitesstencil();
#endif

#ifdef CUDAACCL_NL
    floatingpoint *gpu_rminmax = NULL;
    int *gpu_numpairs = NULL;
    void assigncudavars();
    void freecudavars();
    floatingpoint* getdistancesCUDA() { return gpu_rminmax;}
    int* getpossiblebindingssizeCUDA(){ return gpu_numpairs;}
    int getNLsize(){
        return _neighborLists[_nlIndex]->getNLsize();
    }
    int* getNLCUDA(){
        return _neighborLists[_nlIndex]->getNLCUDA();
    }
    int* getNLsizeCUDA(){
        return  _neighborLists[_nlIndex]->getNLsizeCUDA();
    }
#endif

};

/// Manager for MotorGhost binding
/*!
 *  MotorBindingManager controls motor binding in a compartment.
 *  Manages a multimap of possible binding sites.
 */
class MotorBindingManager : public FilamentBindingManager {

friend class ChemManager;

private:
        //DEPRECATED AS OF 9/22/16
//    vector<int> _unboundIDs = {};
    ///< A vector of unbound motor ID's that are contained in this compartment. This is used
    ///< for tracking binding/unbinding and movement of specific motors.

    float _rMin; ///< Minimum reaction range
    float _rMax; ///< Maximum reaction range
	float _rMinsq = _rMin * _rMin; ///< Minimum reaction range squared
	float _rMaxsq = _rMax * _rMax; ///< Maximum reaction range squared
    vector<floatingpoint> bindingsites1;
    vector<floatingpoint> bindingsites2;
    short HNLID = 1000;//Tells which HneighborList to use
    short _Hbsmidvec[2] = {1000, 1000};//Tells position in HybridBindingSearchManager
    //possible bindings at current state. updated according to neighbor list
    unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>>
    _possibleBindings;

    unordered_map<tuple<CCylinder*, short>, vector<tuple<CCylinder*, short>>> _reversePossibleBindings;


        //possible bindings at current state. updated according to neighbor list stencil
    unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>>
            _possibleBindingsstencil;


    //static neighbor list
    static vector<CylinderCylinderNL*> _neighborLists;

    // The motor species index on each filament type.
    int motorSpeciesIndices_[2] {-1, -1};

public:

    MotorBindingManager(
        ReactionBase* reaction,
        Compartment* compartment,
        short linkerType,
        vector<short> _filamentIDvec,
        int motorSpeciesIndex1,
        int motorSpeciesIndex2,
        float rMax, float rMin);


    //@{
#ifdef NLORIGINAL
    ///add possible binding reactions that could occur
    virtual void addPossibleBindings(CCylinder* cc, short bindingSite);
    virtual void addPossibleBindings(CCylinder* cc);
    //@}

    //@{
    /// Remove all bindings including this cylinder
    virtual void removePossibleBindings(CCylinder* cc, short bindingSite);
    virtual void removePossibleBindings(CCylinder* cc);
    //@}

    ///update all possible binding reactions that could occur
    virtual void updateAllPossibleBindings();

    virtual int numBindingSites() {

        return _possibleBindings.size();
    }

    /// ARAVIND ADDED FEB 22 2016. append possible bindings.
    virtual void appendpossibleBindings(short boundInt, CCylinder* ccyl1, CCylinder* ccyl2, short site1,
    		short site2);
    //get possible bindings.
    virtual unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>> getpossibleBindings(){
        return _possibleBindings;
    }

    /// ARAVIND ADDED FEB 18 2016. clear possible bindings.
    virtual void clearpossibleBindings() {
        floatingpoint oldN=numBindingSites();
        _possibleBindings.clear();
        updateBindingReaction(oldN,0);
    }

    /// Choose random binding sites based on current state
    vector<tuple<CCylinder*, short>> chooseBindingSites();

    virtual void printbindingsites();
#endif
    //@{
    /// Getters for distances
    float getRMin() {return _rMin;}
    float getRMax() {return _rMax;}
    //@}

    // Get motor species indices on each filament type.
    auto& getMotorSpeciesIndices() const { return motorSpeciesIndices_; }

    virtual bool isConsistent();
#if defined(NLSTENCILLIST) || defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
    virtual void addPossibleBindingsstencil(CCylinder* cc);
    virtual void addPossibleBindingsstencil(CCylinder* cc, short bindingSite);
    ///update all possible binding reactions that could occur using stencil NL
    virtual void updateAllPossibleBindingsstencil();
    virtual void appendpossibleBindingsstencil(short boundInt, CCylinder* ccyl1,
                                               CCylinder* ccyl2, short site1, short site2);
	virtual void clearpossibleBindingsstencil();
    virtual int numBindingSitesstencil();
    virtual unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>>
    getpossibleBindingsstencil(){
        return _possibleBindingsstencil;
    }
    vector<tuple<CCylinder*, short>> chooseBindingSitesstencil();
    virtual void removePossibleBindingsstencil(CCylinder* cc, short bindingSite);
    virtual void removePossibleBindingsstencil(CCylinder* cc);
    virtual void crosscheck();
    virtual void replacepairPossibleBindingsstencil
            (unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>> _pbs){
//        _possibleBindingsstencil.clear();
        _possibleBindingsstencil = _pbs;
        int oldN = _bindingSpecies->getN();
        int newN = numBindingSitesstencil();
        updateBindingReaction(oldN, newN);
    }
    virtual void setHNLID(short id, short idvec[2]){
        HNLID = id;
        _Hbsmidvec[0] = idvec[0];
        _Hbsmidvec[1] = idvec[1];
    };
    virtual void printbindingsitesstencil();
#endif

#ifdef CUDAACCL_NL
    floatingpoint *gpu_rminmax = NULL;
    int *gpu_numpairs = NULL;
    void assigncudavars();
    void freecudavars();
    floatingpoint* getdistancesCUDA() { return gpu_rminmax;}
    int* getpossiblebindingssizeCUDA(){ return gpu_numpairs;}
    int getNLsize(){
        return _neighborLists[_nlIndex]->getNLsize();
    }
    int* getNLCUDA(){
        return _neighborLists[_nlIndex]->getNLCUDA();
    }
    int* getNLsizeCUDA(){
        return  _neighborLists[_nlIndex]->getNLsizeCUDA();
    }
#endif
};

} // namespace medyan

#endif

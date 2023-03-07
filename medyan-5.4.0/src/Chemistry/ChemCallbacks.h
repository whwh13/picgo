
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

#ifndef MEDYAN_ChemCallbacks_h
#define MEDYAN_ChemCallbacks_h

#include "common.h"
#include "utility.h"

#include "SubSystem.h"
#include "CompartmentGrid.h"
#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"
#include "Bubble.h"
#include "Linker.h"
#include "MotorGhost.h"
#include "BranchingPoint.h"
#include "Boundary.h"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Structure/SurfaceMesh/MembraneRegion.hpp"

#include "BindingManager.h"

#include "Controller/GController.h"
#include "MathFunctions.h"
#include "SysParams.h"
#include "Rand.h"
#include "DissipationTracker.h"

#include <chrono>
#include "CUDAcommon.h"
#include "Util/Io/Log.hpp"


namespace medyan {
using namespace mathfunc;

/// @note - A NOTE TO ALL DEVELOPERS:
///
///         When creating a RSpecies or ReactionBase callback, be sure to not include
///         in the callback structure any variables that may change in the duration
///         of the simulation (i.e. compartments, binding managers, etc). This information
///         needs to be found dynamically to ensure correct behavior - the alternative may
///         cause brutal bugs that are difficult to track in simulation.
///

/// Callback to update the compartment-local binding species based on
/// a change of copy number for an empty site.
struct UpdateBrancherBindingCallback {
    
    Cylinder* _cylinder; ///< cylinder to update
    
    short _bindingSite;  ///< binding site to update

    chrono::high_resolution_clock::time_point mins, mine;

    //Constructor, sets members
    UpdateBrancherBindingCallback(Cylinder* cylinder, short bindingSite)
    
    : _cylinder(cylinder), _bindingSite(bindingSite) {}
    
    //callback
    void operator() (RSpecies *r, int delta) {
	    mins = chrono::high_resolution_clock::now();

        //update this cylinder
        Compartment* c = _cylinder->getCompartment();
        
        for(auto &manager : c->getFilamentBindingManagers()) {
            
            if(dynamic_cast<BranchingManager*>(manager.get())) {
                
                CCylinder* cc = _cylinder->getCCylinder();
                
                //update binding sites
                if(delta == +1) {
#ifdef NLORIGINAL
                    manager->addPossibleBindings(cc, _bindingSite);
#endif
#if defined(NLSTENCILLIST) || defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
                    manager->addPossibleBindingsstencil(cc, _bindingSite);
#endif
                }
                
                else /* -1 */{
#ifdef NLORIGINAL
                    manager->removePossibleBindings(cc, _bindingSite);
#endif
#if defined(NLSTENCILLIST) || defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
                    manager->removePossibleBindingsstencil(cc, _bindingSite);
#endif
                }
            }
        }
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
		CUDAcommon::ctime.tUpdateBrancherBindingCallback += elapsed_time.count();
	    CUDAcommon::ccount.cUpdateBrancherBindingCallback++;
    }
};


/// In scenarios where you have multiple linkers in the system, they cannot bind to the
// same binding site simultaneously. If linker type x is bound to a site, no other linker
// (of any type) can bind to it till this one unbinds.
struct UpdateLinkerBindingCallback {
    
    Cylinder* _cylinder; ///< cylinder to update
    
    short _bindingSite;  ///< binding site to update

    chrono::high_resolution_clock::time_point mins, mine;

    //Constructor, sets members
    UpdateLinkerBindingCallback(Cylinder* cylinder, short bindingSite)
    
    : _cylinder(cylinder), _bindingSite(bindingSite) {}
    
    //callback
    void operator() (RSpecies *r, int delta) {
	    mins = chrono::high_resolution_clock::now();

        //update this cylinder
        Compartment* c = _cylinder->getCompartment();
        for(auto &manager : c->getFilamentBindingManagers()) {
            
            if(dynamic_cast<LinkerBindingManager*>(manager.get())) {
                
                CCylinder* cc = _cylinder->getCCylinder();
                
                //update binding sites
                if(delta == +1) {
#ifdef NLORIGINAL
                    manager->addPossibleBindings(cc, _bindingSite);
#endif
#if defined(NLSTENCILLIST) || defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
                    manager->addPossibleBindingsstencil(cc, _bindingSite);
#endif
                }

                else /* -1 */{
#ifdef NLORIGINAL
                    manager->removePossibleBindings(cc, _bindingSite);
#endif
#if defined(NLSTENCILLIST) || defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
                    manager->removePossibleBindingsstencil(cc, _bindingSite);
#endif
                }
            }
        }
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tUpdateLinkerBindingCallback += elapsed_time.count();
	    CUDAcommon::ccount.cUpdateLinkerBindingCallback++;
    }
};

struct UpdateMotorBindingCallback {
    
    Cylinder* _cylinder; ///< cylinder to update
    
    short _bindingSite;  ///< binding site to update

    chrono::high_resolution_clock::time_point mins, mine;

    //Constructor, sets members
    UpdateMotorBindingCallback(Cylinder* cylinder, short bindingSite)
    
    : _cylinder(cylinder), _bindingSite(bindingSite) {}
    
    //callback
    void operator() (RSpecies *r, int delta) {
	    mins = chrono::high_resolution_clock::now();
        //update this cylinder
        Compartment* c = _cylinder->getCompartment();
        
        for(auto &manager : c->getFilamentBindingManagers()) {
            
            if(dynamic_cast<MotorBindingManager*>(manager.get())) {
                
                CCylinder* cc = _cylinder->getCCylinder();
                
                //update binding sites
                if(delta == +1) {
#ifdef NLORIGINAL
                    manager->addPossibleBindings(cc, _bindingSite);
#endif
#if defined(NLSTENCILLIST) || defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
                    manager->addPossibleBindingsstencil(cc, _bindingSite);
#endif
                }

                else /* -1 */{
#ifdef NLORIGINAL
                    manager->removePossibleBindings(cc, _bindingSite);
#endif
#if defined(NLSTENCILLIST) || defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
                    manager->removePossibleBindingsstencil(cc, _bindingSite);
#endif
                }
            }
        }
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tUpdateMotorBindingCallback += elapsed_time.count();
	    CUDAcommon::ccount.cUpdateMotorBindingCallback++;
    }
};

struct UpdateMotorIDCallback{
    
    int _motorType; ///< type of motor to find its binding manager

    chrono::high_resolution_clock::time_point mins, mine;

    //Constructor, sets members
    UpdateMotorIDCallback(int motorType) : _motorType(motorType) {};
    
    //callback
    void operator() (RSpecies *r, int delta) {
	    mins = chrono::high_resolution_clock::now();
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tUpdateMotorIDCallback += elapsed_time.count();
	    CUDAcommon::ccount.cUpdateMotorIDCallback++;
    }
};




/// Callback to extend the plus end of a Filament after a polymerization
/// Reaction occurs in the system.
struct FilamentExtensionPlusEndCallback {
    
    Cylinder* _cylinder;
    
    short _plusEnd; ///< Plus end species to mark

    chrono::high_resolution_clock::time_point mins, mine;

    //Constructor, sets members
    FilamentExtensionPlusEndCallback(Cylinder* cylinder, short plusEnd)
    : _cylinder(cylinder), _plusEnd(plusEnd){};
    
    //Callback
    void operator() (ReactionBase *r){
	    mins = chrono::high_resolution_clock::now();
        //extend the front
        Filament* f = (Filament*)(_cylinder->getParent());
        f->extendPlusEnd(_plusEnd);
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tFilamentExtensionPlusEndCallback += elapsed_time.count();
	    CUDAcommon::ccount.cFilamentExtensionPlusEndCallback++;
    }
};

/// Callback to extend the minus end of a Filament after a polymerization
/// Reaction occurs in the system.
struct FilamentExtensionMinusEndCallback {
    
    Cylinder* _cylinder;
    
    short _minusEnd; ///< Minus end species to mark

    chrono::high_resolution_clock::time_point mins, mine;

    //Constructor, sets members
    FilamentExtensionMinusEndCallback(Cylinder* cylinder, short minusEnd)
    : _cylinder(cylinder), _minusEnd(minusEnd){};
    //Callback
    void operator() (ReactionBase *r){
	    mins = chrono::high_resolution_clock::now();
        //extend the back
        Filament* f = (Filament*)(_cylinder->getParent());
        f->extendMinusEnd(_minusEnd);
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tFilamentExtensionMinusEndCallback += elapsed_time.count();
	    CUDAcommon::ccount.cFilamentExtensionMinusEndCallback++;
    }
};

/// Callback to retract the plus end of a Filament after a depolymerization
/// Reaction occurs in the system.
struct FilamentRetractionPlusEndCallback {
    
    Cylinder* _cylinder;

    //Constructor, sets members
    FilamentRetractionPlusEndCallback(Cylinder* cylinder)
    : _cylinder(cylinder) {}
    // Callback
    void operator() (ReactionBase *r) const {
	    const auto mins = chrono::high_resolution_clock::now();
        Filament* f =(Filament*)( _cylinder->getParent());
        f->retractPlusEnd();

        // Warning: beyond this point (after retraction), this callback object is deleted with the corresponding reaction. Any further access to member variables results in UB.
	    const auto mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tFilamentRetractionPlusEndCallback += elapsed_time.count();
	    CUDAcommon::ccount.cFilamentRetractionPlusEndCallback++;
    }
};

/// Callback to retract the minus end of a Filament after a depolymerization
/// Reaction occurs in the system.
struct FilamentRetractionMinusEndCallback {
    
    Cylinder* _cylinder;

    //Constructor, sets members
    FilamentRetractionMinusEndCallback(Cylinder* cylinder)
    : _cylinder(cylinder) {}
    // Callback
    void operator() (ReactionBase *r) const {
	    const auto mins = chrono::high_resolution_clock::now();
        Filament* f = (Filament*)(_cylinder->getParent());
        f->retractMinusEnd();

        // Warning: beyond this point (after retraction), this callback object is deleted with the corresponding reaction. Any further access to member variables results in UB.
	    const auto mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tFilamentRetractionMinusEndCallback += elapsed_time.count();
	    CUDAcommon::ccount.cFilamentRetractionMinusEndCallback++;
    }
};

/// Callback to polymerize the plus end of a Filament after a polymerization
/// Reaction occurs in the system.
struct FilamentPolymerizationPlusEndCallback {
    
    Cylinder* _cylinder;

    chrono::high_resolution_clock::time_point mins, mine;

    //Constructor, sets members
    FilamentPolymerizationPlusEndCallback(Cylinder* cylinder)
    : _cylinder(cylinder) {};
    //Callback
    void operator() (ReactionBase *r){
	    mins = chrono::high_resolution_clock::now();
        Filament* f = (Filament*)(_cylinder->getParent());
        f->polymerizePlusEnd();
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tFilamentPolymerizationPlusEndCallback += elapsed_time.count();
	    CUDAcommon::ccount.cFilamentPolymerizationPlusEndCallback++;
    }
};

/// Callback to polymerize the minus end of a Filament after a polymerization
/// Reaction occurs in the system.
struct FilamentPolymerizationMinusEndCallback {
    
    Cylinder* _cylinder;

    chrono::high_resolution_clock::time_point mins, mine;

    //Constructor, sets members
    FilamentPolymerizationMinusEndCallback(Cylinder* cylinder)
    : _cylinder(cylinder) {};
    //Callback
    void operator() (ReactionBase *r){
	    mins = chrono::high_resolution_clock::now();
        Filament* f = (Filament*)(_cylinder->getParent());
        f->polymerizeMinusEnd();
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tFilamentPolymerizationMinusEndCallback += elapsed_time.count();
	    CUDAcommon::ccount.cFilamentPolymerizationMinusEndCallback++;
    }
};

/// Callback to depolymerize the plus end of a Filament after a depolymerization
/// Reaction occurs in the system.
struct FilamentDepolymerizationPlusEndCallback {
    
    Cylinder* _cylinder;

    chrono::high_resolution_clock::time_point mins, mine;

    //Constructor, sets members
    FilamentDepolymerizationPlusEndCallback(Cylinder* cylinder)
    : _cylinder(cylinder) {};
    //Callback
    void operator() (ReactionBase *r){
	    mins = chrono::high_resolution_clock::now();
        Filament* f = (Filament*)(_cylinder->getParent());
        
        f->depolymerizePlusEnd();
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tFilamentDepolymerizationPlusEndCallback += elapsed_time.count();
	    CUDAcommon::ccount.cFilamentDepolymerizationPlusEndCallback++;
    }
};

/// Callback to depolymerize the back of a Filament after a depolymerization
/// Reaction occurs in the system.
struct FilamentDepolymerizationMinusEndCallback {
    
    Cylinder* _cylinder;

    chrono::high_resolution_clock::time_point mins, mine;

    //Constructor, sets members
    FilamentDepolymerizationMinusEndCallback(Cylinder* cylinder)
    : _cylinder(cylinder) {};
    //Callback
    void operator() (ReactionBase *r){
	    mins = chrono::high_resolution_clock::now();
        Filament* f = (Filament*)(_cylinder->getParent());
        f->depolymerizeMinusEnd();
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tFilamentDepolymerizationMinusEndCallback += elapsed_time.count();
	    CUDAcommon::ccount.cFilamentDepolymerizationMinusEndCallback++;
    }
};

/// Callback to unbind a BranchingPoint from a Filament
struct BranchingPointUnbindingCallback {
    
    SubSystem* _ps;
    BranchingPoint* _branchingPoint;

    BranchingPointUnbindingCallback(BranchingPoint* b, SubSystem* ps)
    : _ps(ps), _branchingPoint(b) {}
    
    void operator() (ReactionBase *r) const {
	    const auto mins = chrono::high_resolution_clock::now();

        //remove the branching point
        _ps->removeTrackable<BranchingPoint>(_branchingPoint);
        delete _branchingPoint;

        // Warning: beyond this point (after removing branching point), this callback object is deleted with the corresponding reaction. Any further access to member variables results in UB.
	    const auto mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tBranchingPointUnbindingCallback += elapsed_time.count();
	    CUDAcommon::ccount.cBranchingPointUnbindingCallback++;
    }
    
    
};


/// Callback to create a BranchingPoint on a Filament
struct BranchingCallback {
    
    SubSystem* _ps;        ///< ptr to subsystem
    
    BranchingManager* _bManager; ///< Branching manager for this compartment
    
    short _plusEnd;        ///< Plus end marker of new cylinder
    
    float _onRate;         ///< Rate of the binding reaction
    float _offRate;        ///< Rate of the unbinding reaction

    chrono::high_resolution_clock::time_point mins, mine;

    BranchingCallback(BranchingManager* bManager, short plusEnd,
                      float onRate, float offRate, SubSystem* ps)
    
    : _ps(ps), _bManager(bManager),
    _plusEnd(plusEnd), _onRate(onRate), _offRate(offRate) {}
    
    void operator() (ReactionBase *r) {
	    mins = chrono::high_resolution_clock::now();
        BranchingPoint* b;
        float frate;
        short branchType = _bManager->getBoundInt();
        
        //choose a random binding site from manager
        tuple<CCylinder*, short> site;
#ifdef NLORIGINAL
        site = _bManager->chooseBindingSite();
#endif
#if defined(NLSTENCILLIST) || defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
        site = _bManager->chooseBindingSitesstencil();
#endif
        
        //get info of mother filament from site
        Cylinder* c1 = get<0>(site)->getCylinder();
        short filType = c1->getType();
        
        floatingpoint pos = floatingpoint(get<1>(site)) / SysParams::Geometry().cylinderNumMon[filType];

        if(SysParams::RUNSTATE==true) {

            //Get a position and direction of a new filament
            auto x1 = c1->getFirstBead()->vcoordinate();
            auto x2 = c1->getSecondBead()->vcoordinate();
            
            //get original direction of cylinder
            auto p= midPointCoordinate(x1, x2, pos);
            vector<floatingpoint> n = twoPointDirection(x1, x2);
            
            //get branch projection
            //use mechanical parameters
            floatingpoint l, t;
            if(SysParams::Mechanics().BrStretchingL.size() != 0) {
                l = SysParams::Mechanics().BrStretchingL[branchType];
                t = SysParams::Mechanics().BrBendingTheta[branchType];
            }
            else {
                cout << "Branching initialization cannot occur unless mechanical parameters are specified."
                << " Using default values for Arp2/3 complex - l=10.0nm, theta=70.7deg"
                << endl;
                l = 10.0;
                t = 1.22;
            }
            floatingpoint s = SysParams::Geometry().monomerSize[filType];

            //n direction vector of mother filament
            //p coordinate of brancher
            //l stretching equilibrium length
            //s monomer size
            //t bending theta.
            int tries = 0;
            constexpr int triesShiftParam = 10;
            constexpr int triesWarning    = 10;
            bool inboundary = false;
            frate=_offRate;
            float cutofffactor = pow(c1->getCompartment()->getVolumeFrac(),1.0/3.0);
            floatingpoint boundary_cutoff_distance = cutofffactor*SysParams::Boundaries()
                    .BoundaryCutoff/4.0 ;

            while(inboundary == false) {

                // We try to find a point at t radians with mother filament. We have triesShiftParam
                //number of trials to find an appropriate one. If not, we choose a random
                //theta and let minimization algorithm take care of bringing it back to t
                // radians.
                const double theta = (tries >= triesShiftParam ? Rand::randfloatingpoint(0.1, 3.04) : t);
                auto branchPosDir = branchProjection(n, p, l, s, theta);
                auto bd = get<0>(branchPosDir);//branch direction
                auto bp = get<1>(branchPosDir);//branch position
                const auto b2 = vector2Vec< 3, floatingpoint >(bp) + s * vector2Vec< 3, floatingpoint >(bd);

                //Check if the branch will be within boundary
                auto projlength = SysParams::Geometry().cylinderSize[filType] / 10;
                auto pos2 = nextPointProjection(bp, projlength, bd);
                //check if within cutoff of boundary
                auto regionInMembrane = _ps->getRegionInMembrane();
                if(
                    (regionInMembrane && (
                        regionInMembrane->contains(*_ps, vector2Vec< 3, floatingpoint >(bp)) &&
                        regionInMembrane->contains(*_ps, b2)
                    )) ||
                    (!regionInMembrane && (
                        _ps->getBoundary()->distance(bp) >= boundary_cutoff_distance &&
                        _ps->getBoundary()->distance(pos2) >= boundary_cutoff_distance
                    ))
                ) {
                    inboundary = true;
                }

                ++tries;
                if(tries >= triesWarning && inboundary == false) {
                    LOG(WARNING) << "Cannot find a branching point in region. Trial " << tries << ": "
                                 << "dir (" << bd[0] << ' ' << bd[1] << ' ' << bd[2] << ") "
                                 << "pos (" << bp[0] << ' ' << bp[1] << ' ' << bp[2] << ')';
                    inboundary = true;
                }

                if(inboundary) {
                    //create a new daughter filament
                    Filament *f = _ps->addTrackable<Filament>(_ps, filType, bp, bd, true, true);

                    //mark first cylinder
                    Cylinder *c = f->getCylinderVector().front();
                    c->getCCylinder()->getCMonomer(0)->speciesPlusEnd(_plusEnd)->up();

                    //create new branch
                    b = _ps->addTrackable<BranchingPoint>(c1, c, branchType, pos);
                }
            }

        }
        else {
            CCylinder* c = nullptr;
            bool check = false;
            vector<tuple<tuple<CCylinder*, short>, tuple<CCylinder*, short>>> BrT=_bManager->getbtuple();
/*	        cout<<"Looking for cylinder with Idx "<<c1->getStableIndex()<<" and pos "
                                                                   ""<<pos<<endl;*/
            for(auto T:BrT){
            	//Mother cylinder
                CCylinder* cx=get<0>(get<0>(T));
                //Mother binding site
                floatingpoint p = floatingpoint(get<1>(get<0>(T)))/ floatingpoint(SysParams::Geometry().cylinderNumMon[filType]);
//	            cout<<"Branch tuple "<<cx->getCylinder()->getStableIndex()<<" "<<p<<endl;
                if(cx->getCylinder()->getId()==c1->getId() && abs(p-pos)/pos < 0.01){
                    c=get<0>(get<1>(T));
                    check = true;
                    break;
                }}
            if(check){
            b= _ps->addTrackable<BranchingPoint>(c1, c->getCylinder(), branchType, pos);
                
            frate=0.0;
            }
            else {
                b = nullptr;
                LOG(ERROR)<<"Brancher Error. Cannot find binding Site in the list. Cannot "
					  "complete restart. Exiting." <<endl;
                exit(EXIT_FAILURE);
            }
        }
        
        //create off reaction
        auto cBrancher = b->getCBranchingPoint();
        
        cBrancher->setRates(_onRate, frate);
        cBrancher->createOffReaction(r, _ps);
        cBrancher->getOffReaction()->setBareRate(SysParams::BUBBareRate[branchType]);

        b -> updateReactionRates();

	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tBranchingCallback += elapsed_time.count();
	    CUDAcommon::ccount.cBranchingCallback++;
    }
};

/// Callback to unbind a Linker from a Filament
struct LinkerUnbindingCallback {
    
    SubSystem* _ps;
    Linker* _linker;

    LinkerUnbindingCallback(Linker* l, SubSystem* ps) : _ps(ps), _linker(l) {}
    
    void operator() (ReactionBase *r) const {
#ifdef OPTIMOUT
	    CUDAcommon::tmin.linkerunbindingcalls++;
#endif
	    const auto mins = chrono::high_resolution_clock::now();

        //remove the linker
        _ps->removeTrackable<Linker>(_linker);
        delete _linker;
	    const auto mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tLinkerUnbindingCallback += elapsed_time.count();
	    CUDAcommon::ccount.cLinkerUnbindingCallback++;
    }
};

/// Callback to bind a Linker to Filament
struct LinkerBindingCallback {
    
    SubSystem* _ps;               ///< ptr to subsystem
    
    LinkerBindingManager* _lManager; ///< Linker binding manager for this compartment
    
    float _onRate;                ///< Rate of the binding reaction
    float _offRate;               ///< Rate of the unbinding reaction

    DissipationTracker* _dt;

    chrono::high_resolution_clock::time_point mins, mine;

    LinkerBindingCallback(LinkerBindingManager* lManager,
                          float onRate, float offRate, SubSystem* ps, DissipationTracker* dt)
    
    : _ps(ps), _lManager(lManager), _onRate(onRate), _offRate(offRate), _dt(dt) {}
    
    void operator() (ReactionBase *r) {
#ifdef OPTIMOUT
	    CUDAcommon::tmin.linkerbindingcalls++;
#endif
	    mins = chrono::high_resolution_clock::now();

        //get a random binding
        short linkerType = _lManager->getBoundInt();
        auto& linkerSpeciesIndices = _lManager->getLinkerSpeciesIndices();
        float f;
        
        //choose a random binding site from manager
        vector<tuple<CCylinder*, short>> site;
#ifdef NLORIGINAL
        site = _lManager->chooseBindingSites();
#endif
#if defined(NLSTENCILLIST) || defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
        site = _lManager->chooseBindingSitesstencil();
#endif
        
        Cylinder* c1 = get<0>(site[0])->getCylinder();
        Cylinder* c2 = get<0>(site[1])->getCylinder();
        
        const auto filType1 = c1->getType();
        const auto filType2 = c2->getType();
        
        // Create a linker
        const int cylinderNumMon1 = SysParams::Geometry().cylinderNumMon[filType1];
        const int cylinderNumMon2 = SysParams::Geometry().cylinderNumMon[filType2];
        
        const auto pos1 = floatingpoint(get<1>(site[0])) / cylinderNumMon1;
        const auto pos2 = floatingpoint(get<1>(site[1])) / cylinderNumMon2;
        
        Linker* l = _ps->addTrackable<Linker>(c1, c2, linkerType, linkerSpeciesIndices[0], linkerSpeciesIndices[1], pos1, pos2);
        
        //create off reaction
        auto cLinker = l->getCLinker();
        f=_offRate;
        //@
        cLinker->setRates(_onRate, f);//offRate during restart is controlled by
        // RESTARTPHASESWITCH option in setRateMulFactor
        cLinker->createOffReaction(r, _ps);
        
        if(SysParams::Chemistry().eventTracking){
            _dt->recordLinkerBinding(l);
        }

        l->updateReactionRates();

	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tLinkerBindingCallback += elapsed_time.count();
	    CUDAcommon::ccount.cLinkerBindingCallback++;
    }
};

/// Callback to unbind a MotorGhost from a Filament
struct MotorUnbindingCallback {
    
    SubSystem* _ps;
    MotorGhost* _motor;

    MotorUnbindingCallback(MotorGhost* m, SubSystem* ps) :
    
    _ps(ps), _motor(m) {}
    
    void operator() (ReactionBase *r) const {
#ifdef OPTIMOUT
    	CUDAcommon::tmin.motorunbindingcalls++;
#endif
	    const auto mins = chrono::high_resolution_clock::now();

        //remove the motor
        _ps->removeTrackable<MotorGhost>(_motor);
        delete _motor;
	    const auto mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tMotorUnbindingCallback += elapsed_time.count();
	    CUDAcommon::ccount.cMotorUnbindingCallback++;
    }
};


/// Callback to bind a MotorGhost to Filament
struct MotorBindingCallback {
    
    SubSystem* _ps;               ///< ptr to subsystem
    
    MotorBindingManager* _mManager;///< Motor binding manager for this compartment
    
    float _onRate;                ///< Rate of the binding reaction
    float _offRate;               ///< Rate of the unbinding reaction

    chrono::high_resolution_clock::time_point mins, mine;

    MotorBindingCallback(MotorBindingManager* mManager,
                         float onRate, float offRate, SubSystem* ps)
    
    : _ps(ps), _mManager(mManager), _onRate(onRate), _offRate(offRate) {}
    
    void operator() (ReactionBase *r) {
#ifdef OPTIMOUT
	    CUDAcommon::tmin.motorbindingcalls++;
#endif
	    mins = chrono::high_resolution_clock::now();

        //get a random binding
        short motorType = _mManager->getBoundInt();
        auto& motorSpeciesIndices = _mManager->getMotorSpeciesIndices();
        float f;
        //choose a random binding site from manager
        vector<tuple<CCylinder*, short>> site;
#ifdef NLORIGINAL
        site = _mManager->chooseBindingSites();
#endif
#if defined(NLSTENCILLIST) || defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
        site = _mManager->chooseBindingSitesstencil();
#endif
        
        Cylinder* c1 = get<0>(site[0])->getCylinder();
        Cylinder* c2 = get<0>(site[1])->getCylinder();
        
        const auto filType1 = c1->getType();
        const auto filType2 = c2->getType();
        
        // Create a motor
        const int cylinderNumMon1 = SysParams::Geometry().cylinderNumMon[filType1];
        const int cylinderNumMon2 = SysParams::Geometry().cylinderNumMon[filType2];
        
        const auto pos1 = floatingpoint(get<1>(site[0])) / cylinderNumMon1;
        const auto pos2 = floatingpoint(get<1>(site[1])) / cylinderNumMon2;

        
        MotorGhost* m = _ps->addTrackable<MotorGhost>(c1, c2, motorType, motorSpeciesIndices[0], motorSpeciesIndices[1], pos1, pos2, _onRate, _offRate);
        
        //attach an ID to the motor based on the transfer ID
        //DEPRECATED AS OF 9/22/16
//        m->setID(MotorGhost::_motorGhosts.getTransferID());
        
        //create off reaction
        auto cMotorGhost = m->getCMotorGhost();
        //aravind June 24, 2016.
//        if(SysParams::RUNSTATE==false){
//        f=0.0;
//        }
//        else
//            f=_offRate;
        //@
        f = _offRate;
        cMotorGhost->setRates(_onRate, f);
        cMotorGhost->createOffReaction(r, _ps);
        
        //reset the associated walking reactions
        m->updateReactionRates();

	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tMotorBindingCallback += elapsed_time.count();
	    CUDAcommon::ccount.cMotorBindingCallback++;

    }
};

/// Callback to walk a MotorGhost on a Filament
struct MotorWalkingCallback {
    
    Cylinder* _c;        ///< Cylinder this callback is attached to
    
    short _oldPosition;  ///< Old position of motor head
    short _newPosition;  ///< New position of motor head
    
    int speciesMotorIndex; // Index of the motor species in the corresponding filament.
    short _boundType;    ///< Type of bound this motor took place of
    
    SubSystem* _ps;      ///< Ptr to subsystem

    chrono::high_resolution_clock::time_point mins, mine;
    DissipationTracker* _dt;

    MotorWalkingCallback(Cylinder* c,
                         short oldPosition, short newPosition,
                         short speciesMotorIndex, short boundType, SubSystem* ps, DissipationTracker* dt)
    
    :_c(c), _oldPosition(oldPosition), _newPosition(newPosition),
    speciesMotorIndex(speciesMotorIndex), _boundType(boundType), _ps(ps), _dt(dt) {}
    
    void operator() (ReactionBase* r) {
#ifdef OPTIMOUT
	    CUDAcommon::tmin.motorwalkingcalls++;
#endif
	    mins = chrono::high_resolution_clock::now();
//        cout<<"Motor walking begins"<<endl;
        //get species
        CCylinder* cc = _c->getCCylinder();
        CMonomer* monomer = cc->getCMonomer(_oldPosition);
        SpeciesBound* sm1 = monomer->speciesMotor(speciesMotorIndex);
        
        short filType = _c->getType();
        
        //get motor
        MotorGhost* m = ((CMotorGhost*)sm1->getCBound())->getMotorGhost();

        int cylinderSize = SysParams::Geometry().cylinderNumMon[filType];

        floatingpoint oldpos = floatingpoint(_oldPosition) / cylinderSize;
        floatingpoint newpos = floatingpoint(_newPosition) / cylinderSize;
	    m->moveMotorHead(_c, oldpos, newpos, speciesMotorIndex, _boundType, _ps);
        
        //reset the associated reactions
        m->updateReactionRates();

        
        if(SysParams::Chemistry().eventTracking){
            _dt->recordWalk(m);
        }
        
        
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tMotorWalkingCallback += elapsed_time.count();
	    CUDAcommon::ccount.cMotorWalkingCallback++;
    }
};

/// Callback to walk a MotorGhost on a Filament to a new Cylinder
struct MotorMovingCylinderCallback {
    
    Cylinder* _oldC;        ///< Old cylinder the motor is attached to
    Cylinder* _newC;        ///< New cylinder motor will be attached to
    
    short _oldPosition;     ///< Old position of motor head
    short _newPosition;     ///< New position of motor head
    
    int speciesMotorIndex;  // Index of the motor species in the corresponding filament.
    short _boundType;       ///< Type of bound this motor is taking place of
    
    SubSystem* _ps;         ///< Ptr to subsystem

    chrono::high_resolution_clock::time_point mins, mine;

    MotorMovingCylinderCallback(Cylinder* oldC, Cylinder* newC,
                                short oldPosition, short newPosition,
                                short speciesMotorIndex, short boundType, SubSystem* ps)
    
    :_oldC(oldC), _newC(newC), _oldPosition(oldPosition), _newPosition(newPosition),
    speciesMotorIndex(speciesMotorIndex), _boundType(boundType), _ps(ps) {}
    
    void operator() (ReactionBase* r) {
#ifdef OPTIMOUT
	    CUDAcommon::tmin.motorwalkingcalls++;
#endif

	    mins = chrono::high_resolution_clock::now();
//        cout<<"Motor moving cylinder begins"<<endl;
        //get species
        CCylinder* oldCC = _oldC->getCCylinder();
        CMonomer* monomer = oldCC->getCMonomer(_oldPosition);
        SpeciesBound* sm1 = monomer->speciesMotor(speciesMotorIndex);
        short filType = _oldC->getType();

        //get motor
        MotorGhost* m = ((CMotorGhost*)sm1->getCBound())->getMotorGhost();
        
        int cylinderSize = SysParams::Geometry().cylinderNumMon[filType];
/*        cout<<"filament Type "<<filType<<endl;
        cout<<"cylinder size "<<cylinderSize<<endl;
        cout<<"Cylinder oldpos "<<_oldC->getID()<<" "<<_oldPosition<<endl;
        cout<<"Cylinder newpos "<<_newC->getID()<<" "<<_newPosition<<endl;
        cout<<"-----------"<<endl;*/
        floatingpoint oldpos = floatingpoint(_oldPosition) / cylinderSize;
        floatingpoint newpos = floatingpoint(_newPosition) / cylinderSize;
        
        m->moveMotorHead(_oldC, _newC, oldpos, newpos, speciesMotorIndex, _boundType, _ps);
        
        //reset the associated reactions
        m->updateReactionRates();

	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tMotorMovingCylinderCallback += elapsed_time.count();
	    CUDAcommon::ccount.cMotorMovingCylinderCallback++;
    }
};

/// Struct to create a new filament based on a given reaction
struct FilamentCreationCallback {
    
    //@{
    /// Integer specifying the chemical properties of first cylinder
    short _plusEnd;
    short _minusEnd;
    short _filament;
    //@}
    
    ///Filament type to create
    short _filType;
    
    Compartment* _compartment; ///< compartment to put this filament in
    SubSystem* _ps; ///< Ptr to the subsystem

    chrono::high_resolution_clock::time_point mins, mine;

    FilamentCreationCallback(short plusEnd, short minusEnd, short filament,
                             short filType, SubSystem* ps, Compartment* c = nullptr)
    
    : _plusEnd(plusEnd), _minusEnd(minusEnd), _filament(filament),
    _filType(filType), _compartment(c), _ps(ps) {}
    
    void operator() (ReactionBase* r) {
	    mins = chrono::high_resolution_clock::now();

        Compartment* c;
        
        //no compartment was set, pick a random one
        if(_compartment == nullptr)
            c = &GController::getRandomCompartment();
        else
            c = _compartment;
        
        //set up a random initial position and direction
        vector<floatingpoint> position;
        vector<floatingpoint> direction;

        if(_ps->getBoundary()->getShape() == BoundaryShape::Cylinder) {
            while(true) {
                position = mathfunc::vec2Vector(_ps->getCompartmentGrid()->getRandomCenterCoordinatesIn(*c));
                
                //getting random numbers between -1 and 1
                
                direction = {Rand::randfloatingpoint(-1,1), Rand::randfloatingpoint(-1,1), 0};
                normalize(direction);
                
                auto npp = nextPointProjection(position,
                                               SysParams::Geometry().cylinderSize[_filType], direction);
                
                bool inbubble = false;
                //check if within bubble
                for(auto& bb : _ps->bubbles) {
                    auto radius = bb.getRadius();

                    if((twoPointDistancesquared(vec2Vector(bb.coord), position) < (radius * radius)) ||
                       (twoPointDistancesquared(vec2Vector(bb.coord), npp) < (radius * radius))){
                        inbubble = true;
                        break;
                    }
                }

                if(inbubble) continue;

                //check if within boundary && if within REPULSIONEXPIN region (set as 125nm)
                if(_ps->getBoundary()->within(position) &&
                   _ps->getBoundary()->within(npp)){
                    if(_ps->getBoundary()->distance(position) > 125 &&
                       _ps->getBoundary()->distance(npp) > 125)
                    	break;
                }

            }
            
            //create filament, set up ends and filament species
            Filament* f = _ps->addTrackable<Filament>(_ps, _filType, position, direction, true, false);
            
            //initialize the nucleation
            f->nucleate(_plusEnd, _filament, _minusEnd);
        }
        else {
            while(true) {
                position = mathfunc::vec2Vector(_ps->getCompartmentGrid()->getRandomCoordinatesIn(*c));
                
                //getting random numbers between -1 and 1
                
                direction = {Rand::randfloatingpoint(-1,1), Rand::randfloatingpoint(-1,1), Rand::randfloatingpoint(-1,1)};
                normalize(direction);
                
                auto npp = nextPointProjection(position,
                                               SysParams::Geometry().cylinderSize[_filType], direction);
                
                bool inbubble = false;
                //check if within bubble
                for(auto& bb : _ps->bubbles) {
                    auto radius = bb.getRadius();

                    if((twoPointDistancesquared(vec2Vector(bb.coord), position) < (radius * radius)) ||
                       (twoPointDistancesquared(vec2Vector(bb.coord), npp) < (radius * radius))){
                        inbubble = true;
                        break;
                    }
                }

                if(inbubble) continue;

                //check if within boundary and membrane[0]
                if(
                    _ps->getBoundary()->within(position) &&
                    _ps->getBoundary()->within(npp) &&
                    c->isActivated() &&
                    (c->getVolumeFrac() >= 1.0 ||
                        (medyan::contains(*_ps, _ps->membranes.begin()->getMesh(), vector2Vec<3, floatingpoint>(position)) &&
                        medyan::contains(*_ps, _ps->membranes.begin()->getMesh(), vector2Vec<3, floatingpoint>(npp))))
                ) break;
            }
            
            //create filament, set up ends and filament species
            Filament* f = _ps->addTrackable<Filament>(_ps, _filType, position, direction, true, false);
            
            //initialize the nucleation
            f->nucleate(_plusEnd, _filament, _minusEnd);
        }
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tFilamentCreationCallback += elapsed_time.count();
	    CUDAcommon::ccount.cFilamentCreationCallback++;
        }

    
};

///Struct to sever a filament based on a reaction
struct FilamentSeveringCallback {
    
    Cylinder* _c1;  ///< Filament severing point

    chrono::high_resolution_clock::time_point mins, mine;

    FilamentSeveringCallback(Cylinder* c1) : _c1(c1) {}
    
    void operator() (ReactionBase* r) {
	    mins = chrono::high_resolution_clock::now();

        //reactants should be re-marked
        for(int i = 0; i < SEVERINGREACTANTS + 1; i++)
            r->rspecies()[i]->up();
        
        //sever the filament at given position
        Filament* f = (Filament*)(_c1->getParent());
        
        //create a new filament by severing
        Filament* newF = f->sever(_c1->getPosition());
        
        //if severing did not occur, return
        if(newF == nullptr) return;
	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tFilamentSeveringCallback += elapsed_time.count();
	    CUDAcommon::ccount.cFilamentSeveringCallback++;
    }
};

/// Struct to destroy a filament based on a reaction
struct FilamentDestructionCallback {
    
    Cylinder* _c; ///< Cylinder to destroy
    
    SubSystem* _ps; ///< SubSystem ptr

    FilamentDestructionCallback(Cylinder* c, SubSystem* ps) : _c(c), _ps(ps) {}
    
    void operator() (ReactionBase* r) {
        const auto mins = chrono::high_resolution_clock::now();

        Filament* f = (Filament*)(_c->getParent());
        
        _ps->removeTrackable<Filament>(f);
        delete f;

        // Warning: beyond this point (after destruction), this callback object is deleted with the corresponding reaction. Any further access to member variables results in UB.
        const auto mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> elapsed_time(mine - mins);
	    CUDAcommon::ctime.tFilamentDestructionCallback += elapsed_time.count();
	    CUDAcommon::ccount.cFilamentDestructionCallback++;
    }
};

} // namespace medyan

#endif

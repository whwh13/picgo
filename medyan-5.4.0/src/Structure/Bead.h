
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

#ifndef MEDYAN_Bead_h
#define MEDYAN_Bead_h

#include <algorithm> // find
#include <vector>
#include <list>

#include "common.h"
#include "MathFunctions.h" // vec2Vector
#include "Structure/Database.h"
#include "Component.h"
#include "Composite.h"
#include "Trackable.h"
#include "Movable.h"
#include "DynamicNeighbor.h"
#include "SysParams.h"
#include "Util/Math/Vec.hpp"

namespace medyan {

//FORWARD DECLARATIONS
class Compartment;
class Filament;

/// Represents a single coordinate between [Cylinders](@ref Cylinder), and holds forces
/// needed for mechanical equilibration.
/*!
 *  Beads are the "hinges" between [Cylinders](@ref Cylinder). In the minimization 
 *  algorithms, beads are moved corresponding to external forces, for example, Filament 
 *  stretching and bending. The bead class contains currernt coordinates and forces, and 
 *  has functions to calculate dot products for the minimization algorithms.
 *
 *  Extending the Movable class, the positions of all instances can 
 *  be updated by the SubSystem.
 *
 *  Extending the DynamicNeighbor class, all instances can be kept in 
 *  [NeighborLists](@ref NeighborList).
 */

class Bead : public Component, public Trackable, public Movable, public DynamicNeighbor,
    public Database< Bead, true > {
    
public:
    using DatabaseType = Database< Bead, true >;

    medyan::Vec< 3, floatingpoint > coord;
    medyan::Vec< 3, floatingpoint > force;

    ///@note - all vectors are in x,y,z coordinates.
    vector<floatingpoint> coordinateP; ///< Prev coordinates of bead in CG minimization

                          ///< Forces should always correspond to current coordinates.
    
    vector<floatingpoint> brforce; //boundary repulsion force
    vector<floatingpoint> pinforce;

    vector<floatingpoint> loadForcesP;
    vector<floatingpoint> loadForcesM;
    ///< The force on this bead due to an external load
    ///< This is not a vector (x,y,z) value, but a list of
    ///< force magnitudes in the direction of polymerization with
    ///< monomer increments (future values).
    ///< These are then used to propagate load forces in between
    ///< mechanical force calculations.
    ///< (Edited 20180216) Different angle between the cylinder and the
    ///< boundary would result in different effective monomer size in the
    ///< calculation of the Brownian Ratchet model. To simply computation, we
    ///< include that factor in our loadForces here. As a clarification, the
    ///< actual physical load force should not have that factor.
    
    short lfip = 0;
    short lfim = 0;  ///< Index which saves which load force to use

    // The monomer serial at (slightly plus side of) this bead in a certain
    // filament.
    //
    //                                 |
    //                      -----+-----|-----+-----
    //  minus end <----       M  |  M  | (M) |  M        ----> plus end
    //                      -----+-----|-----+-----
    //                                 |
    //                                 ^ This bead is indicated by the line.
    //
    // The serial of monomer with parenthesis (M) will be indicated by this
    // variable.
    //
    // If the bead is itself at the plus end, this index indicates the serial
    // number of the next imaginary monomer.
    //
    // This monomer serial, as well as the zero offset information, are managed
    // by the associated filament.
    int monomerSerial = 0;
    
    /// The bead can be pinned to a certain position in the simulation volume.
    /// These parameters describe the pinning. Adding the Bead to the list of pinned
    /// Beads is done by a corresponding special protocol. (see executeSpecialProtocols() in Controller)
    vector<floatingpoint> pinnedPosition;
    
    bool isStatic = false;
    
    ///Main constructor
    Bead (vector<floatingpoint> v, Composite* parent, int position);
    
    ///Default constructor
    Bead(Composite* parent, int position);

    auto& coordinate() { return coord; }
    const auto& coordinate() const { return coord; }
    // Temporary compromise
    auto vcoordinate() const { return mathfunc::vec2Vector(coord); }

    /// Get Compartment
    Compartment* getCompartment() {return _compartment;}
    
    /// Get position
    int getPosition() {return _position;}
    
    /// Get the birth time
    float getBirthTime() {return _birthTime;}
    
    //@{
    /// SubSystem management, inherited from Trackable
    // only takes care of pinned bead removal
    virtual void addToSubSystem() override {}
    virtual void removeFromSubSystem() {
        //remove if pinned
        if(_isPinned) removeAsPinned();
    }
    //@}
    
    /// Get all instances of this class from the SubSystem
    static const vector<Bead*>& getBeads() {
        return getElements();
    }
    
    /// Add this bead as a pinned bead
    void addAsPinned() {
        _isPinned = true;
        _pinnedBeads.push_back(this);
    }
    
    /// Remove this bead as pinned. Will remove from pinnedBeads DB
    /// @note - only usually called upon the destruction of a Bead.
    void removeAsPinned() {
        
        _isPinned = false;
        auto it = std::find(_pinnedBeads.begin(), _pinnedBeads.end(), this);
        if(it != _pinnedBeads.end()) _pinnedBeads.erase(it);
    }
    
    const vector<floatingpoint>& getPinPosition() { return pinnedPosition;}
    // Remove all pinned beads.
    void resetAllPinned() {

        _isPinned = false;
        _pinnedBeads.clear();
    }
    /// Get all pinned beads from subsystem
    static const vector<Bead*>& getPinnedBeads() {
        
        return _pinnedBeads;
    }
    
    bool isPinned() {return _isPinned;}
    
    /// Get the number of beads in this system
    static int numBeads() {
        return getElements().size();
    }
    
    /// Update the position, inherited from Movable
    virtual void updatePosition();
    
    virtual void printSelf()const;
    
    //GetType implementation just returns type of parent
    virtual int getType() {return getParent()->getType();}
    //Aravind return static
    bool getstaticstate() {return isStatic;}
    //Aravind set static
    void setstaticstate(bool index) {isStatic = index;}
    //@{
    /// Auxiliary method for CG minimization
    inline double FDotF() {
        return magnitude2(force);
    }
//    inline double FDotF() {
//        return force1[0]*force1[0] +
//        force1[1]*force1[1] +
//        force1[2]*force1[2];
//    }
    
    //Qin add brFDotbrF
    inline floatingpoint brFDotbrF() {
        return brforce[0]*brforce[0] +
        brforce[1]*brforce[1] +
        brforce[2]*brforce[2];
    }
    //add pinFDotpinF
    inline floatingpoint pinFDotpinF() {
        return pinforce[0]*pinforce[0] +
        pinforce[1]*pinforce[1] +
        pinforce[2]*pinforce[2];
    }
    //@}
    
    ///Helper functions for load forces
    
    floatingpoint getLoadForcesP();
    
    void printLoadForcesP() {
        
        cout << "loadP =";
        
        for (int i = 0; i < loadForcesP.size(); i++) {
            
            cout << " " << loadForcesP[i] << " ";
            
        }
        cout << endl;
    }
    
    floatingpoint getLoadForcesM();
 
    void printLoadForcesM()  {
        
        cout << "loadM =";
        
        for (int i = 0; i < loadForcesM.size(); i++) {
            
            cout << " " << loadForcesM[i] << " ";
            
        }
        cout << endl;
    }

    //To be used exclusively within restart protocol.
    void overrideParentInfo(Composite* parent, int p){
        //set position in the filament
        _position  = p;
        //Add this bead as the child
        parent->addChild(unique_ptr<Component>(this));
    }

private:
    Compartment* _compartment = nullptr; ///< Pointer to the compartment that this bead is in
    
    int _position;     ///< Position on structure
    float _birthTime;  ///< Time of birth
    bool _isPinned = false;
    
    static std::vector<Bead*> _pinnedBeads; ///< Collection of pinned beads in SubSystem
                                         ///< (attached to some element in SubSystem)
};

} // namespace medyan

#endif

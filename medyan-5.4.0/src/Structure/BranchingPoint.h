
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

#ifndef MEDYAN_BranchingPoint_h
#define MEDYAN_BranchingPoint_h

#include "common.h"

#include "MBranchingPoint.h"
#include "CBranchingPoint.h"

#include "Database.h"
#include "Trackable.h"
#include "Movable.h"
#include "Component.h"
#include "Reactable.h"
#include "RateChangerImpl.h"

namespace medyan {
//FORWARD DECLARATIONS
class Compartment;
class Cylinder;

/// A container to store a MBranchingPoint and CBranchingPoint.
/*!
 *  BranchingPoint class is used to manage and store a MBranchingPoint and
 *  CBranchingPoint. Upon intialization, both of these components are created.
 *
 *  Extending the Movable class, the positions of all instances 
 *  can be updated by the SubSystem.
 */
class BranchingPoint : public Component, public Trackable, public Movable, public Reactable,
    public Database< BranchingPoint, false > {
    
    friend class DRController;
    
private:
    unique_ptr<MBranchingPoint> _mBranchingPoint; ///< Pointer to mech branch point
    unique_ptr<CBranchingPoint> _cBranchingPoint; ///< Pointer to chem branch point
    
    Cylinder* _c1; ///< Mother cylinder
    Cylinder* _c2; ///< Branching cylinder
    
    floatingpoint _position;  ///< Position on mother cylinder
    
    short _branchType; ///< Integer specifying the type
    
    float _birthTime;  ///<Birth time
    
    Compartment* _compartment; ///< Where this branch point is
    
    ///Helper to get coordinate
    void updateCoordinate();

    ///For dynamic rate unbinding
    static vector<BranchRateChanger*> _unbindingChangers;

    string diffusingactinspeciesname = "";
    
public:
    vector<floatingpoint> coordinate; ///< coordinate of midpoint,
                               ///< updated with updatePosition()
    
    BranchingPoint(Cylinder* c1, Cylinder* c2,
                   short branchType, floatingpoint position = 0.5);
    virtual ~BranchingPoint() noexcept;
    
    //@{
    ///Get attached cylinder
    Cylinder* getFirstCylinder() const {return _c1;}
    Cylinder* getSecondCylinder() const {return _c2;}
    //@}
    
    /// Set chem branch point
    void setCBranchingPoint(CBranchingPoint* cBranchingPoint) {
        _cBranchingPoint = unique_ptr<CBranchingPoint>(cBranchingPoint);
    }
    /// Get chem branch point
    CBranchingPoint* getCBranchingPoint() {return _cBranchingPoint.get();}
    
    /// Get mech branch point
    MBranchingPoint* getMBranchingPoint() {return _mBranchingPoint.get();}
    
    //@{
    /// Position management
    floatingpoint getPosition() {return _position;}
    void setPosition(floatingpoint position) {_position = position;}
    //@}
    
    //@{
    /// Get branch parameter
    virtual int getType() {return _branchType;}
    auto getType() const { return _branchType; }
    //@}
    
    /// Get compartment
    Compartment* getCompartment() {return _compartment;}
    
    /// Get the birth time
    float getBirthTime() {return _birthTime;}
    
    //@{
    /// SubSystem management, inherited from Trackable
    // Does nothing
    virtual void addToSubSystem() { }
    virtual void removeFromSubSystem() {}
    //@}
    
    /// Get all instances of this class from the SubSystem
    static const vector<BranchingPoint*>& getBranchingPoints() {
        return getElements();
    }
    /// Get the number of branching points in this system
    static int numBranchingPoints() {
        return getElements().size();
    }
    
    virtual void printSelf()const;
    
    /// Count the number of brancher species with a given name in the system
    static species_copy_t countSpecies(const string& name);
    
    /// Update the position, inherited from Movable
    /// @note - changes compartment if needed
    virtual void updatePosition();
    
    //Qin ------
    /// Update the reaction rates, inherited from Reactable
    virtual void updateReactionRates();

    void initializerestart(floatingpoint eqLength){ _mBranchingPoint->initializerestart
                (eqLength);};

    void setdiffusingactinspeciesname(string _diffusingactinspeciesname){
        diffusingactinspeciesname = _diffusingactinspeciesname; }

    string getdiffusingactinspeciesname(){
        return diffusingactinspeciesname; }

};

} // namespace medyan

#endif

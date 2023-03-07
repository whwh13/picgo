
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

#ifndef MEDYAN_Linker_h
#define MEDYAN_Linker_h

#include "common.h"

#include "Composite.h"
#include "CLinker.h"
#include "MLinker.h"

#include "Database.h"
#include "Histogram.h"
#include "Trackable.h"
#include "Movable.h"
#include "Reactable.h"
#include "RateChangerImpl.h"

namespace medyan {
//FORWARD DECLARATIONS
class Cylinder;
class DRController;

/// A container to store a MLinker and CLinker.
/*!
 *  Linker class is used to manage and store a MLinker and CLinker. Upon intialization,
 *  both of these components are created.
 *
 *  Extending the Movable class, the positions of all instances 
 *  can be updated by the SubSystem.
 *
 *  Extending the Reactable class, the reactions associated with all
 *  instances can be updated by the SubSystem.
 */
class Linker : public Component, public Trackable, public Movable, public Reactable,
    public Database< Linker, false > {

friend class DRController;
    
private:
    MLinker mLinker_; // Mechanical information of the linker.
    unique_ptr<CLinker> _cLinker; ///< Pointer to chem linker
    
    Cylinder* _c1; ///< First cylinder the linker is bound to
    Cylinder* _c2; ///< Second cylinder the linker is bound to
    
    floatingpoint _position1; ///< Position on first cylinder
    floatingpoint _position2; ///< Position on second cylinder
    
    short _linkerType; ///< Integer specifying the type
    
    float _birthTime; ///Birth time
    
    Compartment* _compartment; ///< Where this linker is
    
    //@{
    ///Histogram data
    static Histogram* _lifetimes;
    //@}
    
    ///For dynamic rate unbinding
    static vector<LinkerRateChanger*> _unbindingChangers;
    
    ///Helper to get coordinate
    void updateCoordinate();
    
    
public:
    vector<floatingpoint> coordinate;
    ///< coordinate of midpoint, updated with updatePosition()
    
    Linker(
        Cylinder* c1, Cylinder* c2,
        short linkerType,
        int linkerSpeciesIndex1, int linkerSpeciesIndex2,
        floatingpoint position1 = 0.5, floatingpoint position2 = 0.5);
    
    virtual ~Linker() noexcept;
    
    //@{
    ///Get attached cylinder
    Cylinder* getFirstCylinder() const {return _c1;}
    Cylinder* getSecondCylinder() const {return _c2;}
    //@}
    
    /// Set chem linker
    void setCLinker(CLinker* cLinker) {_cLinker = unique_ptr<CLinker>(cLinker);}
    /// Get chem linker
    CLinker* getCLinker() {return _cLinker.get();}
    
    /// Get mech linker
    MLinker* getMLinker() {return &mLinker_;}
    
    //@{
    /// Position management
    floatingpoint getFirstPosition() const {return _position1;}
    void setFirstPosition(floatingpoint position1) {_position1 = position1;}
    
    floatingpoint getSecondPosition() const {return _position2;}
    void setSecondPosition(floatingpoint position2) {_position2 = position2;}
    //@}
    
    //@{
    /// Get linker parameter
    virtual int getType() {return _linkerType;}
    auto getType() const { return _linkerType; }
    //@}
    
    /// Get the birth time
    float getBirthTime() {return _birthTime;}
    
    //@{
    /// SubSystem management, inherited from Trackable
    // Does nothing
    virtual void addToSubSystem() override { }
    virtual void removeFromSubSystem() override {}
    //@}
    
    /// Get all instances of this class from the SubSystem
    static const vector<Linker*>& getLinkers() {
        return getElements();
    }
    /// Get the number of linkers in this system
    static int numLinkers() {
        return getElements().size();
    }
    
    /// Get the lifetimes
    static Histogram* getLifetimes() {return _lifetimes;}
    
    /// Update the position, inherited from Movable
    /// @note - changes compartment if needed
    virtual void updatePosition();
    
    /// Update the reaction rates, inherited from Reactable
    virtual void updateReactionRates();
    
    virtual void printSelf()const;
    
    /// Count the number of linker species with a given name in the system
    static species_copy_t countSpecies(const string& name);

    void initializerestart(floatingpoint eqLength) {
        mLinker_.initializerestart(eqLength);
    }
    
};

} // namespace medyan

#endif 

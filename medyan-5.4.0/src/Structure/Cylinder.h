
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

#ifndef MEDYAN_Cylinder_h
#define MEDYAN_Cylinder_h

#include <iostream>

#include "common.h"

#include "MCylinder.h"
#include "CCylinder.h"
#include "RateChanger.h"

#include "Database.h"
#include "Trackable.h"
#include "Movable.h"
#include "Reactable.h"
#include "DynamicNeighbor.h"
#include "Component.h"
#include "Bead.h"
#include "Structure/CellList.hpp"
#include "Util/Math/Vec.hpp"
#include <fstream>

namespace medyan {

//FORWARD DECLARATIONS
class Filament;
class Compartment;
class Bin;

struct CylinderInfo {
    int filamentId = -1;
    int positionOnFilament = -1;
    int filamentFirstEntry = 0;//This variable is updated in CylinderInfoData.
    int compartmentId = -1;
    std::size_t beadIndices[2];
    medyan::Vec< 3, floatingpoint > coord;
    short type = -1;
    int id = -1;
    CCylinder* chemCylinder;
};

using CylinderInfoVec = std::vector< CylinderInfo >;


/// A container to store a MCylinder and CCylinder.
/*!
 *  Cylinder class is used to manage and store a MCylinder and CCylinder.
 *  Upon intialization, both of these components are created.
 *
 *  Extending the Movable class, the positions of all instances
 *  can be updated by the SubSystem.
 *
 *  Extending the Reactable class, the reactions associated with
 *  all instances can be updated by the SubSystem.
 *
 *  Extending the DynamicNeighbor class, all instances can be
 *  kept in [NeighborLists](@ref NeighborList).
 */
class Cylinder : public Component, public Trackable, public Movable,
                                   public Reactable, public DynamicNeighbor,
                                   public Database< Cylinder, true, CylinderInfoVec > {

friend class CController;
friend class DRController;

private:

    chrono::high_resolution_clock::time_point mins, mine;

    Bead* _b1;  ///< Pointer to the first bead.
    Bead* _b2; ///< Pointer to the end bead.

    unique_ptr<MCylinder> _mCylinder; ///< Pointer to mech cylinder
    unique_ptr<CCylinder> _cCylinder; ///< Pointer to chem cylinder

    bool _plusEnd = false;  ///< If the cylinder is at the plus end
    bool _minusEnd = false; ///< If the cylinder is at the minus end

    short _type; ///< Type of cylinder, either corresponding to Filament or other

	int _position;          ///< Position on structure

    medyan::CellListElementUser< Cylinder*, Compartment* > _cellElement;

    Cylinder* _branchingCylinder = nullptr; ///< ptr to a branching cylinder

    ///For dynamic polymerization rate
    static vector<FilamentRateChanger*> _polyChanger;

    inline static medyan::ChemManager* _chemManager = nullptr; ///< A pointer to the ChemManager,
                                      ///< intiailized by CController

    ///Helper to get coordinate
    void updateCoordinate();


    /// ID of filament
    int _filID;


public:
    using DatabaseType = Database< Cylinder, true, CylinderInfoVec >;

    // Coordinates of midpoint, updated with updatePosition().
    vector<floatingpoint> coordinate;
    vector<Bin*> _binvec; //vector of bins. binID corresponding to each binGrid.
    // Pointer to the Bin used by the hybrid neighbor list.
    Bin* hbin = nullptr;

    static bool setpositionupdatedstate; //Setter to check if position has been updated

    // Update CylinderInfoData using newest information in the system
    static void updateAllData() {
        // Update data for all cylinders
        for(auto c : getCylinders()) c->updateData();
    }
    void updateData(); // Update data for this cylinder. TODO: make it const

    /// Constructor, initializes a cylinder
    Cylinder(Composite* parent, Bead* b1, Bead* b2, short type, int position,
             bool extensionFront = false,
             bool extensionBack  = false,
             bool initialization = false,
             floatingpoint eqLength = -1.0);

    virtual ~Cylinder() noexcept;

    const auto& getCoordinate() const { return getDbData()[getStableIndex()].coord; }
    auto      & getCoordinate()       { return getDbData()[getStableIndex()].coord; }

    /// Get mech cylinder
    MCylinder* getMCylinder() const {return _mCylinder.get();}

    /// Get chem cylinder
    CCylinder* getCCylinder() {return _cCylinder.get();}
    /// set chem cylinder
    /// @note: since this is a unique ptr, will implicitly delete old chem cylinder
    void setCCylinder(CCylinder* c) {_cCylinder = unique_ptr<CCylinder>(c);}

    /// Get cylinder type
    virtual int getType();
    int getType() const { return _type; }

    //@{
    /// Get beads
    Bead* getFirstBead() const {return _b1;}
    Bead* getSecondBead() const {return _b2;}
    //@}

    //@{
    /// Set beads
    void setFirstBead(Bead* b) {_b1 = b;}
    void setSecondBead(Bead* b) {_b2 = b;}
    //@}

    /// Get compartment
    Compartment* getCompartment() const { return _cellElement.manager->getHead(_cellElement); }

    //@{
    /// Branching cylinder management
    Cylinder* getBranchingCylinder() const {return _branchingCylinder;}
    void setBranchingCylinder(Cylinder* c) {_branchingCylinder = c;}
    //@}

    ///@{
    /// Set plus and minus end boolean markers
    bool isPlusEnd() {return _plusEnd;}
    void setPlusEnd(bool plusEnd) {_plusEnd = plusEnd;}

    bool isMinusEnd() {return _minusEnd;}
    void setMinusEnd(bool minusEnd) {_minusEnd = minusEnd;}
    //@}

    int getPosition() {return _position;}

    //@{
    /// SubSystem management, inherited from Trackable
    // Does nothing
    virtual void addToSubSystem() override {}
    virtual void removeFromSubSystem() override {}
    //@}

    /// Get all instances of this class from the SubSystem
    static const vector<Cylinder*>& getCylinders() {
        return getElements();
    }
    /// Get the number of cylinders in this system
    static int numCylinders() {
        return getElements().size();
    }

    /// Update the position, inherited from Movable
    /// @note - changes compartment if needed
    virtual void updatePosition();

    /// Update the reaction rates, inherited from Reactable
    virtual void updateReactionRates();

    /// Check if this cylinder is grown to full length
    bool isFullLength();

    virtual void printSelf()const;

                                       
    void setFilID(int filID){
       _filID = filID;
    };

    int getFilID(){
       return _filID;
    };

    /// Returns whether a cylinder is within a certain distance from another
    /// Uses the closest point between the two cylinders
    virtual bool within(Cylinder* other, floatingpoint dist);


    //Function to initialize Cylinders properly during restart
	void initializerestart(int nummonomers, int firstmonomer, int lastmonomer,
						bool minusendstatus, bool plusendstatus, short minusendtype,
						short plusendtype);

	//Adjusts the _alpha value for partial cylinders. Determines the ratio of
	// D(minusend-bindingsite)/D(minusend-plusend) given
	// _alpha = D(0thmonomer-bindingsite)/D(fullcylinder)
	//Note. Use this function only to determine mechanical coordinate.
	//Refer Docs/Design/PartialCylinderAlpha.pdf
    floatingpoint adjustedrelativeposition(floatingpoint _alpha, bool verbose = false);

    static floatingpoint timecylinder1;
	static floatingpoint timecylinder2;
	static floatingpoint timecylinderchem;
	static floatingpoint timecylindermech;
	static ofstream _crosscheckdumpFile;
    //@}
};

} // namespace medyan

#endif

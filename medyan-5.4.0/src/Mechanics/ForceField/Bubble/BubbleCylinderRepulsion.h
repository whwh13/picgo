
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

#ifndef MEDYAN_BubbleCylinderRepulsion_h
#define MEDYAN_BubbleCylinderRepulsion_h

#include <vector>

#include "common.h"

#include "Mechanics/ForceField/Bubble/BubbleCylinderRepulsionExp.h"
#include "Mechanics/ForceField/ForceField.h"
#include "NeighborListImpl.h"
#include "Structure/SubSystem.h"
#include "SysParams.h"

namespace medyan {
//FORWARD DECLARATIONS
class Bead;

/// Represents a repulsive interaction between a Bubble and Cylinder.
class BubbleCylinderRepulsion : public ForceField {
public:
    struct PairInteraction {
        Index bubbleCoordIndex = 0;
        Index beadCoordIndex = 0;
        floatingpoint krep = 0.0;
        floatingpoint slen = 0.0;
        floatingpoint radius = 0.0;
    };
private:
    BubbleBeadRepulsionExp _FFType;

    SubSystem* ps_ = nullptr;
    std::vector<PairInteraction> pairInteractions_;
public:

    /// Constructor
    BubbleCylinderRepulsion() = default;
    
    virtual void vectorize(const FFCoordinateStartingIndex&) override;
    
    virtual floatingpoint computeEnergy(floatingpoint *coord) override;
    virtual void computeForces(floatingpoint *coord, floatingpoint *f) override;
    
    virtual void computeLoadForces() override;
    virtual void computeLoadForce(SubSystem& sys, Cylinder* c, LoadForceEnd end) const override;
    virtual void whoIsCulprit() override {}
    
    /// Get the neighbor list for this interaction
    virtual std::vector<NeighborList*> getNeighborLists() override {return {};}
    
    virtual std::string getName() override {return "BubbleCylinderRepulsion";}
};

} // namespace medyan

#endif
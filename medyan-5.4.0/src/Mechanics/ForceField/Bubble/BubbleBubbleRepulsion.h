
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

#ifndef MEDYAN_BubbleBubbleRepulsion_h
#define MEDYAN_BubbleBubbleRepulsion_h

#include <vector>

#include "common.h"

#include "Mechanics/ForceField/Bubble/BubbleBubbleRepulsionExp.h"
#include "Mechanics/ForceField/ForceField.h"
#include "NeighborListImpl.h"
#include "Structure/SubSystem.h"
#include "SysParams.h"

namespace medyan {
//FORWARD DECLARATIONS
class Bead;

/// Represents a repulsive interaction between two [Bubbles](@ref Bubble).
class BubbleBubbleRepulsion : public ForceField {
public:
    struct PairInteraction {
        Index coordIndex1 = 0;
        Index coordIndex2 = 0;
        floatingpoint krep = 0.0;
        floatingpoint slen = 0.0;
        floatingpoint radius1 = 0.0;
        floatingpoint radius2 = 0.0;
    };
private:
    BubbleBubbleRepulsionExp _FFType;

    SubSystem* ps_ = nullptr;
    std::vector<PairInteraction> pairInteractions_;

public:
    
    /// Constructor
    BubbleBubbleRepulsion() = default;
    
    virtual void vectorize(const FFCoordinateStartingIndex&) override;
    
    virtual floatingpoint computeEnergy(floatingpoint *coord) override;
    virtual void computeForces(floatingpoint *coord, floatingpoint *f) override;
    
    virtual void computeLoadForces() override {}
    virtual void whoIsCulprit() override {}
    
    /// Get the neighbor list for this interaction
    virtual std::vector<NeighborList*> getNeighborLists() override {return {};}
    
    virtual std::string getName() override {return "BubbleBubbleRepulsion";}
};

} // namespace medyan

#endif

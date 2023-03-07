
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

#ifndef MEDYAN_MTOCBending_h
#define MEDYAN_MTOCBending_h

#include <vector>

#include "common.h"
#include "SysParams.h"
#include "Mechanics/ForceField/Bubble/MTOCBendingCosine.h"
#include "Mechanics/ForceField/ForceField.h"
#include "Structure/SubSystem.h"

namespace medyan {

/// Represents an attachment potential of a MTOC.
class MTOCBending : public ForceField {
public:
    struct Interaction {
        Index bubbleCoordIndex = 0;
        Index beadCoordIndex1 = 0;
        Index beadCoordIndex2 = 0;
        floatingpoint kbend = 0;
    };

private:
    MTOCBendingCosine _FFType;

    SubSystem* ps_ = nullptr;
    std::vector<Interaction> interactions_;
    
public:
        
    virtual void vectorize(const FFCoordinateStartingIndex&) override;
    
    virtual floatingpoint computeEnergy(floatingpoint *coord) override;
    virtual void computeForces(floatingpoint *coord, floatingpoint *f) override;
    
    virtual void computeLoadForces() override {}
    virtual void whoIsCulprit() override {}
    
    /// Get the neighbor list for this interaction
    virtual std::vector<NeighborList*> getNeighborLists() override {return {};}
    
    virtual std::string getName() override {return "MTOCAttachment";}
};

} // namespace medyan

#endif


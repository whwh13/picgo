
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

#ifndef MEDYAN_BoundaryBubbleRepulsion_h
#define MEDYAN_BoundaryBubbleRepulsion_h

#include <vector>

#include "common.h"

#include "Mechanics/ForceField/ForceField.h"
#include "NeighborListImpl.h"
#include "Structure/SubSystem.h"
#include "SysParams.h"

namespace medyan {
//FORWARD DECLARATIONS
class BoundaryElement;
class Bead;

/// Represents a repulsive interaction between a BoundaryElement and Bubble.
template <class BRepulsionInteractionType>
class BoundaryBubbleRepulsion : public ForceField {
    
private:
    BRepulsionInteractionType _FFType;

    SubSystem* ps_ = nullptr;
    std::vector<int> beadSet;

    ///Array describing the constants in calculation
    std::vector<FP> krep;
	std::vector<FP> slen;
    ///Array describing the number of neighbors for each boundary element (num boundary elements long)
    std::vector<int> nneighbors;
public:
    
    const static int n = 1;

    /// Constructor
    BoundaryBubbleRepulsion() = default;
    
    virtual FP computeEnergy(FP *coord) override;
   
    virtual void computeForces(FP *coord, FP *f) override;
    
    virtual std::string getName() override { return "BoundaryBubbleRepulsion"; }

    virtual void vectorize(const FFCoordinateStartingIndex& si, const SimulConfig&) override;

};

} // namespace medyan

#endif

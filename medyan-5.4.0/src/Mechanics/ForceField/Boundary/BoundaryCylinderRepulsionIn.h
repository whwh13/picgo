
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

#ifndef MEDYAN_BoundaryCylinderRepulsionIn_h
#define MEDYAN_BoundaryCylinderRepulsionIn_h

#include <vector>

#include "common.h"
#include "Mechanics/ForceField/ForceField.h"
#include "NeighborListImpl.h"

#include "SysParams.h"

namespace medyan {
//FORWARD DECLARATIONS
class BoundaryElement;
class Bead;
class Cylinder;

/// Represents a repulsive interaction between a BoundaryElement and Cylinder.
template <class BRepulsionInteractionType>
class BoundaryCylinderRepulsionIn : public ForceField {
    
private:
    BRepulsionInteractionType _FFType;
    std::unique_ptr<BoundaryCylinderNL> _neighborList; ///<Neighbor list of BoundaryElement - Cylinder
    
    std::vector<int> beadSet;
    
    ///Array describing the constants in calculation
    std::vector<FP> krep;
    std::vector<FP> slen;
    ///Array describing the number of neighbors for each boundary element (num boundary elements long)
    std::vector<int> nneighbors;
    
public:
    
    ///Array describing indexed set of interactions
    ///For filaments, this is a 1-bead potential
    const static int n = 1;
    
    /// Constructor
    BoundaryCylinderRepulsionIn(const SimulConfig& conf) {
        _neighborList = std::make_unique<BoundaryCylinderNL>(conf.boundParams.BoundaryCutoff);
    }
    
    virtual void vectorize(const FFCoordinateStartingIndex&, const SimulConfig&) override;
    
    virtual FP computeEnergy(FP *coord) override;
    //@{
    /// This repulsive force calculation also updates load forces
    /// on beads within the interaction range.
    virtual void computeForces(FP *coord, FP *f) override;
    
    virtual void computeLoadForces();
    virtual void computeLoadForce(SubSystem& sys, Cylinder* c, LoadForceEnd end) const override;
    //@}
    
    /// Get the neighbor list for this interaction
    virtual std::vector<NeighborList*> getNeighborLists() override {
        return { _neighborList.get() };
    }
    
    virtual std::string getName() override { return "BoundaryCylinderRepulsionIn"; }
};

} // namespace medyan

#endif

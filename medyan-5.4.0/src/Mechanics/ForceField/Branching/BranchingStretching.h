
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

#ifndef MEDYAN_BranchingStretching_h
#define MEDYAN_BranchingStretching_h

#include "common.h"
#ifdef CUDAACCL
#include "CUDAcommon.h"
#endif

#include "Mechanics/ForceField/ForceField.h"

namespace medyan {
//FORWARD DECLARATIONS
class BranchingPoint;

/// Represents an interaction fixing a Cylinder anchored by a BranchingPoint on the parent.
template <class BStretchingInteractionType>
class BranchingStretching : public ForceField {
    
private:
    BStretchingInteractionType _FFType;
    
    std::vector<int> beadSet;
    
    ///Array describing the constants in calculation
    std::vector<FP> kstr;
    std::vector<FP> eql;
    std::vector<FP> pos;
    std::vector<FP> stretchforce;

#ifdef CUDAACCL
    int * gpu_beadSet;
    floatingpoint * gpu_kstr;
    floatingpoint *gpu_eql;
    floatingpoint *gpu_pos;
    int * gpu_params;
    CUDAvars cvars;
    floatingpoint *F_i;
#endif
    
public:
    
    ///Array describing indexed set of interactions
    ///this is a 3-bead potential
    const static int n = 3;
    
    virtual void vectorize(const FFCoordinateStartingIndex&, const SimulConfig&) override;
    
    virtual FP computeEnergy(FP *coord) override;
    virtual void computeForces(FP *coord, FP *f) override;

    virtual void assignforcemags() override;
    
    virtual std::string getName() override { return "BranchingStretching"; }
};

} // namespace medyan

#endif

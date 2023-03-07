
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

#ifndef MEDYAN_BranchingPosition_h
#define MEDYAN_BranchingPosition_h

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
class BranchingPosition : public ForceField {
    
private:
    BStretchingInteractionType _FFType;
    
    std::vector<unsigned int> beadSet;
    
    ///Array describing the constants in calculation
    std::vector<FP> kpos;
    std::vector<FP> pos;
    std::vector<FP> stretchforce;
#ifdef CUDAACCL
    int * gpu_beadSet;
    floatingpoint * gpu_kpos;
    floatingpoint *gpu_pos;
    int * gpu_params;
    CUDAvars cvars;
    floatingpoint *F_i;
#endif
    
public:
    
    ///Array describing indexed set of interactions
    ///For filaments, this is a 3-bead potential
    const static int n = 3;
    
    virtual void vectorize(const FFCoordinateStartingIndex&, const SimulConfig&) override;
    
    virtual FP computeEnergy(FP *coord) override;
    virtual void computeForces(FP *coord, FP *f) override;

    virtual void assignforcemags() override;
    
    virtual std::string getName() override { return "BranchingPosition"; }
};

} // namespace medyan

#endif

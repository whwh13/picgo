
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

#ifndef MEDYAN_BranchingBending_h
#define MEDYAN_BranchingBending_h

#include "common.h"
#ifdef CUDAACCL
#include "CUDAcommon.h"
#endif
#include "Mechanics/ForceField/ForceField.h"

namespace medyan {

/// Represents an interaction maintaining a BranchingPoint angle (~70 for Arp2/3)
template <class BBendingInteractionType>
class BranchingBending : public ForceField {
    
private:
    BBendingInteractionType _FFType;
    
    std::vector<int> beadSet;
    
    ///Array describing the constants in calculation
    std::vector<FP> kbend;
    std::vector<FP> eqt;
    std::vector<FP> stretchforce;
#ifdef CUDAACCL
    int * gpu_beadSet;
    floatingpoint * gpu_kbend;
    floatingpoint *gpu_eqt;
    int * gpu_params;
    CUDAvars cvars;
    floatingpoint *F_i;
#endif
public:
    
    ///Array describing indexed set of interactions
    ///this is a 4-bead potential
    const static int n = 4;
    
    virtual void vectorize(const FFCoordinateStartingIndex&, const SimulConfig&) override;
    
    virtual FP computeEnergy(FP *coord) override;
    virtual void computeForces(FP *coord, FP *f) override;
    
    virtual std::string getName() override { return "BranchingBending"; }

    virtual void assignforcemags() override;
};

} // namespace medyan

#endif

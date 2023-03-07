
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

#ifndef MEDYAN_FilamentBending_h
#define MEDYAN_FilamentBending_h

#include "common.h"
#ifdef CUDAACCL
#include "CUDAcommon.h"
#endif
#include "Mechanics/ForceField/ForceField.h"

namespace medyan {
//FORWARD DECLARATIONS
class Filament;

/// Represents a Filament bending interaction
template <class FBendingInteractionType>
class FilamentBending : public ForceField {
    
private:
    FBendingInteractionType _FFType;
    
    // Cache of vectorized data
    Size _numInteractions;
    std::vector<int> beadSet;
    
    ///Array describing the constants in calculation
    std::vector<FP> kbend;
    std::vector<FP> eqt;

#ifdef CUDAACCL
    int * gpu_beadSet;
    floatingpoint *gpu_kbend;
    floatingpoint *gpu_eqt;
    int * gpu_params;
    CUDAvars cvars;
    floatingpoint *F_i;
    cudaStream_t stream = NULL;
#endif
public:
    
    ///Array describing indexed set of interactions
    ///For filaments, this is a 3-bead potential
    const static int n = 3;
    
    virtual void vectorize(const FFCoordinateStartingIndex&, const SimulConfig&) override;
    
    virtual FP computeEnergy(FP *coord) override;
    virtual void computeForces(FP *coord, FP *f) override;
    
    virtual std::string getName() override {return "FilamentBending";}
};

} // namespace medyan

#endif

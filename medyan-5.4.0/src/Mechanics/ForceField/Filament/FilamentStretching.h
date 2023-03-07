
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

#ifndef MEDYAN_FilamentStretching_h
#define MEDYAN_FilamentStretching_h
#include "Filament.h"
#include "Cylinder.h"
#include "common.h"
#include "Mechanics/ForceField/ForceField.h"
#ifdef CUDAACCL
#include "CUDAcommon.h"
#endif

namespace medyan {
/// Represents a Filament stretching interaction
template <class FStretchingInteractionType>
class FilamentStretching : public ForceField {
    
private:
    FStretchingInteractionType _FFType; 
    
    std::vector<int> beadSet;
    ///Array describing the constants in calculation
    std::vector<FP> kstr;
    std::vector<FP> eql;

#ifdef CUDAACCL
    int * gpu_beadSet;
    floatingpoint * gpu_kstr;
    floatingpoint *gpu_eql;
    int * gpu_params;
    CUDAvars cvars;
    floatingpoint *F_i;
    cudaStream_t stream = NULL;
#endif


public:
    
    ///Array describing indexed set of interactions
    ///For filaments, this is a 2-bead potential
    const static int n = 2;
    
    virtual void vectorize(const FFCoordinateStartingIndex&, const SimulConfig&) override;
    
    virtual FP computeEnergy(FP *coord) override;
    virtual void computeForces(FP *coord, FP *f) override;
    
    virtual std::string getName() override {return "FilamentStretching";}

};

} // namespace medyan

#endif


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

#ifndef MEDYAN_MotorGhostStretching_h
#define MEDYAN_MotorGhostStretching_h

#include "common.h"
#ifdef CUDAACCL
#include "CUDAcommon.h"
#endif

#include "Mechanics/ForceField/ForceField.h"

namespace medyan {

/// Represents a MotorGhost stretching interaction.
template <class MStretchingInteractionType>
class MotorGhostStretching : public ForceField {
    
private:
    MStretchingInteractionType _FFType;

    std::vector<int> beadSet;
    
    ///Array describing the constants in calculation
    std::vector<FP> kstr;
    std::vector<FP> eql;
    std::vector<FP> pos1;
    std::vector<FP> pos2;
    std::vector<FP> stretchforce;

#ifdef CUDAACCL
    int * gpu_beadSet;
    floatingpoint * gpu_kstr;
    floatingpoint *gpu_eql;
    int * gpu_params;
    floatingpoint *gpu_pos1;
    floatingpoint *gpu_pos2;
    floatingpoint *F_i;
    floatingpoint *gpu_Mstretchforce;
    cudaStream_t stream = NULL;
#endif
    
public:
    ///Array describing indexed set of interactions
    ///For linkers, this is a 4-bead potential
    const static int n = 4;
    
    ///< Constructor
    MotorGhostStretching () = default;
    
    virtual void vectorize(const FFCoordinateStartingIndex&, const SimulConfig&) override;
    
    
    virtual FP computeEnergy(FP *coord) override;
    virtual void computeForces(FP *coord, FP *f) override;


    virtual std::string getName() override { return "MotorStretching"; }

    virtual void assignforcemags() override;

};

} // namespace medyan

#endif

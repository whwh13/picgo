
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

#ifndef MEDYAN_BranchingDihedral_h
#define MEDYAN_BranchingDihedral_h

#include <vector>

#include "common.h"
#ifdef CUDAACCL
#include "CUDAcommon.h"
#endif
#include "Mechanics/ForceField/ForceField.h"

namespace medyan {
//FORWARD DECLARATIONS
class BranchingPoint;

/// Represents an interaction keeping BranchingPoint in dihedral plane
template <class BDihedralInteractionType>
class BranchingDihedral : public ForceField {
    
private:
    BDihedralInteractionType _FFType;

    std::vector< unsigned > beadSet;
    
    ///Array describing the constants in calculation
    std::vector<FP> kdih;
    std::vector<FP> pos;
    std::vector<FP> stretchforce;
#ifdef CUDAACCL
    int * gpu_beadSet;
    floatingpoint * gpu_kdih;
    floatingpoint *gpu_pos;
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

    virtual void assignforcemags() override;
    
    virtual std::string getName() override { return "BranchingDihedral"; }
};

} // namespace medyan

#endif

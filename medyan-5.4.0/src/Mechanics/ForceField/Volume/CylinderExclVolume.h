
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

#ifndef MEDYAN_CylinderExclVolume_h
#define MEDYAN_CylinderExclVolume_h

#include "common.h"

#include "Mechanics/ForceField/ForceField.h"
#include "Mechanics/ForceField/Volume/CylinderExclVolRepulsion.h"
#include "NeighborListImpl.h"
#include "HybridNeighborListImpl.h"

#include "SysParams.h"
#ifdef CUDAACCL
#include "CUDAcommon.h"
#endif

namespace medyan {
//FORWARD DECLARATIONS
class Cylinder;

/// Represents an excuded volume interaction between two [Cylinders](@ref Cylinder).
template <class CVolumeInteractionType>
class CylinderExclVolume : public ForceField {

private:
    CVolumeInteractionType _FFType;
    CylinderCylinderNL* _neighborList = nullptr;  ///< Neighbor list of cylinders
#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
    HybridCylinderCylinderNL* _HneighborList;
    vector<short> _HnlIDvec;
#else
    short _HnlID;
#endif
    ///Array describing the constants in calculation
    std::vector<int> beadSet;
    std::vector<FP> krep;
    std::vector<floatingpoint> vecEqLength;

#ifdef CUDAACCL
    int * gpu_beadSet = NULL;
    floatingpoint * gpu_krep = NULL;
    int * gpu_params = NULL;
    CUDAvars cvars;
    floatingpoint *F_i;
    cudaStream_t stream = NULL;
#endif
public:
    ///Array describing indexed set of interactions
    ///For volume, this is a 4-bead potential
    const static int n = 4;
    int numInteractions = 0;

    ///Constructor
    CylinderExclVolume(const SimulConfig& conf) {
        //If Hybrid NeighborList is not preferred, neighborList is created using Original
        // framework.
#if !defined(HYBRID_NLSTENCILLIST) || !defined(SIMDBINDINGSEARCH)
        //Not a full list as it does not pass a full variable.
        _neighborList = new CylinderCylinderNL(conf.mechParams.VolumeCutoff);
#endif
#ifdef CUDAACCL_NL
        _neighborList->cudacpyforces = true;
#endif
    }

    virtual std::string getName() override { return "CylinderExcludeVolume"; }

    virtual void vectorize(const FFCoordinateStartingIndex&, const SimulConfig&) override;
    
    virtual FP computeEnergy(FP *coord) override;
    //@{
    /// This repulsive force calculation also updates load forces
    /// on beads within the interaction range.
    virtual void computeForces(FP *coord, FP *f) override;

    /// Get the neighbor list for this interaction
    virtual std::vector<NeighborList*> getNeighborLists() override {
        return { _neighborList };
    }


    void setHNeighborLists(HybridCylinderCylinderNL* Hnl) {
        _HneighborList = Hnl;
#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
        short filtypes = SysParams::Geometry().cylinderNumMon.size();
        for (int i = 0; i < filtypes; i++){
            for (int j = i; j < filtypes; j++){
                _HnlIDvec.push_back(Hnl->setneighborsearchparameters(i,j,true,false,SysParams::Mechanics()
                                                                     .VolumeCutoff,0.0));
            }
            
        }
#else
        _HnlID = Hnl->setneighborsearchparameters(0,0,true,false,SysParams::Mechanics()
                                                                .VolumeCutoff,0.0);
#endif
    };

    vector<tuple<floatingpoint, int, vector<tuple<floatingpoint*,floatingpoint*,floatingpoint*,floatingpoint*, floatingpoint>>>> getCylEnergies() {
        return _FFType.getCylEnergies();
    };
    
    void clearCylEnergies(){
        _FFType.clearCylEnergies();
    }
};

} // namespace medyan

#endif

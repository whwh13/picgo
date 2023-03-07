
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

#include "BranchingDihedral.h"

#include "Mechanics/ForceField/Branching/BranchingDihedralCosine.h"
#include "Mechanics/ForceField/Branching/BranchingDihedralQuadratic.hpp"
#include "Mechanics/ForceField/Branching/BranchingDihedralCosineV2.h"
#include "Mechanics/ForceField/Branching/BranchingDihedralQuadraticV2.h"

#include "BranchingPoint.h"
#include "Cylinder.h"
#include "Bead.h"
#include "cross_check.h"
#ifdef CUDAACCL
#include "nvToolsExt.h"
#endif
#include "Mechanics/CUDAcommon.h"

namespace medyan {
template <class BDihedralInteractionType>
void BranchingDihedral<BDihedralInteractionType>::vectorize(const FFCoordinateStartingIndex& si, const SimulConfig& conf) {

    CUDAcommon::tmin.numinteractions[6] += BranchingPoint::getBranchingPoints().size();
    beadSet.resize(n * BranchingPoint::getBranchingPoints().size());
    kdih.assign(BranchingPoint::getBranchingPoints().size(), 0);
    pos.assign(BranchingPoint::getBranchingPoints().size(), 0);
    stretchforce.assign(3*BranchingPoint::getBranchingPoints().size(), 0);

    int i = 0;

    for (auto b: BranchingPoint::getBranchingPoints()) {

        beadSet[n * i] = b->getFirstCylinder()->getFirstBead()->getIndex() * 3 + si.bead;
        beadSet[n * i + 1] = b->getFirstCylinder()->getSecondBead()->getIndex() * 3 + si.bead;
        beadSet[n * i + 2] = b->getSecondCylinder()->getFirstBead()->getIndex() * 3 + si.bead;
        beadSet[n * i + 3] = b->getSecondCylinder()->getSecondBead()->getIndex() * 3 + si.bead;

        kdih[i] = b->getMBranchingPoint()->getDihedralConstant();
        pos[i] = b->getFirstCylinder()->adjustedrelativeposition(b->getPosition());
        for(int j = 0; j < 3; j++)
            stretchforce[3*i + j] = 0.0;

        // Reset branch force as a side effect.
        b->getMBranchingPoint()->branchForce = { 0.0, 0.0, 0.0 };

        i++;
    }
    //CUDA
#ifdef CUDAACCL
//    cudaEvent_t start, stop;
//    CUDAcommon::handleerror(cudaEventCreate( &start));
//    CUDAcommon::handleerror(cudaEventCreate( &stop));
//    CUDAcommon::handleerror(cudaEventRecord( start, 0));

    int numInteractions =BranchingPoint::getBranchingPoints().size();
    _FFType.optimalblocksnthreads(numInteractions);
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_beadSet, n * numInteractions * sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_beadSet, beadSet.data(), n * numInteractions * sizeof(int),
                                       cudaMemcpyHostToDevice));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_kdih, numInteractions * sizeof(floatingpoint)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_kdih, kdih, numInteractions * sizeof(floatingpoint), cudaMemcpyHostToDevice));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_pos, numInteractions * sizeof(floatingpoint)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_pos, pos, numInteractions * sizeof(floatingpoint), cudaMemcpyHostToDevice));

    vector<int> params;
    params.push_back(int(n));
    params.push_back(numInteractions);

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, 2 * sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpy(gpu_params, params.data(), 2 * sizeof(int), cudaMemcpyHostToDevice));

#endif
}

template<class BDihedralInteractionType>
void BranchingDihedral<BDihedralInteractionType>::assignforcemags() {
    for(auto b:BranchingPoint::getBranchingPoints()){
        //Using += to ensure that the stretching forces are additive.

        for(int j = 0; j < 3; j++)
            b->getMBranchingPoint()->branchForce[j] += stretchforce[3*b->getIndex() + j];
    }
}


template <class BDihedralInteractionType>
FP BranchingDihedral<BDihedralInteractionType>::computeEnergy(FP *coord) {

#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    floatingpoint * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    floatingpoint * gpu_force=CUDAcommon::getCUDAvars().gpu_force;
    floatingpoint * gpu_d = CUDAcommon::getCUDAvars().gpu_lambda;


//    if(d == 0.0){
//        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kdih, gpu_pos, gpu_params);
//
//    }
//    else{
        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kdih, gpu_pos, gpu_d,
                            gpu_params);
//    }

#endif

    FP U_ii = _FFType.energy(coord, BranchingPoint::getBranchingPoints().size(), beadSet.data(), kdih.data(), pos.data());

#if defined(SERIAL_CUDACROSSCHECK) && defined(DETAILEDOUTPUT_ENERGY)
    floatingpoint U_i[1];
    floatingpoint* gU_i;

    CUDAcommon::handleerror(cudaDeviceSynchronize(),"ForceField", "ForceField");
    floatingpoint cuda_energy[1];
    if(gU_i == NULL)
        cuda_energy[0] = 0.0;
    else {
        CUDAcommon::handleerror(cudaMemcpy(cuda_energy, gU_i, sizeof(floatingpoint),
                                           cudaMemcpyDeviceToHost));
    }
    std::cout<<getName()<<" Serial Energy "<<U_ii<<" Cuda Energy "<<cuda_energy[0]<<endl;
#endif
    return U_ii;

}

template <class BDihedralInteractionType>
void BranchingDihedral<BDihedralInteractionType>::computeForces(FP *coord, FP *f) {
#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    floatingpoint * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;

    floatingpoint * gpu_force;

    if(cross_checkclass::Aux){

        gpu_force=CUDAcommon::getCUDAvars().gpu_forceAux;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kdih, gpu_pos, gpu_params);

    }
    else {


        gpu_force = CUDAcommon::getCUDAvars().gpu_force;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kdih, gpu_pos, gpu_params);

    }
#endif
    _FFType.forces(coord, f, BranchingPoint::getBranchingPoints().size(), beadSet.data(),
            kdih.data(), pos.data(), stretchforce.data());

}

// Template instantiations.
template class BranchingDihedral< BranchingDihedralCosine >;
template class BranchingDihedral< BranchingDihedralCosineV2 >;
template class BranchingDihedral< BranchingDihedralQuadratic >;
template class BranchingDihedral< BranchingDihedralQuadraticV2 >;

} // namespace medyan

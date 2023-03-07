
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

#include "MotorGhostStretching.h"

#include "MotorGhostStretchingHarmonic.h"

#include "MotorGhost.h"
#include "Cylinder.h"
#include "Bead.h"
#include "cross_check.h"
#ifdef CUDAACCL
#include "nvToolsExt.h"
#endif
#include "Mechanics/CUDAcommon.h"

namespace medyan {

template <class MStretchingInteractionType>
void MotorGhostStretching<MStretchingInteractionType>::vectorize(const FFCoordinateStartingIndex& si, const SimulConfig& conf) {

    CUDAcommon::tmin.numinteractions[3] += MotorGhost::getMotorGhosts().size();
    beadSet.assign(n * MotorGhost::getMotorGhosts().size(), 0);
    kstr.assign(MotorGhost::getMotorGhosts().size(), 0);
    eql.assign(MotorGhost::getMotorGhosts().size(), 0);
    pos1.assign(MotorGhost::getMotorGhosts().size(), 0);
    pos2.assign(MotorGhost::getMotorGhosts().size(), 0);
    stretchforce.assign(MotorGhost::getMotorGhosts().size(), 0);

    int i = 0;

    for (auto m: MotorGhost::getMotorGhosts()) {
        beadSet[n * i] = m->getFirstCylinder()->getFirstBead()->getIndex() * 3 + si.bead;
        beadSet[n * i + 1] = m->getFirstCylinder()->getSecondBead()->getIndex() * 3 + si.bead;
        beadSet[n * i + 2] = m->getSecondCylinder()->getFirstBead()->getIndex() * 3 + si.bead;
        beadSet[n * i + 3] = m->getSecondCylinder()->getSecondBead()->getIndex() * 3 + si.bead;

        kstr[i] = m->getMMotorGhost()->getStretchingConstant();
        eql[i] = m->getMMotorGhost()->getEqLength();
        pos1[i] = m->getFirstCylinder()->adjustedrelativeposition(m->getFirstPosition());
        pos2[i] = m->getSecondCylinder()->adjustedrelativeposition(m->getSecondPosition());
        stretchforce[i] = 0.0;

        // Reset motor stretch force as a side effect.
        m->getMMotorGhost()->stretchForce = 0;

        i++;
    }

    //CUDA
#ifdef CUDAACCL
    //CUDA stream create
    if(stream == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamCreate(&stream));
//    F_i = new floatingpoint[3 * Bead::getBeads().size()];
//    cudaEvent_t start, stop;
//    CUDAcommon::handleerror(cudaEventCreate( &start));
//    CUDAcommon::handleerror(cudaEventCreate( &stop));
//    CUDAcommon::handleerror(cudaEventRecord( start, 0));

    int numInteractions = MotorGhost::getMotorGhosts().size();
    _FFType.optimalblocksnthreads(numInteractions, stream);
//    blocksnthreads.clear();
//    blocksnthreads.push_back(numInteractions/THREADSPERBLOCK + 1);
//
//    if(blocksnthreads[0]==1) blocksnthreads.push_back( numInteractions);
////    if(blocksnthreads[0]==1) blocksnthreads.push_back( 32*(int(numInteractions/32 +1)) );
//    else blocksnthreads.push_back(THREADSPERBLOCK);
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_beadSet, n * numInteractions * sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_beadSet, beadSet, n * numInteractions *
                                                sizeof(int),
                                       cudaMemcpyHostToDevice, stream));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_kstr, numInteractions * sizeof(floatingpoint)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_kstr, kstr, numInteractions * sizeof
                            (floatingpoint), cudaMemcpyHostToDevice, stream));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_eql, numInteractions * sizeof(floatingpoint)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_eql, eql, numInteractions * sizeof(floatingpoint),
                                        cudaMemcpyHostToDevice, stream));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_pos1, numInteractions * sizeof(floatingpoint)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_pos1, pos1, numInteractions * sizeof
                            (floatingpoint), cudaMemcpyHostToDevice, stream));

//    floatingpoint checkpos1[numInteractions];
//    cudaMemcpy(checkpos1, gpu_pos1, numInteractions * sizeof(floatingpoint), cudaMemcpyDeviceToHost);
//    for(auto i=0;i<numInteractions;i++) std::cout<<pos1[i]<<" "<<checkpos1[i]<<endl;

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_pos2, numInteractions * sizeof(floatingpoint)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_pos2, pos2, numInteractions * sizeof
                           (floatingpoint), cudaMemcpyHostToDevice, stream));
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_Mstretchforce, numInteractions *
                                                                     sizeof(floatingpoint)),"cuda data transfer",
                            "MotorGhostStretching.cu");

    vector<int> params;
    params.push_back(int(n));
    params.push_back(numInteractions);
    params.push_back(CUDAcommon::cudavars.offset_E);
    //set offset
    CUDAcommon::cudavars.offset_E += numInteractions;
//    std::cout<<"offset "<<getName()<<" "<<CUDAcommon::cudavars.offset_E<<endl;

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, 3 * sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_params, params.data(), 3 * sizeof(int),
                                       cudaMemcpyHostToDevice, stream));
//    CUDAcommon::cudavars.motorparams = gpu_params;
#endif

    //
}

template<class MStretchingInteractionType>
void MotorGhostStretching<MStretchingInteractionType>::assignforcemags() {
    for(auto m: MotorGhost::getMotorGhosts()){
        //Using += to ensure that the stretching forces are additive.
        m->getMMotorGhost()->stretchForce += stretchforce[m->getIndex()];
    }
}


template <class MStretchingInteractionType>
FP MotorGhostStretching<MStretchingInteractionType>::computeEnergy(FP* coord){
#ifdef CUDAACCL

    //has to be changed to accomodate aux force
    floatingpoint * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    floatingpoint * gpu_force=CUDAcommon::getCUDAvars().gpu_force;
    floatingpoint * gpu_d = CUDAcommon::getCUDAvars().gpu_lambda;


//    if(d == 0.0){
//        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos1, gpu_pos2, gpu_params);
//
//    }
//    else{
        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos1, gpu_pos2, gpu_d,
                            gpu_params);
//    }



#endif

    FP U_ii = _FFType.energy(coord, beadSet.data(), kstr.data(), eql.data(), pos1.data(), pos2.data());


    return U_ii;
}

template <class MStretchingInteractionType>
void MotorGhostStretching<MStretchingInteractionType>::computeForces(FP *coord, FP *f) {
#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    floatingpoint * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;

    floatingpoint * gpu_force;
    if(cross_checkclass::Aux){
        gpu_force=CUDAcommon::getCUDAvars().gpu_forceAux;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos1,
                       gpu_pos2, gpu_params, gpu_Mstretchforce);
    }
    else {
        gpu_force = CUDAcommon::getCUDAvars().gpu_force;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos1,
                       gpu_pos2, gpu_params, gpu_Mstretchforce);
    }
#endif
    _FFType.forces(coord, f, beadSet.data(), kstr.data(), eql.data(), pos1.data(), pos2.data(), stretchforce.data());
#ifdef DETAILEDOUTPUT
    floatingpoint maxF = 0.0;
    floatingpoint mag = 0.0;
    for(int i = 0; i < CGMethod::N/3; i++) {
        mag = 0.0;
        for(int j = 0; j < 3; j++)
            mag += f[3 * i + j]*f[3 * i + j];
        mag = sqrt(mag);
//        std::cout<<"SL "<<i<<" "<<mag*mag<<" "<<forceAux[3 * i]<<" "<<forceAux[3 * i + 1]<<" "<<forceAux[3 * i +
//                2]<<endl;
        if(mag > maxF) maxF = mag;
    }
    std::cout<<"max "<<getName()<<" "<<maxF<<endl;
#endif
}

// Explicit instantiation.
template class MotorGhostStretching<MotorGhostStretchingHarmonic>;

} // namespace medyan

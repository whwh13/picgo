
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

#include "LinkerStretching.h"

#include "LinkerStretchingHarmonic.h"

#include "Cylinder.h"
#include "Linker.h"
#include "Bead.h"
#include "cross_check.h"
#ifdef CUDAACCL
#include "nvToolsExt.h"
#endif
#include "Mechanics/CUDAcommon.h"

namespace medyan {

template <class LStretchingInteractionType>
void LinkerStretching<LStretchingInteractionType>::vectorize(const FFCoordinateStartingIndex& si, const SimulConfig& conf) {
    CUDAcommon::tmin.numinteractions[2] += Linker::getLinkers().size();
    beadSet.assign(n * Linker::getLinkers().size(), 0);//stableindex of the bead
    kstr.assign(Linker::getLinkers().size(), 0);
    eql.assign(Linker::getLinkers().size(), 0);
    pos1.assign(Linker::getLinkers().size(), 0);
    pos2.assign(Linker::getLinkers().size(), 0);
    stretchforce.assign(Linker::getLinkers().size(), 0);

    int i = 0;
    //Filling stage
    for (auto l: Linker::getLinkers()) {
        beadSet[n * i] = l->getFirstCylinder()->getFirstBead()->getIndex() * 3 + si.bead;
        beadSet[n * i + 1] = l->getFirstCylinder()->getSecondBead()->getIndex() * 3 + si.bead;
        beadSet[n * i + 2] = l->getSecondCylinder()->getFirstBead()->getIndex() * 3 + si.bead;
        beadSet[n * i + 3] = l->getSecondCylinder()->getSecondBead()->getIndex() * 3 + si.bead;

        kstr[i] = l->getMLinker()->getStretchingConstant();
        eql[i] = l->getMLinker()->getEqLength();
        pos1[i] = l->getFirstCylinder()->adjustedrelativeposition(l->getFirstPosition());
        pos2[i] = l->getSecondCylinder()->adjustedrelativeposition(l->getSecondPosition());
        stretchforce[i] = 0.0;

        // Reset linker stretch force as a side effect.
        l->getMLinker()->stretchForce = 0;

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

    int numInteractions =Linker::getLinkers().size();
    _FFType.optimalblocksnthreads(numInteractions, stream);
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_beadSet, n * numInteractions * sizeof(int)),"cuda data transfer",
                                       "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_beadSet, beadSet, n * numInteractions *
                                                sizeof(int),
                                       cudaMemcpyHostToDevice, stream),"cuda data transfer", "LinkerStretching.cu");

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_kstr, numInteractions * sizeof(floatingpoint)),"cuda data transfer",
                                       "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_kstr, kstr, numInteractions * sizeof
                            (floatingpoint), cudaMemcpyHostToDevice, stream),
                            "cuda data transfer", "LinkerStretching.cu");

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_eql, numInteractions * sizeof(floatingpoint)),"cuda data transfer",
                            "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_eql, eql, numInteractions * sizeof(floatingpoint),
                            cudaMemcpyHostToDevice, stream), "cuda data transfer", "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_pos1, numInteractions * sizeof(floatingpoint)),"cuda data transfer",
                            "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_pos1, pos1, numInteractions * sizeof
                            (floatingpoint), cudaMemcpyHostToDevice, stream),
                            "cuda data transfer", "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_pos2, numInteractions * sizeof(floatingpoint)),"cuda data transfer",
                            "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_pos2, pos2, numInteractions * sizeof
                            (floatingpoint), cudaMemcpyHostToDevice, stream),
                            "cuda data transfer", "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_Lstretchforce, numInteractions *
                                                                     sizeof(floatingpoint)),"cuda data transfer",
                            "LinkerStretching.cu");
    vector<int> params;
    params.push_back(int(n));
    params.push_back(numInteractions);
    params.push_back(CUDAcommon::cudavars.offset_E);
    //set offset
    CUDAcommon::cudavars.offset_E += numInteractions;
//    std::cout<<"offset "<<getName()<<" "<<CUDAcommon::cudavars.offset_E<<endl;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, 3 * sizeof(int)),"cuda data"
                                    " transfer",
                            "LinkerStretching.cu");
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_params, params.data(), 3 * sizeof(int),
                                       cudaMemcpyHostToDevice, stream),
                            "cuda data transfer", "LinkerStretching.cu");
#endif
}

template<class LStretchingInteractionType>
void LinkerStretching<LStretchingInteractionType>::assignforcemags() {
    for(auto l:Linker::getLinkers()){
        //Using += to ensure that the stretching forces are additive.
        l->getMLinker()->stretchForce += stretchforce[l->getIndex()];
    }
}


template <class LStretchingInteractionType>
FP LinkerStretching<LStretchingInteractionType>::computeEnergy(FP* coord){

#ifdef CUDAACCL
//    std::cout<<"Linker size "<<Linker::getLinkers().size()<<endl;
	floatingpoint* gU_i;
    //has to be changed to accomodate aux force
    floatingpoint * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    floatingpoint * gpu_force=CUDAcommon::getCUDAvars().gpu_force;
    floatingpoint * gpu_d = CUDAcommon::getCUDAvars().gpu_lambda;


//    if(d == 0.0){
//        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos1, gpu_pos2, gpu_params);
//    }
//    else{
        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos1, gpu_pos2, gpu_d,
                            gpu_params);
//    }
#endif
    FP U_ii = _FFType.energy(coord, beadSet.data(), kstr.data(), eql.data(), pos1.data(), pos2.data());


    return U_ii;
}

template <class LStretchingInteractionType>
void LinkerStretching<LStretchingInteractionType>::computeForces(FP *coord, FP *f) {
#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    floatingpoint * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    floatingpoint * gpu_force;
    if(cross_checkclass::Aux){
        gpu_force=CUDAcommon::getCUDAvars().gpu_forceAux;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos1,
                       gpu_pos2, gpu_params, gpu_Lstretchforce);
    }
    else {
        gpu_force = CUDAcommon::getCUDAvars().gpu_force;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_pos1,
                       gpu_pos2, gpu_params , gpu_Lstretchforce);
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
template class LinkerStretching<LinkerStretchingHarmonic>;

} // namespace medyan

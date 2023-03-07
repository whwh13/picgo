
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

#include "FilamentStretching.h"

#include "FilamentStretchingHarmonic.h"
#include "Bead.h"
#include "cross_check.h"
#include "Mechanics/CUDAcommon.h"
#ifdef CUDAACCL
#include "nvToolsExt.h"
#endif

namespace medyan {
template <class FStretchingInteractionType>
void FilamentStretching<FStretchingInteractionType>::vectorize(const FFCoordinateStartingIndex& si, const SimulConfig& conf) {
    CUDAcommon::tmin.numinteractions[0] += Cylinder::getCylinders().size();
    beadSet.assign(n * Cylinder::getCylinders().size(), 0);
    kstr.assign(Cylinder::getCylinders().size(), 0);
    eql.assign(Cylinder::getCylinders().size(), 0);

    int i = 0;

    for (auto c: Cylinder::getCylinders()) {
        beadSet[n * i] = c->getFirstBead()->getIndex() * 3 + si.bead;
        beadSet[n * i + 1] = c->getSecondBead()->getIndex() * 3 + si.bead;
        kstr[i] = c->getMCylinder()->getStretchingConst();
        eql[i] = c->getMCylinder()->getEqLength();
/*        std::cout<<"Filstretching with cindex "<<c->_dcIndex<<" and ID "
                ""<<c->getID()<<" with bindices "<<c->getFirstBead()
                         ->getIndex()<<" "<<c->getSecondBead()->getIndex()<<endl;*/
        i++;
    }
    //CUDA
#ifdef CUDAACCL
    int numInteractions = Cylinder::getCylinders().size();
    _FFType.optimalblocksnthreads(numInteractions, stream);

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_beadSet, n * numInteractions * sizeof(int)),"cuda data "
            "transfer", " FilamentStretching.cu");
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_beadSet, beadSet, n * numInteractions *
                            sizeof(int), cudaMemcpyHostToDevice, stream),"cuda data transfer", " FilamentStretching.cu");

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_kstr, numInteractions * sizeof(floatingpoint)),"cuda data transfer",
                            " FilamentStretching.cu");
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_kstr, kstr, numInteractions * sizeof
                            (floatingpoint), cudaMemcpyHostToDevice, stream), "cuda data transfer", " FilamentStretching.cu");

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_eql, numInteractions * sizeof(floatingpoint)),"cuda data transfer",
                            " FilamentStretching.cu");
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_eql, eql, numInteractions * sizeof(floatingpoint),
                            cudaMemcpyHostToDevice, stream), "cuda data transfer", " FilamentStretching.cu");

    vector<int> params;
    params.push_back(int(n));
    params.push_back(numInteractions);
    params.push_back(CUDAcommon::cudavars.offset_E);
    //set offset
    CUDAcommon::cudavars.offset_E += numInteractions;
//    std::cout<<"offset "<<getName()<<" "<<CUDAcommon::cudavars.offset_E<<endl;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, 3 * sizeof(int)),"cuda data"
                                    " transfer", " FilamentStretching.cu");
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_params, params.data(), 3 * sizeof(int),
                                       cudaMemcpyHostToDevice, stream),
            "cuda data transfer", " FilamentStretching.cu");
#endif

}


template <class FStretchingInteractionType>
FP FilamentStretching<FStretchingInteractionType>::computeEnergy(FP* coord){

#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    floatingpoint * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    floatingpoint * gpu_force=CUDAcommon::getCUDAvars().gpu_force;
    floatingpoint * gpu_d = CUDAcommon::getCUDAvars().gpu_lambda;
//    std::cout<<"Fil Stretching Forces"<<endl;


//    if(d == 0.0){
//        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_params);
//
//    }
//    else{
        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_d,
                            gpu_params);
//    }
#endif

    FP U_ii = _FFType.energy(coord, beadSet.data(), kstr.data(), eql.data());

    return U_ii;
}

template <class FStretchingInteractionType>
void FilamentStretching<FStretchingInteractionType>::computeForces(FP *coord, FP *f) {
#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    floatingpoint * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;

    floatingpoint * gpu_force;
    if(cross_checkclass::Aux){
        gpu_force=CUDAcommon::getCUDAvars().gpu_forceAux;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_params);
    }
    else {
        gpu_force = CUDAcommon::getCUDAvars().gpu_force;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kstr, gpu_eql, gpu_params);
    }
#endif
    _FFType.forces(coord, f, beadSet.data(), kstr.data(), eql.data());
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
    std::cout<<"max "<<getName()<<" "<<maxF<<" nint "<<Cylinder::getCylinders().size()
             <<endl;
#endif
}

///Temlate specializations
template class FilamentStretching<FilamentStretchingHarmonic>;

} // namespace medyan

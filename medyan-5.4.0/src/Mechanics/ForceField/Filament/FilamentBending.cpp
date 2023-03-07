
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
#include "FilamentBending.h"

#include "FilamentBendingHarmonic.h"
#include "FilamentBendingCosine.h"

#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"
#include "SysParams.h"
#ifdef CUDAACCL
#include "nvToolsExt.h"
#endif
#include "cross_check.h"

namespace medyan {

template <class FBendingInteractionType>
void FilamentBending<FBendingInteractionType>::vectorize(const FFCoordinateStartingIndex& si, const SimulConfig& conf) {

    // Count number of interactions
    _numInteractions = 0;
    for(auto f : Filament::getFilaments())
        if(f->getCylinderVector().size() > 1) _numInteractions += f->getCylinderVector().size() - 1;

    beadSet.assign(n * _numInteractions, 0);
    kbend.assign(_numInteractions, 0);
    eqt.assign(_numInteractions, 0);

    int i = 0;

    int istr = 0;

    for (auto f: Filament::getFilaments()) {

        if (f->getCylinderVector().size() > 1){

            for (auto it = f->getCylinderVector().begin()+1;
                 it != f->getCylinderVector().end(); it++){

                auto it2 = it - 1;
                beadSet[n * i] = (*it2)->getFirstBead()->getIndex() * 3 + si.bead;
                beadSet[n * i + 1] = (*it)->getFirstBead()->getIndex() * 3 + si.bead;
                beadSet[n * i + 2] = (*it)->getSecondBead()->getIndex() * 3 + si.bead;
//                std::cout<<f->getCylinderVector().size()<<" "<<(*it2)->getFirstBead()
//                        ->getIndex()<<" "<<(*it)->getFirstBead()
//                        ->getIndex()<<" "<<(*it)->getSecondBead()->getIndex()<<endl;
                kbend[i] = (*it)->getMCylinder()->getBendingConst();
                eqt[i]  = (*it)->getMCylinder()->getEqTheta();

                i++;
            }
        }
    }

    //CUDA
#ifdef CUDAACCL
    //CUDA stream create
    if(stream == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamCreate(&stream));

//    F_i = new floatingpoint[3 * Bead::getBeads().size()];

    _FFType.optimalblocksnthreads(_numInteractions, stream);

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_beadSet, n * _numInteractions * sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_beadSet, beadSet, n * _numInteractions *
                                                sizeof(int),
                                       cudaMemcpyHostToDevice, stream));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_kbend, _numInteractions * sizeof(floatingpoint)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_kbend, kbend, _numInteractions * sizeof
                            (floatingpoint), cudaMemcpyHostToDevice, stream));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_eqt, _numInteractions * sizeof(floatingpoint)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_eqt, eqt, _numInteractions * sizeof(floatingpoint),
                                        cudaMemcpyHostToDevice, stream));

    vector<int> params;
    params.push_back(int(n));
    params.push_back(_numInteractions);
    params.push_back(CUDAcommon::cudavars.offset_E);
    //set offset
    CUDAcommon::cudavars.offset_E += _numInteractions;
//    std::cout<<"offset "<<getName()<<" "<<CUDAcommon::cudavars.offset_E<<endl;
//    std::cout<<"offset "<<getName()<<" "<<CUDAcommon::cudavars.offset_E<<endl;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, 3 * sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_params, params.data(), 3 * sizeof(int),
                                       cudaMemcpyHostToDevice, stream));

#endif
}


template <class FBendingInteractionType>
FP FilamentBending<FBendingInteractionType>::computeEnergy(FP *coord){

#ifdef CUDAACCL
    floatingpoint* gU_i;
    //has to be changed to accomodate aux force
    floatingpoint * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    floatingpoint * gpu_force=CUDAcommon::getCUDAvars().gpu_force;
    floatingpoint * gpu_d = CUDAcommon::getCUDAvars().gpu_lambda;


//    if(d == 0.0){
//        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kbend, gpu_eqt, gpu_params);
//
//    }
//    else{
        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_kbend, gpu_eqt, gpu_d,
                            gpu_params);
//    }
#endif

    FP U_ii = _FFType.energy(coord, _numInteractions, beadSet.data(), kbend.data(), eqt.data());

    return U_ii;

}

template <class FBendingInteractionType>
void FilamentBending<FBendingInteractionType>::computeForces(FP *coord, FP *f) {
#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    floatingpoint * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    floatingpoint * gpu_force;
    if(cross_checkclass::Aux){
        gpu_force=CUDAcommon::getCUDAvars().gpu_forceAux;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kbend, gpu_eqt, gpu_params);
    }
    else {
        gpu_force = CUDAcommon::getCUDAvars().gpu_force;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_kbend, gpu_eqt, gpu_params);
    }
#endif
    _FFType.forces(coord, f, _numInteractions, beadSet.data(), kbend.data(), eqt.data());
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
template class FilamentBending<FilamentBendingHarmonic>;
template class FilamentBending<FilamentBendingCosine>;

} // namespace medyan

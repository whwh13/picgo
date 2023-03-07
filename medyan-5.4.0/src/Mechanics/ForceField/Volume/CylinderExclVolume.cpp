
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

#include "CylinderExclVolume.h"

#include "CylinderExclVolRepulsion.h"

#include "Cylinder.h"
#include "Bead.h"

#include "MathFunctions.h"
#include "cross_check.h"
#ifdef CUDAACCL
#include "nvToolsExt.h"
#endif
#include "Mechanics/CUDAcommon.h"
#include "Structure/DofSerializer.hpp"


namespace medyan {
using namespace mathfunc;

template <class CVolumeInteractionType>
void CylinderExclVolume<CVolumeInteractionType>::vectorize(const FFCoordinateStartingIndex& si, const SimulConfig& conf) {
    //count interactions
    int nint = 0;

    for(auto ci : Cylinder::getCylinders()) {

#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
        for (int ID = 0; ID < _HnlIDvec.size(); ID ++){
            auto neighbors = _HneighborList->getNeighborsstencil(_HnlIDvec[ID], ci);
            for(auto &cn : neighbors)
            {
                if(cn->getBranchingCylinder() == ci) continue;
                nint++;
            }
        }
#else
        auto neighbors = _neighborList->getNeighbors(ci);

        for(auto &cn : neighbors)
        {

            if(cn->getBranchingCylinder() == ci) continue;

            nint++;
        }
#endif
    }

    numInteractions = nint;
    CUDAcommon::tmin.numinteractions[8] += numInteractions;
//    std::cout<<"NINT1 "<<nint<<endl;
    beadSet.assign(n * nint, 0);
    krep.assign(nint, 0);
    vecEqLength.resize(2 * nint);


    int nc = Cylinder::getCylinders().size();
    int i = 0;
    int Cumnc=0;
    for (i = 0; i < nc; i++) {
        auto ci = Cylinder::getCylinders()[i];

#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
        
        for (int ID = 0; ID < _HnlIDvec.size(); ID ++){
            auto neighbors = _HneighborList->getNeighborsstencil(_HnlIDvec[ID], ci);
            int nn = neighbors.size();
            //        std::cout<<"Cylinder "<<i<<" "<<nn<<endl;
            for (int ni = 0; ni < nn; ni++) {
                
                auto cin = neighbors[ni];
                if(cin->getBranchingCylinder() == ci) continue;
                beadSet[n * (Cumnc)] = findBeadCoordIndex(*ci->getFirstBead(), si);
                beadSet[n * (Cumnc) + 1] = findBeadCoordIndex(*ci->getSecondBead(), si);
                beadSet[n * (Cumnc) + 2] = findBeadCoordIndex(*cin->getFirstBead(), si);
                beadSet[n * (Cumnc) + 3] = findBeadCoordIndex(*cin->getSecondBead(), si);

                vecEqLength[2 * Cumnc    ] = ci ->getMCylinder()->getEqLength();
                vecEqLength[2 * Cumnc + 1] = cin->getMCylinder()->getEqLength();
                
                //Get KRepuls based on filament type
                if(ci->getType() != cin->getType()){
                    auto ki = ci->getMCylinder()->getExVolConst();
                    auto kin = cin->getMCylinder()->getExVolConst();
                    krep[Cumnc] = max(ki, kin);
                }
                else{
                    krep[Cumnc] = ci->getMCylinder()->getExVolConst();
                }
                
                Cumnc++;
            }
        }
                
#else
        auto neighbors = _neighborList->getNeighbors(ci);
        
        int nn = neighbors.size();
//        std::cout<<"Cylinder "<<i<<" "<<nn<<endl;
        for (int ni = 0; ni < nn; ni++) {

            auto cin = neighbors[ni];
            if(cin->getBranchingCylinder() == ci) continue;
            beadSet[n * (Cumnc)] = findBeadCoordIndex(*ci->getFirstBead(), si);
            beadSet[n * (Cumnc) + 1] = findBeadCoordIndex(*ci->getSecondBead(), si);
            beadSet[n * (Cumnc) + 2] = findBeadCoordIndex(*cin->getFirstBead(), si);
            beadSet[n * (Cumnc) + 3] = findBeadCoordIndex(*cin->getSecondBead(), si);

            vecEqLength[2 * Cumnc    ] = ci ->getMCylinder()->getEqLength();
            vecEqLength[2 * Cumnc + 1] = cin->getMCylinder()->getEqLength();

            //Get KRepuls based on filament type
            if(ci->getType() != cin->getType()){
                auto ki = ci->getMCylinder()->getExVolConst();
                auto kin = cin->getMCylinder()->getExVolConst();
                krep[Cumnc] = max(ki, kin);
            }
            else{
                krep[Cumnc] = ci->getMCylinder()->getExVolConst();
            }

            Cumnc++;
            //std::cout<<"CV"<<ci->getID()<<" "<<cin->getID()<<endl;
        }
#endif
    }
    //CUDA
#ifdef CUDAACCL

    //CUDA stream create
    if(stream == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamCreate(&stream),"cuda stream", "CylinderExclVolume.cu");
    _FFType.optimalblocksnthreads(numInteractions, stream);

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_beadSet, n * numInteractions * sizeof(int)),"cuda data transfer",
                            "CylinderExclVolume.cu");
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_beadSet, beadSet, n * numInteractions *
                                    sizeof(int), cudaMemcpyHostToDevice, stream),
                            "cuda data transfer", "CylinderExclVolume.cu");
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_krep, numInteractions * sizeof(floatingpoint)),"cuda data transfer",
                            "CylinderExclVolume.cu");
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_krep, krep, numInteractions * sizeof
                            (floatingpoint), cudaMemcpyHostToDevice, stream),
                            "cuda data transfer", "CylinderExclVolume.cu");
    vector<int> params;
    params.push_back(int(n));
    params.push_back(numInteractions);
    params.push_back(CUDAcommon::cudavars.offset_E);
//    std::cout<<"CUDA exvol offsetE "<<CUDAcommon::cudavars.offset_E<<endl;
    //set offset
    CUDAcommon::cudavars.offset_E += nint;
//    std::cout<<"offset "<<getName()<<" "<<CUDAcommon::cudavars.offset_E<<endl;
    //TODO make sure not using cudafree here is fine.
    if(gpu_params != NULL )
        CUDAcommon::handleerror(cudaFree(gpu_params),"cudaFree", "CylinderExclVolume.cu");
    if(nint > 0) {
        CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, 3 * sizeof(int)), "cuda"
                                        " data transfer",
                                "CylinderExclVolume.cu");
        CUDAcommon::handleerror(cudaMemcpyAsync(gpu_params, params.data(), 3 * sizeof(int),
                                           cudaMemcpyHostToDevice, stream),
                                "cuda data transfer", "CylinderExclVolume.cu");
    }
#endif
    //
}



template <class CVolumeInteractionType>
FP CylinderExclVolume<CVolumeInteractionType>::computeEnergy(FP *coord) {

#ifdef CUDAACCL


    //has to be changed to accomodate aux force
    floatingpoint * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    floatingpoint * gpu_force = CUDAcommon::getCUDAvars().gpu_force;
    floatingpoint * gpu_d = CUDAcommon::getCUDAvars().gpu_lambda;

//    if(d == 0.0){
//        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_krep, gpu_params);
//    }
//    else{
        gU_i=_FFType.energy(gpu_coord, gpu_force, gpu_beadSet, gpu_krep, gpu_d, gpu_params);
//    }

#endif

    FP U_ii = _FFType.energy(coord, beadSet.data(), krep.data(), vecEqLength.data(), numInteractions);

    return U_ii;
}

template <class CVolumeInteractionType>
void CylinderExclVolume<CVolumeInteractionType>::computeForces(FP *coord, FP *f) {
#ifdef CUDAACCL
    //has to be changed to accomodate aux force
    floatingpoint * gpu_coord=CUDAcommon::getCUDAvars().gpu_coord;
    floatingpoint * gpu_force;
    if(cross_checkclass::Aux) {

        gpu_force = CUDAcommon::getCUDAvars().gpu_forceAux;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_krep, gpu_params);

    }
    else {

        gpu_force = CUDAcommon::getCUDAvars().gpu_force;
        _FFType.forces(gpu_coord, gpu_force, gpu_beadSet, gpu_krep, gpu_params);
    }
#endif
    _FFType.forces(coord, f, beadSet.data(), krep.data(), vecEqLength.data(), numInteractions);
#ifdef DETAILEDOUTPUT
    floatingpoint maxF = 0.0;
    floatingpoint mag = 0.0;
    for(int i = 0; i < CGMethod::N/3; i++) {
        mag = 0.0f;
        for(int j = 0; j < 3; j++)
            mag += f[3 * i + j]*f[3 * i + j];
        mag = sqrt(mag);
//        std::cout<<"SL "<<i<<" "<<mag*mag<<" "<<forceAux[3 * i]<<" "<<forceAux[3 * i + 1]<<" "<<forceAux[3 * i +
//                2]<<endl;
        if(mag > maxF) maxF = mag;
    }
    std::cout<<"max "<<maxF<<endl;
#endif
}

// Explicit instantiation.
template class CylinderExclVolume<CylinderExclVolRepulsion>;

} // namespace medyan

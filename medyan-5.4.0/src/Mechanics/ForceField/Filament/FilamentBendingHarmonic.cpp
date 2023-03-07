
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

#include "FilamentBendingHarmonic.h"
#include "FilamentBendingHarmonicCUDA.h"
#include "FilamentBending.h"

#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"
#include "MathFunctions.h"
#ifdef CUDAACCL
#include <cuda.h>
#include <cuda_runtime.h>
#include "nvToolsExt.h"
#endif

namespace medyan {
using namespace mathfunc;
#ifdef CUDAACCL
void FilamentBendingHarmonic::deallocate(){
//    if(!(CUDAcommon::getCUDAvars().conservestreams))
//        CUDAcommon::handleerror(cudaStreamDestroy(stream));
    CUDAcommon::handleerror(cudaFree(gU_i));
    CUDAcommon::handleerror(cudaFree(gU_sum));
    CUDAcommon::handleerror(cudaFree(gFF));
    CUDAcommon::handleerror(cudaFree(ginteraction));
}
void FilamentBendingHarmonic::optimalblocksnthreads( int nint, cudaStream_t stream_pass){
//    //CUDA stream create
//    if(stream == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
//        CUDAcommon::handleerror(cudaStreamCreate(&stream));
    stream = stream_pass;
    blocksnthreadse.clear();
    blocksnthreadsez.clear();
    blocksnthreadsf.clear();
    int blockSize;   // The launch configurator returned block size
    int minGridSize; // The minimum grid size needed to achieve the
    // maximum occupancy for a full device launch
    if(nint>0) {
        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       FilamentBendingHarmonicenergy, blockToSmemFB, 0);
        blocksnthreadse.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadse.push_back(blockSize);
        blockSize = 0;

        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       FilamentBendingHarmonicenergyz, blockToSmemFB2, 0);
        blocksnthreadsez.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadsez.push_back(blockSize);
        blockSize = 0;

        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       FilamentBendingHarmonicforces, blockToSmemFB, 0);
        blocksnthreadsf.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadsf.push_back(blockSize);
        //get addition vars
        bntaddvec2.clear();
        bntaddvec2 = getaddred2bnt(nint);
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, bntaddvec2.at(0)*sizeof(floatingpoint)));
        CUDAcommon::handleerror(cudaMemset(gU_i, 0, bntaddvec2.at(0) * sizeof(floatingpoint)));
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_sum, sizeof(floatingpoint)));
        char a[] = "FilamentFF";
        char b[] = "Filament Bending Harmonic";
        CUDAcommon::handleerror(cudaMalloc((void **) &gFF, 100 * sizeof(char)));
        CUDAcommon::handleerror(cudaMalloc((void **) &ginteraction, 100 * sizeof(char)));
        CUDAcommon::handleerror(cudaMemcpy(gFF, a, 100 * sizeof(char), cudaMemcpyHostToDevice));
        CUDAcommon::handleerror(cudaMemcpy(ginteraction, b, 100 * sizeof(char), cudaMemcpyHostToDevice));

    }
    else{
        blocksnthreadse.push_back(0);
        blocksnthreadse.push_back(0);
        blocksnthreadsez.push_back(0);
        blocksnthreadsez.push_back(0);
        blocksnthreadsf.push_back(0);
        blocksnthreadsf.push_back(0);
    }

}
floatingpoint* FilamentBendingHarmonic::energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                      floatingpoint *kbend, floatingpoint *eqt, int *params) {
    if(blocksnthreadse[1]>0) {
        FilamentBendingHarmonicenergy<<<blocksnthreadse[0], blocksnthreadse[1], (9 * blocksnthreadse[1]) * sizeof
                (floatingpoint), stream>>> (coord, f, beadSet, kbend, eqt, params, gU_i, CUDAcommon::getCUDAvars().gculpritID,
                CUDAcommon::getCUDAvars().gculpritFF,
                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
        CUDAcommon::handleerror( cudaGetLastError() ,"FilamentBendingHarmonicenergy", "FilamentBendingHarmonic.cu");
        floatingpoint* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        cudaStreamSynchronize(stream);
//        addvectorred<<<1,200,200*sizeof(floatingpoint),stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        cudaStreamSynchronize(stream);
//        std::cout<<"bntaddvec "<<bntaddvec2.at(0)<<" "<<bntaddvec2.at(1)<<" "<<bntaddvec2.at(0)<<" "
//                ""<<bntaddvec2.at(2)<<" "<<bntaddvec2.at(3)<<endl;
        resetfloatingpointvariableCUDA<<<1,1,0,stream>>>(gU_sum);
        addvectorred2<<<bntaddvec2.at(2),bntaddvec2.at(3), bntaddvec2.at(3) * sizeof(floatingpoint),stream>>>(gU_i,
                params, gU_sum, gpu_Utot);
        CUDAcommon::handleerror( cudaGetLastError() ,"FilamentBendingHarmonicenergy", "FilamentBendingHarmonic.cu");
        return gU_sum;}
    else
        return NULL;
}


floatingpoint* FilamentBendingHarmonic::energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                      floatingpoint *kbend, floatingpoint *eqt, floatingpoint *z, int *params) {
    if(blocksnthreadsez[1]>0) {

        FilamentBendingHarmonicenergyz << < blocksnthreadsez[0], blocksnthreadsez[1], (18 * blocksnthreadsez[1]) *
                                            sizeof(floatingpoint), stream>> > (coord, f, beadSet, kbend, eqt, params, gU_i,
                                            z, CUDAcommon::getCUDAvars().gculpritID,
                                            CUDAcommon::getCUDAvars().gculpritFF,
                                            CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
        CUDAcommon::handleerror(cudaGetLastError(),"FilamentBendingHarmonicenergyz", "FilamentBendingHarmonic.cu");
        floatingpoint* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        cudaStreamSynchronize(stream);
//        addvectorred<<<1,200,200*sizeof(floatingpoint),stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        std::cout<<"bntaddvec "<<bntaddvec2.at(0)<<" "<<bntaddvec2.at(1)<<" "<<bntaddvec2.at(0)<<" "
//                ""<<bntaddvec2.at(2)<<" "<<bntaddvec2.at(3)<<endl;
        resetfloatingpointvariableCUDA<<<1,1,0,stream>>>(gU_sum);
        addvectorred2<<<bntaddvec2.at(2),bntaddvec2.at(3), bntaddvec2.at(3) * sizeof(floatingpoint),stream>>>(gU_i,
                params, gU_sum, gpu_Utot);
//        cudaStreamSynchronize(stream);
        CUDAcommon::handleerror(cudaGetLastError(),"FilamentBendingHarmonicenergyz", "FilamentBendingHarmonic.cu");
        return gU_sum;
    }else
        return NULL;
}

void FilamentBendingHarmonic::forces(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                   floatingpoint *kbend, floatingpoint *eqt, int *params){
    if(blocksnthreadsf[1]>0) {
        FilamentBendingHarmonicforces << < blocksnthreadsf[0], blocksnthreadsf[1], (9 * blocksnthreadsf[1]) *
                                                                                 sizeof(floatingpoint), stream >> > (coord, f, beadSet, kbend, eqt, params);
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
        CUDAcommon::handleerror(cudaGetLastError(),"FilamentBendingHarmonicforces", "FilamentBendingHarmonic.cu");
    }
}
void FilamentBendingHarmonic::checkforculprit() {
    CUDAcommon::printculprit("FilamentBending","FilamentBendingCosine");
    Filament* fil;
    int i = 0;
    bool found = false;
    for (auto f: Filament::getFilaments()) {

        if (f->getCylinderVector().size() > 1){
            i = i + 2 * f->getCylinderVector().size() - 2;
            if(i > CUDAcommon::getCUDAvars().culpritID[0] && !found){
                found = true;
                fil = (Filament*)(Cylinder::getCylinders()[i]->getParent());
            }
        }
    }
    cout<<"Printing culprit Filament information."<<endl;
    fil->printSelf();
    exit(EXIT_FAILURE);
}
#endif
floatingpoint FilamentBendingHarmonic::energy(floatingpoint *coord, size_t nint, int *beadSet,
                                       floatingpoint *kbend, floatingpoint *eqt){

    int n = FilamentBending<FilamentBendingHarmonic>::n;

    floatingpoint *coord1, *coord2, *coord3, U_i, L1, L2, L1L2, l1l2;

    floatingpoint U = 0.0;

    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[beadSet[n * i]];
        coord2 = &coord[beadSet[n * i + 1]];
        coord3 = &coord[beadSet[n * i + 2]];

        L1 = sqrt(scalarProduct(coord1, coord2,
                                coord1, coord2));
        L2 = sqrt(scalarProduct(coord2, coord3,
                                coord2, coord3));

        L1L2 = L1*L2;
        l1l2 = scalarProduct(coord1, coord2,
                             coord2, coord3);


        U_i = kbend[i] * ( 1 - l1l2 / L1L2 );

        U += U_i;
    }

    return U;
}


void FilamentBendingHarmonic::forces(floatingpoint *coord, floatingpoint *f, size_t nint, int *beadSet,
                                     floatingpoint *kbend, floatingpoint *eqt){

    int n = FilamentBending<FilamentBendingHarmonic>::n;

    floatingpoint *coord1, *coord2, *coord3, L1, L2, l1l2, invL1, invL2, A,B,C, k;
    floatingpoint *force1, *force2, *force3;

    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[beadSet[n * i]];
        coord2 = &coord[beadSet[n * i + 1]];
        coord3 = &coord[beadSet[n * i + 2]];

        force1 = &f[beadSet[n * i]];
        force2 = &f[beadSet[n * i + 1]];
        force3 = &f[beadSet[n * i + 2]];

        L1 = sqrt(scalarProduct(coord1, coord2,
                                coord1, coord2));
        L2 = sqrt(scalarProduct(coord2, coord3,
                                coord2, coord3));

        l1l2 = scalarProduct(coord1, coord2,
                             coord2, coord3);

        invL1 = 1/L1;
        invL2 = 1/L2;
        A = invL1*invL2;
        B = l1l2*invL1*A*A*L2;
        C = l1l2*invL2*A*A*L1;

        k = kbend[i];

        //force on i-1, f = k*(-A*l2 + B*l1):
        force1[0] +=  k * ((-coord3[0] + coord2[0])*A +
                           (coord2[0] - coord1[0])*B );
        force1[1] +=  k * ((-coord3[1] + coord2[1])*A +
                           (coord2[1] - coord1[1])*B );
        force1[2] +=  k * ((-coord3[2] + coord2[2])*A +
                           (coord2[2] - coord1[2])*B );


        //force on i, f = k*(A*(l1-l2) - B*l1 + C*l2):
        force2[0] +=  k *( (coord3[0] - 2*coord2[0] + coord1[0])*A -
                           (coord2[0] - coord1[0])*B +
                           (coord3[0] - coord2[0])*C );

        force2[1] +=  k *( (coord3[1] - 2*coord2[1] + coord1[1])*A -
                           (coord2[1] - coord1[1])*B +
                           (coord3[1] - coord2[1])*C );

        force2[2] +=  k *( (coord3[2] - 2*coord2[2] + coord1[2])*A -
                           (coord2[2] - coord1[2])*B +
                           (coord3[2] - coord2[2])*C );

        //force on i-1, f = k*(A*l - B*l2):
        force3[0] +=  k *( (coord2[0] - coord1[0])*A -
                           (coord3[0] - coord2[0])*C );

        force3[1] +=  k *( (coord2[1] - coord1[1])*A -
                           (coord3[1] - coord2[1])*C );

        force3[2] +=  k *( (coord2[2] - coord1[2])*A -
                           (coord3[2] - coord2[2])*C );

    }
}

} // namespace medyan

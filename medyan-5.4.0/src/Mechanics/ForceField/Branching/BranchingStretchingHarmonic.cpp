
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

#include "BranchingStretchingHarmonic.h"
#include "BranchingStretchingHarmonicCUDA.h"
#include "BranchingStretching.h"

#include "BranchingPoint.h"
#include "Bead.h"

#include "MathFunctions.h"
#include "Cylinder.h"

#ifdef CUDAACCL
#include <cuda.h>
#include <cuda_runtime.h>
#include "nvToolsExt.h"
#endif

namespace medyan {
using namespace mathfunc;
#ifdef CUDAACCL
void BranchingStretchingHarmonic::deallocate(){
    if(!(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamDestroy(stream));
    CUDAcommon::handleerror(cudaFree(gU_i));
    CUDAcommon::handleerror(cudaFree(gU_sum));
    CUDAcommon::handleerror(cudaFree(gFF));
    CUDAcommon::handleerror(cudaFree(ginteraction));
}
void BranchingStretchingHarmonic::optimalblocksnthreads( int nint){
    //CUDA stream create
    if(stream == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamCreate(&stream));
    blocksnthreadse.clear();
    blocksnthreadsez.clear();
    blocksnthreadsf.clear();
    int blockSize;   // The launch configurator returned block size
    int minGridSize; // The minimum grid size needed to achieve the
    // maximum occupancy for a full device launch
    if(nint>0) {
        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       BranchingStretchingHarmonicenergy, blockToSmemFB, 0);
        blocksnthreadse.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadse.push_back(blockSize);
        blockSize = 0;

        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       BranchingStretchingHarmonicenergyz, blockToSmemFB2, 0);
        blocksnthreadsez.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadsez.push_back(blockSize);
        blockSize = 0;

        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       BranchingStretchingHarmonicforces, blockToSmemFB, 0);
        blocksnthreadsf.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadsf.push_back(blockSize);
        //get addition vars
        bntaddvec2.clear();
        bntaddvec2 = getaddred2bnt(nint);
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, bntaddvec2.at(0)*sizeof(floatingpoint)));
        CUDAcommon::handleerror(cudaMemset(gU_i, 0, bntaddvec2.at(0) * sizeof(floatingpoint)));
//        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, nint*sizeof(floatingpoint)));
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_sum, sizeof(floatingpoint)));

//        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, nint*sizeof(floatingpoint)));
//        CUDAcommon::handleerror(cudaMalloc((void **) &gU_sum, sizeof(floatingpoint)));

        char a[] = "BranchingFF";
        char b[] = "Branching Stretching Harmonic";
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
floatingpoint* BranchingStretchingHarmonic::energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                            floatingpoint *kstr, floatingpoint *eql, floatingpoint *pos, int *params) {
//    if(blocksnthreadse[1]>0) {
//        BranchingStretchingHarmonicenergy<<<blocksnthreadse[0], blocksnthreadse[1], (9 * blocksnthreadse[1]) * sizeof
//                (floatingpoint), stream>>> (coord, f, beadSet, kstr, eql, pos, params, gU_i, CUDAcommon::getCUDAvars()
//                .gculpritID, CUDAcommon::getCUDAvars().gculpritFF, CUDAcommon::getCUDAvars().gculpritinteraction,
//                gFF, ginteraction);
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//        CUDAcommon::handleerror( cudaGetLastError() ,"BranchingStretchingHarmonicenergy", "BranchingStretchingHarmonic.cu");
//        floatingpoint* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror( cudaGetLastError() ,"BranchingStretchingHarmonicenergy", "BranchingStretchingHarmonic.cu");
//        return gU_sum;}
//    else
//        return NULL;
}


floatingpoint* BranchingStretchingHarmonic::energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                            floatingpoint *kstr, floatingpoint *eql, floatingpoint *pos, floatingpoint *z, int *params) {
    if(blocksnthreadse[1]>0) {
        BranchingStretchingHarmonicenergy<<<blocksnthreadse[0], blocksnthreadse[1], (9 * blocksnthreadse[1]) * sizeof
                (floatingpoint), stream>>> (coord, f, beadSet, kstr, eql, pos, params, gU_i, z, CUDAcommon::getCUDAvars()
                .gculpritID, CUDAcommon::getCUDAvars().gculpritFF, CUDAcommon::getCUDAvars().gculpritinteraction,
                gFF, ginteraction);
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//        CUDAcommon::handleerror( cudaGetLastError() ,"BranchingStretchingHarmonicenergy", "BranchingStretchingHarmonic.cu");
//        floatingpoint* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror( cudaGetLastError() ,"BranchingStretchingHarmonicenergy", "BranchingStretchingHarmonic.cu");
//        return gU_sum;
    }

    if(blocksnthreadsez[1]>0) {
        BranchingStretchingHarmonicenergyz << < blocksnthreadsez[0], blocksnthreadsez[1], (18 * blocksnthreadsez[1]) *
                                                                                          sizeof(floatingpoint), stream>> >
                (coord, f, beadSet, kstr, eql, pos, params, gU_i, z, CUDAcommon::getCUDAvars().gculpritID,
                CUDAcommon::getCUDAvars().gculpritFF,
                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//        CUDAcommon::handleerror(cudaGetLastError(),"BranchingStretchingHarmonicenergyz", "BranchingStretchingHarmonic"
//                ".cu");
//        floatingpoint* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//
//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror(cudaGetLastError(),"BranchingStretchingHarmonicenergyz", "BranchingStretchingHarmonic"
//                ".cu");
//        return gU_sum;
    }
    if(blocksnthreadse[1]<=0 && blocksnthreadsez[1]<=0)
        return NULL;
    else{
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
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
//        CUDAcommon::handleerror(cudaDeviceSynchronize(),"FilamentBendingCosineenergyz", "FilamentBendingCosine.cu");
        CUDAcommon::handleerror(cudaGetLastError(),"FilamentBendingCosineenergyz", "FilamentBendingCosine.cu");
        return gU_sum;
    }

}

void BranchingStretchingHarmonic::forces(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                         floatingpoint *kstr, floatingpoint *eql, floatingpoint *pos, int *params){
    if(blocksnthreadsf[1]>0) {
        BranchingStretchingHarmonicforces << < blocksnthreadsf[0], blocksnthreadsf[1], (9 * blocksnthreadsf[1]) *
                                                                                       sizeof(floatingpoint), stream >> >
                (coord, f, beadSet, kstr, eql, pos, params);
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
    }
}
void BranchingStretchingHarmonic::checkforculprit() {
    CUDAcommon::printculprit("BranchingStretching","BranchingStretchingHarmonic");
    BranchingPoint* br;
    br = (BranchingPoint::getBranchingPoints()[CUDAcommon::getCUDAvars().culpritID[0]]);
    cout<<"Printing culprit branching point information."<<endl;
    br->printSelf();
    exit(EXIT_FAILURE);
}
#endif
floatingpoint BranchingStretchingHarmonic::energy(floatingpoint *coord, int *beadSet,
                                           floatingpoint *kstr, floatingpoint *eql, floatingpoint *pos){

    int n = BranchingStretching<BranchingStretchingHarmonic>::n;
    int nint = BranchingPoint::getBranchingPoints().size();

    floatingpoint *coord1, *coord2, *coord3, dist, U_i;
    floatingpoint *v1 = new floatingpoint[3];

    floatingpoint U = 0.0;

    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[beadSet[n * i]];
        coord2 = &coord[beadSet[n * i + 1]];
        coord3 = &coord[beadSet[n * i + 2]];

        midPointCoordinate(v1, coord1, coord2, pos[i]);
        dist = twoPointDistance(v1, coord3) - eql[i];

        U_i = 0.5 * kstr[i] * dist * dist;

        U += U_i;
    }
    delete [] v1;
    return U;
}


void BranchingStretchingHarmonic::forces(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                         floatingpoint *kstr, floatingpoint *eql, floatingpoint *pos,
                                         floatingpoint *stretchforce){


    int n = BranchingStretching<BranchingStretchingHarmonic>::n;
    int nint = BranchingPoint::getBranchingPoints().size();

    floatingpoint *coord1, *coord2, *coord3, dist, invL, f0;
    floatingpoint *f1, *f2, *f3;
    floatingpoint *v1 = new floatingpoint[3];
    floatingpoint r1x, r1y, r1z;


    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[beadSet[n * i]];
        coord2 = &coord[beadSet[n * i + 1]];
        coord3 = &coord[beadSet[n * i + 2]];

        f1 = &f[beadSet[n * i]];
        f2 = &f[beadSet[n * i + 1]];
        f3 = &f[beadSet[n * i + 2]];


        midPointCoordinate(v1, coord1, coord2, pos[i]);
        dist = twoPointDistance(v1, coord3);

        invL = 1 / dist;
        f0 = kstr[i] * ( dist - eql[i]) * invL;

        r1x = coord3[0] - v1[0];
        r1y = coord3[1] - v1[1];
        r1z = coord3[2] - v1[2];

        f1[0] +=  -f0 * ( r1x ) * (pos[i] - 1);
        f1[1] +=  -f0 * ( r1y ) * (pos[i] - 1);
        f1[2] +=  -f0 * ( r1z ) * (pos[i] - 1);

        // force i+1
        f2[0] +=  f0 * ( r1x ) * pos[i];
        f2[1] +=  f0 * ( r1y ) * pos[i];
        f2[2] +=  f0 * ( r1z ) * pos[i];

        //force on j
        f3[0] +=  -f0 * ( r1x );
        f3[1] +=  -f0 * ( r1y );
        f3[2] +=  -f0 * ( r1z );

        stretchforce[3*i] = -f0 * ( r1x );
        stretchforce[3*i + 1] = -f0 * ( r1y );
        stretchforce[3*i + 2] = -f0 * ( r1z );

        #ifdef CHECKFORCES_INF_NAN
        if(checkNaN_INF<floatingpoint>(f1, 0, 2)||checkNaN_INF<floatingpoint>(f2,0,2)||checkNaN_INF<floatingpoint>(f3,0,2)){
            cout<<"Branching Stretching Force becomes infinite. Printing data "<<endl;

            auto b = BranchingPoint::getBranchingPoints()[i];
            auto cyl1 = b->getFirstCylinder();
            auto cyl2 = b->getSecondCylinder();
            cout<<"Cylinder IDs "<<cyl1->getId()<<" "<<cyl2->getId()<<" with cIndex "
                <<cyl1->getStableIndex()<<" "<<cyl2->getStableIndex()<<" and bIndex "
                <<cyl1->getFirstBead()->getStableIndex()<<" "
                <<cyl1->getSecondBead()->getStableIndex()<<" "
                <<cyl2->getFirstBead()->getStableIndex()<<" "
                <<cyl2->getSecondBead()->getStableIndex()<<endl;
            cyl1->adjustedrelativeposition(pos[i], true);

            cout<<"Printing coords"<<endl;
            cout<<coord1[0]<<" "<<coord1[1]<<" "<<coord1[2]<<endl;
            cout<<coord2[0]<<" "<<coord2[1]<<" "<<coord2[2]<<endl;
            cout<<coord3[0]<<" "<<coord3[1]<<" "<<coord3[2]<<endl;

            cout<<"Printing force"<<endl;
            cout<<f1[0]<<" "<<f1[1]<<" "<<f1[2]<<endl;
            cout<<f2[0]<<" "<<f2[1]<<" "<<f2[2]<<endl;
            cout<<f3[0]<<" "<<f3[1]<<" "<<f3[2]<<endl;

            cout<<"Printing binary Coords"<<endl;
            printvariablebinary(coord1,0,2);
            printvariablebinary(coord2,0,2);
            printvariablebinary(coord3,0,2);

            cout<<"Printing binary Force"<<endl;
            printvariablebinary(f1,0,2);
            printvariablebinary(f2,0,2);
            printvariablebinary(f3,0,2);

            exit(EXIT_FAILURE);
        }
        #endif

    }
    delete [] v1;
}

} // namespace medyan

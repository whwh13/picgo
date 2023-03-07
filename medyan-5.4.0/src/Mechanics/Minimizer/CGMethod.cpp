
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

#include "Mechanics/Minimizer/CGMethod.hpp"

#include <algorithm> // max

#include "ForceFieldManager.h"

#include "CGMethodCUDA.h"
#include "MathFunctions.h"
#ifdef CUDAACCL
#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#else
#define CUDA_HOSTDEV
#include "nvToolsExt.h"
#endif
#include <cuda.h>
#include <cuda_runtime.h>
#include "CUDAcommon.h"
#endif
#define ARRAY_SIZE 128
//
#include <vector>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include "cross_check.h"

namespace medyan {

void CGMethod::calcCoordLineSearch(const std::vector<floatingpoint>& coord, floatingpoint d) {
    for(size_t i = 0; i < numDof; ++i)
        coordLineSearch[i] = coord[i] + d * searchDir[i];
}


#ifdef CUDAACCL
void CGMethod::CUDAresetlambda(cudaStream_t stream) {
    resetlambdaCUDA<<<1,1,0, stream>>>(CUDAcommon::getCUDAvars().gpu_lambda);
            CUDAcommon::handleerror(cudaGetLastError(), "resetlambdaCUDA", "CGMethod.cu");
}
void CGMethod::CUDAinitializeLambda(cudaStream_t stream, bool *check_in, bool *check_out, bool *Polaksafestate, int
                                    *gpu_state){

//    cudaStream_t  s;
//    CUDAcommon::handleerror(cudaStreamCreate(&s));
//    cudaEvent_t  e;
//    CUDAcommon::handleerror(cudaEventCreate(&e));


//    maxFCUDA<<<1,1, 0, s>>>(CUDAcommon::getCUDAvars().gpu_forceAux, gpu_nint, gpu_fmax);
//    cudaStreamSynchronize(s);

//    maxFCUDAred<<<1,3, 3*sizeof(floatingpoint), s>>>(CUDAcommon::getCUDAvars().gpu_forceAux, gpu_nint, gpu_fmax);
//    cudaStreamSynchronize(s);


//    CUDAcommon::handleerror(cudaDeviceSynchronize());
//    std::cout<<"======"<<endl;
//    CUDAcommon::handleerror(cudaEventRecord(e,s));
//    CUDAcommon::handleerror(cudaGetLastError(), "maxFCUDA", "CGMethod.cu");

////    CUDAcommon::handleerror(cudaStreamWaitEvent(stream,e,0));

////    CUDAcommon::handleerror(cudaEventDestroy(e));

//    CUDAcommon::handleerror(cudaStreamDestroy(s));

//    auto gpu_lambda = CUDAcommon::getCUDAvars().gpu_lambda;
    auto gpu_energy = CUDAcommon::getCUDAvars().gpu_energy;
    initializeLambdaCUDA<<<1,1,0, stream>>>(check_in, check_out, g_currentenergy, gpu_energy, gpu_initlambdalocal, gpu_fmax,
            gpu_params, Polaksafestate, gpu_state);

//    CUDAcommon::handleerror(cudaStreamSynchronize (stream));

    CUDAcommon::handleerror(cudaGetLastError(), "initializeLambdaCUDA", "CGMethod.cu");
}

//void CGMethod::getmaxFCUDA(floatingpoint *gpu_forceAux, int *gpu_nint, floatingpoint *gpu_fmax) {
//    maxFCUDA<<<1,1>>>(CUDAcommon::getCUDAvars().gpu_forceAux, gpu_nint, gpu_fmax);
//    CUDAcommon::handleerror(cudaGetLastError(), "getmaxFCUDA", "CGMethod.cu");
//}
void CGMethod::CUDAfindLambda(cudaStream_t  stream1, cudaStream_t stream2, cudaEvent_t  event, bool *checkin, bool
        *checkout, bool *gpu_safestate, int *gpu_state) {
//ToDo remove stream2 from the list of args.
    auto gpu_energy = CUDAcommon::getCUDAvars().gpu_energy;
    auto gpu_lambda = CUDAcommon::getCUDAvars().gpu_lambda;
    findLambdaCUDA << < 1, 1, 0, stream1 >> > (gpu_energy, g_currentenergy, gpu_FDotFA, gpu_fmax, gpu_lambda,
            gpu_params, checkin, checkout, gpu_safestate, gpu_state);
    CUDAcommon::handleerror(cudaEventRecord(event, stream1));
#ifdef ALLSYNC
    cudaDeviceSynchronize();
#endif
/*    findLambdaCUDA2 << < 1, 1, 0, stream1 >> > (gpu_fmax, gpu_lambda, gpu_params,
            checkin, checkout, gpu_safestate,
            gpu_state);
    CUDAcommon::handleerror(cudaGetLastError(), "findLambdaCUDA", "CGMethod.cu")*/;
}
//void CGMethod::CUDAprepforbacktracking(cudaStream_t stream, bool *check_in, bool *check_out){

//    cudaStream_t  s;
//    CUDAcommon::handleerror(cudaStreamCreate(&s));
//    cudaEvent_t  e;
//    CUDAcommon::handleerror(cudaEventCreate(&e));

//    maxFCUDA<<<1,1, 0, s>>>(CUDAcommon::getCUDAvars().gpu_forceAux, gpu_nint, gpu_fmax);
//    CUDAcommon::handleerror(cudaEventRecord(e,s));
//    CUDAcommon::handleerror(cudaGetLastError());


//    CUDAcommon::handleerror(cudaStreamWaitEvent(stream,e,0));


//    CUDAcommon::handleerror(cudaEventDestroy(e));
//    CUDAcommon::handleerror(cudaStreamDestroy(s));

////    CUDAcommon::handleerror(cudaStreamSynchronize (stream));
//    auto gpu_lambda = CUDAcommon::getCUDAvars().gpu_lambda;
//    auto gpu_energy = CUDAcommon::getCUDAvars().gpu_energy;
//    prepbacktracking<<<1,1,0, stream>>>(check_in, check_out, g_currentenergy, gpu_energy, gpu_lambda, gpu_fmax,
//            gpu_params);
//    CUDAcommon::handleerror(cudaStreamSynchronize (stream));
//    CUDAcommon::handleerror(cudaGetLastError());
//}
//void CGMethod::CUDAprepforsafebacktracking(cudaStream_t stream, bool *check_in, bool *check_out){
//    auto gpu_lambda = CUDAcommon::getCUDAvars().gpu_lambda;
//    auto gpu_energy = CUDAcommon::getCUDAvars().gpu_energy;
//    prepsafebacktracking<<<1,1,0,stream>>>(check_in, check_out, g_currentenergy, gpu_energy, gpu_lambda, gpu_params);
//    CUDAcommon::handleerror(cudaStreamSynchronize (stream));
//    CUDAcommon::handleerror(cudaGetLastError());
//}

void CGMethod::CUDAallFDotF(cudaStream_t stream){

    allFADotFCUDA<<<blocksnthreads[0], blocksnthreads[1],0,stream>>>(CUDAcommon::getCUDAvars().gpu_force,
            CUDAcommon::getCUDAvars().gpu_force ,gpu_g, gpu_nint);
    CUDAcommon::handleerror(cudaGetLastError(), "allFADotFCUDA", "CGMethod.cu");
//    addvector<<<1,1,0,stream>>>(gpu_g, gpu_nint, gpu_FDotF);
//    cudaStreamSynchronize(stream);
//    addvectorred<<<1,200,200 * sizeof(floatingpoint),stream>>>(gpu_g, gpu_nint, gpu_FDotF);
//    floatingpoint Sum[1];
//        CUDAcommon::handleerror(cudaMemcpy(Sum, gpu_FDotF, sizeof(floatingpoint), cudaMemcpyDeviceToHost));
    resetfloatingpointvariableCUDA<<<1,1,0,stream>>>(gpu_FDotF);
    addvectorredcgm<<<bntaddvector.at(2),bntaddvector.at(3), bntaddvector.at(3) * sizeof(floatingpoint),stream>>>(gpu_g,
            gpu_nint, gpu_FDotF);
//    floatingpoint Sum2[1];
//    CUDAcommon::handleerror(cudaMemcpy(Sum2, gpu_FDotF, sizeof(floatingpoint), cudaMemcpyDeviceToHost));
//    std::cout<<Sum[0]<<" "<<Sum2[0]<<endl;
//    cudaStreamSynchronize(stream);
    CUDAcommon::handleerror(cudaGetLastError(), "allFADotFCUDA", "CGMethod.cu");

}
void CGMethod::CUDAallFADotFA(cudaStream_t stream){

    allFADotFCUDA<<<blocksnthreads[0], blocksnthreads[1],0,stream>>>(CUDAcommon::getCUDAvars().gpu_forceAux,
            CUDAcommon::getCUDAvars().gpu_forceAux ,gpu_g, gpu_nint);
    CUDAcommon::handleerror(cudaGetLastError(), "allFADotFCUDA", "CGMethod.cu");
//    addvector<<<1,1,0,stream>>>(gpu_g, gpu_nint, gpu_FADotFA);
//    cudaStreamSynchronize(stream);
//    addvectorred<<<1,200,200* sizeof(floatingpoint),stream>>>(gpu_g, gpu_nint, gpu_FADotFA);
//    cudaStreamSynchronize(stream);
    resetfloatingpointvariableCUDA<<<1,1,0,stream>>>(gpu_FADotFA);
    addvectorredcgm<<<bntaddvector.at(2),bntaddvector.at(3), bntaddvector.at(3) * sizeof(floatingpoint),stream>>>(gpu_g,
            gpu_nint, gpu_FADotFA);
    CUDAcommon::handleerror(cudaGetLastError(), "allFADotFCUDA", "CGMethod.cu");

}
void CGMethod::CUDAallFADotFAP(cudaStream_t stream){

    allFADotFCUDA<<<blocksnthreads[0], blocksnthreads[1],0,stream>>>(CUDAcommon::getCUDAvars().gpu_forceAux,
            CUDAcommon::getCUDAvars().gpu_forceAuxP ,gpu_g, gpu_nint);
    CUDAcommon::handleerror(cudaGetLastError(), "allFADotFCUDA", "CGMethod.cu");
//    addvector<<<1,1,0,stream>>>(gpu_g, gpu_nint, gpu_FADotFAP);
//    cudaStreamSynchronize(stream);
//    addvectorred<<<1,200,200 * sizeof(floatingpoint),stream>>>(gpu_g, gpu_nint, gpu_FADotFAP);
    resetfloatingpointvariableCUDA<<<1,1,0,stream>>>(gpu_FADotFAP);
    addvectorredcgm<<<bntaddvector.at(2),bntaddvector.at(3), bntaddvector.at(3) * sizeof(floatingpoint),stream>>>(gpu_g,
            gpu_nint, gpu_FADotFAP);
//    cudaStreamSynchronize(stream);
    CUDAcommon::handleerror(cudaGetLastError(), "allFADotFCUDA", "CGMethod.cu");

}
void CGMethod::CUDAallFDotFA(cudaStream_t stream){

    allFADotFCUDA<<<blocksnthreads[0], blocksnthreads[1],0,stream>>>(CUDAcommon::getCUDAvars().gpu_force,
            CUDAcommon::getCUDAvars().gpu_forceAux ,gpu_g, gpu_nint);
    CUDAcommon::handleerror(cudaGetLastError(), "allFADotFCUDA", "CGMethod.cu");
//    addvector<<<1,1,0,stream>>>(gpu_g, gpu_nint, gpu_FDotFA);
//    cudaStreamSynchronize(stream);
//    addvectorred<<<1,200,200* sizeof(floatingpoint),stream>>>(gpu_g, gpu_nint, gpu_FDotFA);
//    cudaStreamSynchronize(stream);
    resetfloatingpointvariableCUDA<<<1,1,0,stream>>>(gpu_FDotFA);
    addvectorredcgm<<<bntaddvector.at(2),bntaddvector.at(3), bntaddvector.at(3) * sizeof(floatingpoint),stream>>>(gpu_g,
            gpu_nint, gpu_FDotFA);
    CUDAcommon::handleerror(cudaGetLastError(), "allFADotFCUDA", "CGMethod.cu");

}

void CGMethod::CUDAshiftGradient(cudaStream_t stream, bool *Mcheckin) {
    shiftGradientCUDA<<<blocksnthreads[0], blocksnthreads[1],0, stream>>>(CUDAcommon::getCUDAvars().gpu_force,
            CUDAcommon::getCUDAvars().gpu_forceAux, gpu_nint, gpu_FADotFA, gpu_FADotFAP, gpu_FDotF, Mcheckin);
}

void CGMethod::CUDAshiftGradientifSafe(cudaStream_t stream, bool *Mcheckin, bool *Scheckin){
    shiftGradientCUDAifsafe<<<blocksnthreads[0], blocksnthreads[1],0, stream>>>(CUDAcommon::getCUDAvars().gpu_force, CUDAcommon::getCUDAvars().gpu_forceAux, gpu_nint,
                            Mcheckin, Scheckin);
    CUDAcommon::handleerror(cudaGetLastError(),"CUDAshiftGradientifSafe", "CGMethod.cu");
}

//void CGMethod::CUDAgetPolakvars(bool calc_safestate,cudaStream_t streamcalc, floatingpoint* gpu_GRADTOL, bool *gminstatein,
//                                    bool *gminstateout, bool *gsafestateout, volatile bool *cminstate){
////    state[0] = false;
////    state[1] = false;
//    if(cminstate[0] == true) {

////        maxFCUDA << < 1, 1, 0, streamcalc >> > (CUDAcommon::getCUDAvars().gpu_forceAux, gpu_nint, gpu_fmax);
//        maxFCUDAred<<<1,3, 3*sizeof(floatingpoint), streamcalc>>>(CUDAcommon::getCUDAvars().gpu_forceAux, gpu_nint, gpu_fmax);
////        CUDAcommon::handleerror(cudaDeviceSynchronize());
////        std::cout<<"======"<<endl;
//        CUDAcommon::handleerror(cudaGetLastError(), "maxFCUDA", "CGMethod.cu");

//        getminimizestateCUDA << < 1, 1, 0, streamcalc >> > (gpu_fmax, gpu_GRADTOL, gminstatein, gminstateout);
//        CUDAcommon::handleerror(cudaGetLastError(), "getminimizestateCUDA", "CGMethod.cu");
////        CUDAcommon::handleerror(cudaStreamSynchronize(streamcalc));
//    }
//    if(calc_safestate){
//        CUDAallFDotFA(streamcalc);
//        getsafestateCUDA<<<1,1,0,streamcalc>>>(gpu_FDotFA, gpu_FDotF, gpu_FADotFA, gsafestateout);
//        CUDAcommon::handleerror(cudaGetLastError(), "getsafestateCUDA", "CGMethod.cu");
//    }
//}

void CGMethod::CUDAgetPolakvars(cudaStream_t streamcalc, floatingpoint* gpu_GRADTOL, bool *gminstatein,
                                bool *gminstateout, volatile bool *cminstate){
//    state[0] = false;
//    state[1] = false;
    if(cminstate[0] == true) {

        //@{ V3
        //TODO combine with FADOTFA calculation by making it write it to gpu_maxF before
        // adding.
        allFADotFCUDA<<<blocksnthreads[0], blocksnthreads[1],0,streamcalc>>>(CUDAcommon::getCUDAvars().gpu_forceAux,
                CUDAcommon::getCUDAvars().gpu_forceAux ,gpu_maxF, gpu_nint);
        resetfloatingpointvariableCUDA<<<1,1,0,streamcalc>>>(gpu_fmax);
        resetintvariableCUDA<<<1,1,0,streamcalc>>>(gpu_mutexlock);
        maxFCUDAredv3<<<bntaddvector.at(2),bntaddvector.at(3), bntaddvector.at(3) *
                sizeof(floatingpoint),streamcalc>>>(gpu_maxF, gpu_nint, gpu_fmax, gpu_mutexlock);
        CUDAcommon::handleerror(cudaGetLastError(), "maxFCUDA", "CGMethod.cu");
        getminimizestateCUDA << < 1, 1, 0, streamcalc >> > (gpu_fmax, gpu_GRADTOL, gminstatein, gminstateout);
        CUDAcommon::handleerror(cudaGetLastError(), "getminimizestateCUDA", "CGMethod.cu");

    CUDAcommon::handleerror(cudaGetLastError(),"CUDAgetPolakvars", "CGMethod.cu");
}

void CGMethod::CUDAgetPolakvars2(cudaStream_t streamcalc, bool *gsafestateout){
        CUDAallFDotFA(streamcalc);
        getsafestateCUDA<<<1,1,0,streamcalc>>>(gpu_FDotFA, gpu_FDotF, gpu_FADotFA, gsafestateout);
        CUDAcommon::handleerror(cudaGetLastError(), "getsafestateCUDA", "CGMethod.cu");
}

void CGMethod::CUDAmoveBeads(cudaStream_t stream, bool *gpu_checkin){
    floatingpoint *gpu_lambda = CUDAcommon::getCUDAvars().gpu_lambda;
    floatingpoint *gpu_coord = CUDAcommon::getCUDAvars().gpu_coord;
    floatingpoint *gpu_force = CUDAcommon::getCUDAvars().gpu_force;

    moveBeadsCUDA<<<blocksnthreads[0], blocksnthreads[1],0, stream>>>(gpu_coord, gpu_force, gpu_lambda, gpu_nint,
            gpu_checkin);

    CUDAcommon::handleerror(cudaGetLastError(),"moveBeadsCUDA", "CGMethod.cu");
}

void CGMethod::CUDAinitializePolak(cudaStream_t stream, bool *minstatein, bool *minstateout, bool *safestatein, bool
        *safestateout){
    CUDAallFDotFA(stream);
    initializePolak<<<1,1,0,stream>>>(minstatein, minstateout, safestatein, safestateout);
    CUDAcommon::handleerror(cudaGetLastError(),"CUDAinitializePolak", "CGPolakRibiereMethod.cu");
}

//floatingpoint CGMethod::gpuFDotF(floatingpoint *f1,floatingpoint *f2){
//
//    allFADotFCUDA<<<blocksnthreads[0], blocksnthreads[1]>>>(f1, f2 ,gpu_g, gpu_nint);
//    CUDAcommon::handleerror(cudaGetLastError(),"allFADotFCUDA", "CGMethod.cu");
////    addvector<<<1,1>>>(gpu_g, gpu_nint, gSum);
//    addvectorred<<<1,200,200* sizeof(floatingpoint)>>>(gpu_g, gpu_nint, gSum);
//    CUDAcommon::handleerror(cudaGetLastError(),"allFADotFCUDA", "CGMethod.cu");
//
////    CUDAcommon::handleerror( cudaPeekAtLastError() );
////    CUDAcommon::handleerror(cudaDeviceSynchronize());
//
//    floatingpoint g[1];
//    CUDAcommon::handleerror(cudaMemcpy(g, gSum, sizeof(floatingpoint),
//                                       cudaMemcpyDeviceToHost));
//
//
////    floatingpoint g[N/3];
////    CUDAcommon::handleerror(cudaMemcpy(g, gpu_g, N/3 * sizeof(floatingpoint),
////                                       cudaMemcpyDeviceToHost));
////    CUDAcommon::handleerror(cudaFree(gpu_g));
////    floatingpoint sum=0.0;
////    for(auto i=0;i<N/3;i++)
////        sum+=g[i];
//    return g[0];
//}
#endif
double CGMethod::searchDirDotSearchDir() const {
    double res = 0.0;
    for(auto x : searchDir) {
        res += x * x;
    }
    return res;
}

double CGMethod::forceDotForce() const {
    double res = 0.0;
    for(auto x : force) {
        res += x * x;
    }
    return res;
}

double CGMethod::forceDotForcePrev() const {
    double res = 0.0;
    for(std::size_t i = 0; i < numDof; ++i) {
        res += force[i] * forcePrev[i];
    }
    return res;
}

double CGMethod::searchDirDotForce() const {
    double res = 0.0;
    for(std::size_t i = 0; i < numDof; ++i) {
        res += searchDir[i] * force[i];
    }
    return res;
}


void CGMethod::moveAlongSearchDir(std::vector<floatingpoint>& coord, floatingpoint d)
{
    for(std::size_t i = 0; i < numDof; ++i)
        coord[i] += d * searchDir[i];
}

void CGMethod::shiftSearchDir(double d)
{
    for(std::size_t i = 0; i < numDof; ++i)
        searchDir[i] = force[i] + d * searchDir[i];
}

void CGMethod::startMinimization() {

    // Reset backup coordinates with minimum Energy
    coordMinE.clear();
    minimumE = std::numeric_limits<floatingpoint>::infinity();

	long Ncyl = Cylinder::getCylinders().size();

#ifdef CUDAACCL
    //Start stream
    if(stream_startmin == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamCreate(&stream_startmin));
    int nDevices;
//    cudaDeviceProp prop;
    cudaGetDeviceCount(&nDevices);
    if(nDevices>1){
        cout<<"Code not configured for multiple devices. Exiting..."<<endl;
        exit(EXIT_FAILURE);
    }

    floatingpoint f[N];
    for(auto iter=0;i<N;i++)
        f[iter]=0.0;
    floatingpoint* gpu_coord;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_coord, N*sizeof(floatingpoint)));
    floatingpoint* gpu_lambda;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_lambda, sizeof(floatingpoint)));
    floatingpoint* gpu_force;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_force, N*sizeof(floatingpoint)));
    floatingpoint* gpu_forceAux;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_forceAux, N*sizeof(floatingpoint)));
    floatingpoint* gpu_forceAuxP;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_forceAuxP, N*sizeof(floatingpoint)));
    floatingpoint* gpu_energy;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_energy, sizeof(floatingpoint)));
    bool* gpu_btstate;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_btstate, sizeof(bool)));
    cylinder* gpu_cylindervec;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_cylindervec, Ncyl*sizeof(cylinder)));

    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_initlambdalocal, sizeof(floatingpoint)));

    CUDAcommon::handleerror(cudaMalloc((void **)&gpu_fmax, sizeof(floatingpoint)));
    CUDAcommon::handleerror(cudaMalloc((void **)&g_currentenergy, sizeof(floatingpoint)));
    CUDAcommon::handleerror(cudaMalloc((void **)&gpu_FDotF, sizeof(floatingpoint)));
    CUDAcommon::handleerror(cudaMalloc((void **)&gpu_FADotFA, sizeof(floatingpoint)));
    CUDAcommon::handleerror(cudaMalloc((void **)&gpu_FADotFAP, sizeof(floatingpoint)));
    CUDAcommon::handleerror(cudaMalloc((void **)&gpu_FDotFA, sizeof(floatingpoint)));

    CUDAcommon::handleerror(cudaHostAlloc((void**)&convergencecheck, 3 * sizeof(bool), cudaHostAllocMapped));
    CUDAcommon::handleerror(cudaHostGetDevicePointer(&gpu_convergencecheck, convergencecheck, 0));

    //PING PONG
    CUDAcommon::handleerror(cudaMalloc(&g_stop1, sizeof(bool)));
    CUDAcommon::handleerror(cudaMalloc(&g_stop2, sizeof(bool)));
    CUDAcommon::handleerror(cudaHostAlloc(&h_stop, sizeof(bool), cudaHostAllocDefault));

    //@
    //Store the pointers so they can be tracked while calculating energies.
    CUDAcommon::cudavars.backtrackbools.clear();
    CUDAcommon::cudavars.backtrackbools.push_back(g_stop1);
    CUDAcommon::cudavars.backtrackbools.push_back(g_stop2);

//    CUDAcommon::handleerror(cudaHostAlloc((void**)&convergencecheck, 3 * sizeof(bool), cudaHostAllocDefault));
//    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_convergencecheck, 3 * sizeof(bool)));

//    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_lambda, sizeof(floatingpoint))); REPEAT.
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_coord, coord, N*sizeof(floatingpoint),
                                        cudaMemcpyHostToDevice, stream_startmin));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_force, f, N*sizeof(floatingpoint),
                                        cudaMemcpyHostToDevice, stream_startmin));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_forceAux, f, N*sizeof(floatingpoint),
                                        cudaMemcpyHostToDevice, stream_startmin));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_forceAuxP, f, N*sizeof(floatingpoint),
                                        cudaMemcpyHostToDevice, stream_startmin));
    bool dummy[1];dummy[0] = true;
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_btstate, dummy, sizeof(bool),
                                        cudaMemcpyHostToDevice, stream_startmin));

    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_cylindervec, cylindervec, Ncyl*sizeof
                                                    (cylinder),
                                            cudaMemcpyHostToDevice, stream_startmin));
    int *gculpritID;
    char *gculpritFF;
    char *gculpritinteraction;
    int *culpritID;
    char *culpritFF;
    char *culpritinteraction;
    CUDAcommon::handleerror(cudaHostAlloc((void**)&culpritID, 4 * sizeof(int), cudaHostAllocMapped));
    CUDAcommon::handleerror(cudaHostAlloc((void**)&culpritFF, 100*sizeof(char), cudaHostAllocMapped));
    CUDAcommon::handleerror(cudaHostAlloc((void**)&culpritinteraction, 100*sizeof(char), cudaHostAllocMapped));
    CUDAcommon::handleerror(cudaHostGetDevicePointer(&gculpritID, culpritID, 0));
    CUDAcommon::handleerror(cudaHostGetDevicePointer(&gculpritFF, culpritFF, 0));
    CUDAcommon::handleerror(cudaHostGetDevicePointer(&gculpritinteraction, culpritinteraction, 0));
//    CUDAcommon::handleerror(cudaMalloc((void **) &gculpritID, sizeof(int)));
//    CUDAcommon::handleerror(cudaMalloc((void **) &gculpritFF, 11*sizeof(char)));
//    char a[] = "FilamentFF";
//    CUDAcommon::handleerror(cudaMemcpy(gculpritFF, a, 100 * sizeof(char), cudaMemcpyHostToDevice));
//    CUDAcommon::handleerror(cudaMalloc((void **) &gculpritinteraction, 100*sizeof(char)));

    CUDAvars cvars=CUDAcommon::getCUDAvars();
    cvars.gpu_coord=gpu_coord;
    cvars.gpu_lambda=gpu_lambda;
    cvars.gpu_forceAux = gpu_forceAux;
    cvars.gpu_force = gpu_force;
    cvars.gpu_forceAuxP = gpu_forceAuxP;
    cvars.gpu_energy = gpu_energy;
    cvars.gculpritID = gculpritID;
    cvars.culpritID = culpritID;
    cvars.gculpritinteraction = gculpritinteraction;
    cvars.gculpritFF = gculpritFF;
    cvars.culpritinteraction = culpritinteraction;
    cvars.culpritFF = culpritFF;
    cvars.gpu_btstate = gpu_btstate;
    cvars.gpu_cylindervec = gpu_cylindervec;
    CUDAcommon::cudavars=cvars;
//SET CERTAIN GPU PARAMETERS SET FOR EASY ACCESS DURING MINIMIZATION._
//    int THREADSPERBLOCK;
//    cudaDeviceProp prop;
//    cudaGetDeviceProperties(&prop, 0);
//    THREADSPERBLOCK = prop.maxThreadsPerBlock;
    //@{ Reduction Add variables
    bntaddvector.clear();
    bntaddvector = getaddred2bnt(N/3);
    int M = bntaddvector.at(0);
    vector<floatingpoint> zerovec(M);
    fill(zerovec.begin(),zerovec.begin()+M,0.0);
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_g, M * sizeof(floatingpoint)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_g, zerovec.data(),
                           M * sizeof(floatingpoint), cudaMemcpyHostToDevice, stream_startmin));
    /*CUDAcommon::handleerror(cudaMemsetAsync(gpu_g, 0, M * sizeof(floatingpoint), stream_startmin));*/
    //MaxF
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_maxF, M * sizeof(floatingpoint)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_maxF, zerovec.data(),
                                            M * sizeof(floatingpoint), cudaMemcpyHostToDevice, stream_startmin));
    /*CUDAcommon::handleerror(cudaMemsetAsync(gpu_maxF, 0, M * sizeof(floatingpoint), stream_startmin));*/
    int THREADSPERBLOCK = bntaddvector.at(1);
    //@}

    int nint[1]; nint[0]=CGMethod::N/3;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_nint, sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_nint, nint, sizeof(int),
                                        cudaMemcpyHostToDevice, stream_startmin));
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_state, sizeof(int)));
    blocksnthreads.push_back(CGMethod::N/(3*THREADSPERBLOCK) + 1);
    if(blocksnthreads[0]==1) blocksnthreads.push_back(CGMethod::N/3);
    else blocksnthreads.push_back(THREADSPERBLOCK);
    auto maxthreads = 8 * THREADSPERBLOCK;

    //@{maxFredv3
    int state[1];state[0] = 0;
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_mutexlock, sizeof(int)));
    CUDAcommon::handleerror(cudaMemcpyAsync(gpu_mutexlock, state, sizeof(int),
                                        cudaMemcpyHostToDevice, stream_startmin));
    //Synchronize
    CUDAcommon::handleerror(cudaStreamSynchronize(stream_startmin),"CGMethod.cu",
                            "startMinimization");

#endif
}

void CGMethod::endMinimization() {

#ifdef CUDAACCL

    CUDAcommon::handleerror(cudaMemcpy(coord, CUDAcommon::getCUDAvars().gpu_coord, N *
                            sizeof(floatingpoint), cudaMemcpyDeviceToHost));
    CUDAcommon::handleerror(cudaMemcpy(force, CUDAcommon::getCUDAvars().gpu_force, N *
                            sizeof(floatingpoint), cudaMemcpyDeviceToHost));
//    CUDAcommon::handleerror(cudaMemcpy(forceAux, CUDAcommon::getCUDAvars().gpu_forceAux, N *
//                            sizeof(floatingpoint), cudaMemcpyDeviceToHost));

#endif

}

#ifdef CUDAACCL
floatingpoint CGMethod::backtrackingLineSearchCUDA(ForceFieldManager& FFM, floatingpoint MAXDIST,
                                        floatingpoint LAMBDAMAX, bool *gpu_safestate) {
    //@{ Lambda phase 1
    floatingpoint lambda;
    h_stop[0] = false;
    if(s1 == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamCreate(&s1));
    if(s2 == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamCreate(&s2));
    if(s3 == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamCreate(&s3));
    if(e1 == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaEventCreate(&e1));
    if(e2 == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaEventCreate(&e2));
    sp1 = &s1;
    sp2 = &s2;
    ep1 = &e1;
    ep2 = &e2;
    g_s1 = g_stop1;
    g_s2 = g_stop2;
    //prep for backtracking.
    if(gpu_params == NULL){
        //TODO move gpu_params copy permanently out of the function.
        floatingpoint params[5];
        params[0] = BACKTRACKSLOPE;
        params[1] = LAMBDAREDUCE;
        params[2] = LAMBDATOL;
        params[3] = LAMBDAMAX;
        params[4] = MAXDIST;
        CUDAcommon::handleerror(cudaMalloc((void **) &gpu_params, 5 * sizeof(floatingpoint)));
        CUDAcommon::handleerror(cudaMemcpy(gpu_params, params, 5 * sizeof(floatingpoint),
                                           cudaMemcpyHostToDevice));
    }
    CUDAresetlambda(*sp1);//set lambda to zero.
    if(e == NULL || !(CUDAcommon::getCUDAvars().conservestreams))  {
        CUDAcommon::handleerror(cudaEventCreate(&e));
    }

    CUDAcommon::handleerror(cudaEventRecord(e, *sp1));
    auto cvars = CUDAcommon::getCUDAvars();
    cvars.streamvec.clear();
    CUDAcommon::cudavars = cvars;
    //initialize lambda search
    CUDAinitializeLambda(*sp1, g_s1, g_s2, gpu_safestate, gpu_state);
    //@} Lambda phase 1
    //Calculate current energy.
    floatingpoint currentEnergy = FFM.computeEnergy(coord, force, 0.0);
    //wait for energies to be calculated
    for(auto strm:CUDAcommon::getCUDAvars().streamvec)
        CUDAcommon::handleerror(cudaStreamSynchronize(*strm),"backConvSync","CGMethod.cu");
    if(stream_bt == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamCreate(&stream_bt),"find lambda", "CGMethod.cu");

#ifdef DETAILEDOUTPUT_ENERGY
//    CUDAcommon::handleerror(cudaDeviceSynchronize());
    floatingpoint cuda_energy[1];
    CUDAcommon::handleerror(cudaMemcpy(cuda_energy, CUDAcommon::cudavars.gpu_energy,  sizeof(floatingpoint),
                                       cudaMemcpyDeviceToHost));
    std::cout<<"Total Energy cE pN.nm CUDA "<<cuda_energy[0]<<" SERL "<<currentEnergy<<endl;
    std::cout<<endl;
#endif

    //@{ Lambda phase 1b
    cudaStreamSynchronize(*sp1);
    setcurrentenergy<<<1,1,0,*sp1>>>(CUDAcommon::getCUDAvars().gpu_energy, g_currentenergy, CUDAcommon::getCUDAvars()
            .gpu_lambda, gpu_initlambdalocal);
    CUDAcommon::handleerror(cudaGetLastError(),"setcurrentenergy", "CGMethod.cu");
    cudaStreamSynchronize(*sp1);

    //check if converged.
    //TODO commented coz this line is not needed
//    CUDAcommon::handleerror(cudaStreamWaitEvent(s3, *ep1, 0));
//    CUDAcommon::handleerror(cudaEventRecord(*CUDAcommon::getCUDAvars().event, *sp1));
    CUDAcommon::handleerror(cudaMemcpyAsync(h_stop, g_s2, sizeof(bool), cudaMemcpyDeviceToHost, s3));
//    CUDAcommon::handleerror(cudaStreamSynchronize (*sp1)); CHECK IF NEEDED
    cconvergencecheck = h_stop;
    int iter = 0;
    //@} Lambda phase 1b

    while(!(cconvergencecheck[0])) {
        //@{ Lambda phase 2
        iter++;
        CUDAcommon::handleerror(cudaStreamWaitEvent(*sp2, *ep1, 0));
        CUDAcommon::handleerror(cudaStreamSynchronize(*sp2));
        //ping pong swap
        sps = sp1;
        sp1 = sp2;
        sp2 = sps;
        eps = ep1;
        ep1 = ep2;
        ep2 = eps;
//        g_bs = g_b1;
//        g_b1 = g_b2;
//        g_b2 = g_bs;
        g_ss = g_s1;
        g_s1 = g_s2;
        g_s2 = g_ss;

        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.clear();
//        cvars.event = ep1;
        CUDAcommon::cudavars = cvars;
        //@} Lambda phase 2

#ifdef SERIAL_CUDACROSSCHECK
        floatingpoint cuda_lambda[1];
        CUDAcommon::handleerror(cudaDeviceSynchronize(),"CGPolakRibiereMethod.cu","CGPolakRibiereMethod.cu");
        CUDAcommon::handleerror(cudaMemcpy(cuda_lambda, CUDAcommon::cudavars.gpu_lambda,  sizeof(floatingpoint),
                                           cudaMemcpyDeviceToHost));
        lambda = cuda_lambda[0];
#endif

        //TODO let each forcefield calculate energy IFF conv state = false. That will help
        // them avoid unnecessary iterations.
        //let each forcefield also add energies to two different energy variables.
        floatingpoint energyLambda = FFM.computeEnergy(coord, force, lambda);

        //wait for energies to be calculated
         for(auto strm:CUDAcommon::getCUDAvars().streamvec) {
            CUDAcommon::handleerror(cudaStreamSynchronize(*strm), "backConvsync", "CGMethod.cu");
        }
#ifdef SERIAL_CUDACROSSCHECK
        for(auto strm:CUDAcommon::getCUDAvars().streamvec) {
            CUDAcommon::handleerror(cudaStreamSynchronize(*strm), "backConvsync", "CGMethod.cu");
        }
        CUDAcommon::handleerror(cudaDeviceSynchronize());
        floatingpoint cuda_energy[1];
        CUDAcommon::handleerror(cudaMemcpy(cuda_energy, CUDAcommon::cudavars.gpu_energy,  sizeof(floatingpoint),
                                           cudaMemcpyDeviceToHost));
        std::cout<<"Total Energy EL pN.nm CUDA "<<cuda_energy[0]<<" SERL "
                ""<<energyLambda<<endl;
        std::cout<<endl;
#endif
        //@{ Lambda phase 2
        if(!(cconvergencecheck[0])){
            CUDAcommon::handleerror(cudaStreamSynchronize(stream_bt));
            CUDAfindLambda(*sp1, stream_bt, *ep1, g_s1, g_s2, gpu_safestate, gpu_state);
            CUDAcommon::handleerror(cudaStreamSynchronize(*sp1));
            CUDAcommon::handleerror(cudaStreamSynchronize(stream_bt));
            if(cconvergencecheck[0]  == false){
                CUDAcommon::handleerror(cudaStreamWaitEvent(s3, *ep1, 0));
                CUDAcommon::handleerror(cudaMemcpyAsync(h_stop, g_s2, sizeof(bool), cudaMemcpyDeviceToHost, s3));
            }
        }
        //@Lambda phase 2
    }
    if(!(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaFree(gpu_params), "CudaFree", "CGMethod.cu");

    CUDAcommon::handleerror(cudaStreamSynchronize(stream_bt));
    CUDAcommon::handleerror(cudaStreamSynchronize(s1));
    CUDAcommon::handleerror(cudaStreamSynchronize(s2));
    CUDAcommon::handleerror(cudaStreamSynchronize(s3));
    //@} Lambda phase 3
    if(!(CUDAcommon::getCUDAvars().conservestreams))  {
        CUDAcommon::handleerror(cudaStreamDestroy(s1));
        CUDAcommon::handleerror(cudaStreamDestroy(s2));
        CUDAcommon::handleerror(cudaStreamDestroy(s3));
        CUDAcommon::handleerror(cudaStreamDestroy(stream_bt));
        CUDAcommon::handleerror(cudaEventDestroy(e1));
        CUDAcommon::handleerror(cudaEventDestroy(e2));
    }
    std::cout<<"CUDA lambda determined in "<<iter<< " iterations "<<endl;

    if(cconvergencecheck[0]||sconvergencecheck)
        return lambda;

}
#endif // CUDAACCL




ofstream CGMethod::_crosscheckdumpMechFile;

} // namespace medyan

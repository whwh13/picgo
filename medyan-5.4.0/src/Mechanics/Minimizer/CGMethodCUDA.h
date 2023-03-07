//
// Created by aravind on 11/1/17.
//

#ifndef CUDA_VEC_CGMETHODCUDA_H
#define CUDA_VEC_CGMETHODCUDA_H
#ifdef CUDAACCL
#include <cuda.h>
#include <cuda_runtime.h>
#include "math.h"
#include "utility.h"
#include "assert.h"

__global__ void correctlambdaCUDA(floatingpoint *gpu_lambda, int *gpu_state, floatingpoint *gpu_params){
#ifdef DETAILEDOUTPUT_LAMBDA
    printf("gpu_state %d\n",gpu_state[0]);
#endif
    if(gpu_state[0] == 1 || gpu_state[0] == 3 ) {
        gpu_lambda[0] = gpu_lambda[0] / gpu_params[1];
        printf(" state Y %d lambda %f params %f \n", gpu_state[0], gpu_lambda[0],
               gpu_params[1]);
    }
    else
        printf(" state N %d lambda %f \n", gpu_state[0], gpu_lambda[0]);
#ifdef DETAILEDOUTPUT_LAMBDA
    printf("final lambda %f\n", gpu_lambda[0]);
#endif
}

__global__ void moveBeadsCUDA(floatingpoint *coord, floatingpoint* f, floatingpoint *d,  int *nint, bool *checkin) {
    if(checkin[0] == false) return; //if it is not in minimization state
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;
//    if(thread_idx == 0)
//        printf("%d %f\n", checkin[0], d[0]);

    if (thread_idx < nint[0]) {
        for (auto i = 0; i < 3; i++) {
            coord[3 * thread_idx + i] = coord[3 * thread_idx + i] + d[0] * f[3 * thread_idx + i] ;
        }
    }
}
__global__ void shiftGradientCUDA(floatingpoint *f, floatingpoint* fAux, int * nint, floatingpoint* newGrad, floatingpoint* prevGrad, floatingpoint*
curGrad, bool *Mstate) {
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;
    if(Mstate[0] == false) return;

    floatingpoint d  = fmax(0.0, (newGrad[0] - prevGrad[0]) / curGrad[0]);
#ifdef DETAILEDOUTPUT_BETA
    if(thread_idx == 0) {
        printf("beta CUDA %f\n", d);
        printf(" newGrad %f prevGrad %f, curGrad %f\n", newGrad[0], prevGrad[0], curGrad[0]);
    }
#endif
//    if(thread_idx == 0)
//        printf("Shift Gradient %f %f %f %f\n", d, newGrad[0], prevGrad[0],curGrad[0] );

    if (thread_idx < nint[0] && Mstate[0]) {
        for (auto i = 0; i < 3; i++) {
            f[3 * thread_idx + i] = fAux[3 * thread_idx + i] + d * f[3 * thread_idx + i];
//            if (fabs(f[3 * thread_idx + i]) == __longlong_as_floatingpoint(0x7ff0000000000000) //infinity
//                || f[3 * thread_idx + i] != f[3 * thread_idx + i]) {
//                printf("Force became infinite during gradient shift. \n Force %f Aux Force %f Beta %f\n NewGrad %f "
//                               "PrevGrad %f curGrad %f \n",
//                       f[3 * thread_idx + i], fAux[3 * thread_idx + i], d, newGrad[0], prevGrad[0],curGrad[0]);
//                assert(0);
//            }
        }
    }
}

__global__ void shiftGradientCUDAifsafe(floatingpoint *f, floatingpoint* fAux, int * nint, bool *Mstate, bool *Sstate) {
    if(Mstate[0] == false || Sstate[0] == false) return;//checks for Minimization state and Safe state.
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;
//    if(thread_idx == 0)
//        printf("shiftGradient safe \n");
    if (thread_idx < nint[0]) {
        for (auto i = 0; i < 3; i++) {
            f[3 * thread_idx + i] = fAux[3 * thread_idx + i];
//            if (fabs(f[3 * thread_idx + i]) == __longlong_as_floatingpoint(0x7ff0000000000000) //infinity
//                || f[3 * thread_idx + i] != f[3 * thread_idx + i]) {
//                printf("Force became infinite during SAFE gradient shift. \n Force %f Aux Force %f \n",
//                       f[3 * thread_idx + i], fAux[3 * thread_idx + i]);
//                assert(0);
//            }
        }
    }
}

__global__ void allFADotFCUDA(floatingpoint *f1, floatingpoint *f2, floatingpoint *g, int * nint) {
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;
    if (thread_idx < nint[0]) {
            g[thread_idx] = 0.0;
    }
    __syncthreads();
    if (thread_idx < nint[0]) {
        for (auto i = 0; i < 3; i++) {
           g[thread_idx] +=f1[3 * thread_idx +i] * f2[3 * thread_idx +i];
        }
//        printf("CUDA %d %f %f %f %f\n",thread_idx, g[thread_idx],f1[3 * thread_idx],f1[3 * thread_idx +1],f1[3 *
//               thread_idx +2]);
    }
}

//__global__ void maxFCUDA(floatingpoint *f, int * nint, floatingpoint *fmax) {
//
//    floatingpoint mag;
//    fmax[0]=0.0;
//    for(int i = 0; i < 3 * nint[0]; i++) {
//        mag = sqrt(f[i]*f[i]);
//        if(mag > fmax[0]) {fmax[0] = mag;}
//    }
////    id = id -id%3;
//    printf("Fmaxv1 %f \n", fmax[0]);
//}

__global__ void maxFCUDAred(floatingpoint *f, int * nint, floatingpoint *fmax) {

    int offset = 0;
    extern __shared__ floatingpoint s[];
    floatingpoint *c1 = s;
    floatingpoint temp = -1.0;
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;
    while(thread_idx + offset < 3 * nint[0]){
        if(fabs(f[thread_idx + offset]) > temp)
            temp = fabs(f[thread_idx + offset]);
        offset += blockDim.x; //3
    }
    c1[threadIdx.x] = temp;
    __syncthreads();
    if(threadIdx.x == 0){
        fmax[0] = 0.0;
        for(auto i=0; i <3; i++){
            if(c1[i] > fmax[0])
                fmax[0] = c1[i];
        }
        printf("Fmax %f %f %f %f \n", fmax[0], c1[0], c1[1], c1[2]);
    }
}

__global__ void maxFCUDAredv2(floatingpoint *f, int * nint, floatingpoint *fmax) {

    extern __shared__ floatingpoint s[];
    floatingpoint *c1 = s;
    int start = 0;
    int end = nint[0];
    int factor = nint[0]/blockDim.x;
//    if(threadIdx.x ==0)
//        printf("CMAXF %d %d %d\n", factor, nint[0],blockDim.x);
    if (threadIdx.x > 0)
        start = threadIdx.x * factor;
    if (threadIdx.x < blockDim.x - 1)
        end = (threadIdx.x + 1) * factor;
    c1[threadIdx.x] = -1.0;
    for (auto i = start; i < end; i++) {
        if(fabs(sqrt(f[i])) > c1[threadIdx.x])
            c1[threadIdx.x] = sqrt(fabs(f[i]));
    }
    __syncthreads();

    if (threadIdx.x == 0) {
        fmax[0] = -1.0;
        for (auto i = 0; i < blockDim.x; i++) {
            if(fabs(c1[i]) > fmax[0])
                fmax[0] = c1[i];
        }
//        printf("Fmaxv2 %f f_x %f f_y %f f_z %f \n", fmax[0], c1[0], c1[1], c1[2]);
        printf("Fmax CUDA %f \n", fmax[0]);
    }
}
__global__ void maxFCUDAredv3(floatingpoint *g_idata, int *num, floatingpoint *F_max, int *mutex)
{
    extern __shared__ floatingpoint sdata[];
    unsigned int tid = threadIdx.x;
    auto blockSize = blockDim.x;
    unsigned int i = blockIdx.x*(blockSize*2) + tid;
    unsigned int gridSize = blockSize*2*gridDim.x;
    int n = num[0];
//    if(tid == 0)
//        printf("%d %d\n",num[0],num[1]);
    float ip[3];
    sdata[tid] = 0.0;
    while (i < n) {
//        if(g_idata[i] == -1.0 || g_idata[i+blockSize] == -1.0)
//        {printf("CUDA addition of energies. Energy is infinite\n");assert(0);}
        sdata[tid] = max(max(g_idata[i] , g_idata[i+blockSize]),sdata[tid]);
//        printf("blockID %d sdata[tid] %f Fmax %f \n", blockIdx.x, g_idata[i], F_max[0]);
        i += gridSize;
    }
    __syncthreads();
    if(blockSize >=2048 ) {printf("Cannot handle blocks with threads larger than 2048\n");assert(0);}
    if (blockSize >= 1024) { if (tid < 512) {
//            printf("blockID %d sdata[tid] %f Fmax %f \n", blockIdx.x, sdata[tid], F_max[0]);
            sdata[tid] = max(sdata[tid],sdata[tid + 512]); }__syncthreads(); }
    if (blockSize >= 512) { if (tid < 256) { sdata[tid] = max(sdata[tid],sdata[tid + 256]); } __syncthreads(); }
    if (blockSize >= 256) { if (tid < 128) { sdata[tid] = max(sdata[tid],sdata[tid + 128]); } __syncthreads(); }
    if (blockSize >= 128) { if (tid < 64) { sdata[tid] = max(sdata[tid],sdata[tid + 64]); }__syncthreads(); }
    if (tid < 32) {
        if (blockSize >= 64) {sdata[tid] = max(sdata[tid],sdata[tid + 32]);__syncthreads
                    (); }
        if (blockSize >= 32) {sdata[tid] = max(sdata[tid],sdata[tid + 16]);__syncthreads
                    (); }
        if (blockSize >= 16) {sdata[tid] = max(sdata[tid],sdata[tid + 8]);__syncthreads(); }
        if (blockSize >= 8) {sdata[tid] = max(sdata[tid],sdata[tid + 4]);__syncthreads(); }
        if (blockSize >= 4) {sdata[tid] = max(sdata[tid],sdata[tid + 2]);__syncthreads(); }
        if (blockSize >= 2) {sdata[tid] = max(sdata[tid],sdata[tid + 1]);__syncthreads(); }
    }
    if (tid == 0) {
        while(atomicCAS(mutex,0,1) != 0);  //lock
        F_max[0] = max(sqrt(sdata[0]), F_max[0]);
        atomicExch(mutex, 0);  //unlock
#ifdef DETAILEDOUTPUT_LAMBDA
            printf("Fmax CUDA %f \n", F_max[0]);
#endif
//        atomicAdd(&F_max[0],-F_max[0] + max(sqrt(sdata[0]), F_max[0]));
//        printf("blockID %d sdata[0] %f Fmax %f \n", blockIdx.x, sqrt(sdata[0]), F_max[0]);
    }
}

//__global__ void maxFCUDAredv2(floatingpoint *f, int * nint, floatingpoint *fmax) {
//
//    extern __shared__ floatingpoint s[];
//    floatingpoint *c1 = s;
//    int start = 0;
//    int end = 3 * nint[0];
//    int factor = 3 * nint[0]/blockDim.x;
//    if(threadIdx.x ==0)
//        printf("CMAXF %d %d %d\n", factor, nint[0],blockDim.x);
//    if (threadIdx.x > 0)
//        start = threadIdx.x * factor;
//    if (threadIdx.x < blockDim.x - 1)
//        end = (threadIdx.x + 1) * factor;
//    c1[threadIdx.x] = -1.0;
//    for (auto i = start; i < end; i++) {
//        if(fabs(f[i]) > c1[threadIdx.x])
//            c1[threadIdx.x] = fabs(f[i]);
//    }
//    __syncthreads();
//
//    if (threadIdx.x == 0) {
//        fmax[0] = -1.0;
//        for (auto i = 0; i < blockDim.x; i++) {
//            if(fabs(c1[i]) > fmax[0])
//                fmax[0] = c1[i];
//        }
////        printf("Fmaxv2 %f f_x %f f_y %f f_z %f \n", fmax[0], c1[0], c1[1], c1[2]);
//        printf("Fmax CUDA %f \n", fmax[0]);
//    }
//}

//__global__ void addvector(floatingpoint *U, int *params, floatingpoint *U_sum){
//    U_sum[0] = 0.0;
////    printf("%d \n", params[0]);
//
//    for(auto i=0;i<params[0];i++){
//        U_sum[0]  += U[i];
//    }
//    printf("add1 %f \n", U_sum[0]);
//}

//__global__ void addvectorred(floatingpoint *U, int *params, floatingpoint *U_sum){
//    extern __shared__ floatingpoint s[];
//    floatingpoint *c1 = s;
//    int start = 0;
//    int end = params[0];
//    int factor = params[0]/blockDim.x;
//    if(threadIdx.x > 0)
//        start = threadIdx.x * factor;
//    if(threadIdx.x < blockDim.x - 1)
//        end = (threadIdx.x + 1) * factor;
//    c1[threadIdx.x] = 0.0;
//    for(auto i = start; i < end; i++)
//        c1[threadIdx.x] += U[i];
////    printf("%d \n", params[0]);
//    __syncthreads();
//
//    if(threadIdx.x == 0) {
//        U_sum[0] = 0.0;
//        for (auto i = 0; i < blockDim.x; i++) {
//            U_sum[0] += c1[i];
//        }
////        printf("add2 %f \n", U_sum[0]);
//    }
//}

__global__ void addvectorredcgm(floatingpoint *g_idata, int *num, floatingpoint *g_odata) {
    extern __shared__ floatingpoint sdata[];
    unsigned int tid = threadIdx.x;
    auto blockSize = blockDim.x;
    unsigned int i = blockIdx.x*(blockSize*2) + tid;
    unsigned int gridSize = blockSize*2*gridDim.x;
    int n = num[0];
    sdata[tid] = 0;
    while (i < n) {
        sdata[tid] += g_idata[i] + g_idata[i+blockSize];
        i += gridSize;
    }
    __syncthreads();
    if(blockSize >=2048 ) {printf("Cannot handle blocks with threads larger than 2048\n");assert(0);}
    if (blockSize >= 1024) { if (tid < 512) { sdata[tid] += sdata[tid + 512]; } __syncthreads(); }
    if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
    if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
    if (blockSize >= 128) { if (tid < 64) { sdata[tid] += sdata[tid + 64]; } __syncthreads(); }
    if (tid < 32) {
        if (blockSize >= 64) {sdata[tid] += sdata[tid + 32];__syncthreads(); }
        if (blockSize >= 32) {sdata[tid] += sdata[tid + 16];__syncthreads(); }
        if (blockSize >= 16) {sdata[tid] += sdata[tid + 8];__syncthreads(); }
        if (blockSize >= 8) {sdata[tid] += sdata[tid + 4];__syncthreads(); }
        if (blockSize >= 4) {sdata[tid] += sdata[tid + 2];__syncthreads(); }
        if (blockSize >= 2) {sdata[tid] += sdata[tid + 1];__syncthreads(); }
    }
    if (tid == 0) {atomicAdd(&g_odata[0], sdata[0]);}
}


__global__ void initializeLambdaCUDA(bool *checkin, bool* checkout, floatingpoint *currentEnergy, floatingpoint *energy,
                                     floatingpoint* CUDA_lambda, floatingpoint* fmax, floatingpoint* params, bool *Safestate, int *status){
    checkin[0] = false;
    checkout[0] = false;
    status[0] = 0;
//    printf("lambda_status %d\n", status[0]);
    floatingpoint LAMBDAMAX = params[3];
//    printf("SS %d \n",Safestate[0]);
    if(Safestate[0] == true){//safebacktrackinglinesearch
        CUDA_lambda[0] = LAMBDAMAX;
    }
    else{//backtrackinglinesearch
        floatingpoint MAXDIST = params[4];

        if(fmax[0]==0.0) {
            CUDA_lambda[0] = 0.0;
            checkout[0] = true;
        }
        else
            CUDA_lambda[0] = min(LAMBDAMAX, MAXDIST / fmax[0]);
    }

#ifdef DETAILEDOUTPUT_LAMBDA
    printf("CUDA init lambda %.10e LAMBDAMAX %f MAXDIST %f fmax %f\n",
           CUDA_lambda[0], LAMBDAMAX, params[4], fmax[0]);
#endif
}

__global__ void resetlambdaCUDA (floatingpoint *CUDA_lambda){
    CUDA_lambda[0] = 0.0;
}
//__global__ void prepbacktracking(bool *checkin, bool* checkout, floatingpoint *currentEnergy, floatingpoint *energy,
//                                 floatingpoint* CUDA_lambda, floatingpoint* fmax, floatingpoint* params){
//    //if(checkin[0]) return;
//    checkin[0] = false;
//    checkout[0] = false;
//    currentEnergy[0] = energy[0];
//    checkout[0] = false;
//    floatingpoint LAMBDAMAX = params[3];
//    floatingpoint MAXDIST = params[4];
//
//    if(fmax[0]==0.0) {
//       CUDA_lambda[0] = 0.0;
//        checkout[0] = true;
//    }
//    else
//        CUDA_lambda[0] = fmin(LAMBDAMAX, MAXDIST / fmax[0]);
//    //printf("%f \n", CUDA_lambda[0]);
//}
//
//__global__ void prepsafebacktracking(bool* checkin, bool *checkout, floatingpoint* currentEnergy, floatingpoint* energy, floatingpoint*
//CUDA_lambda, floatingpoint* params) {
//    checkin[0] = false;
//    checkout[0] = false;
//    currentEnergy[0] = energy[0];
//    floatingpoint LAMBDAMAX = params[3];
//    CUDA_lambda[0] = LAMBDAMAX;
//}
__global__ void setcurrentenergy( floatingpoint* energy, floatingpoint* currentenergy, floatingpoint *CUDAlambda, floatingpoint *initlambdalocal){
    currentenergy[0] = energy[0];
    CUDAlambda[0] = initlambdalocal[0];
}

__global__ void findLambdaCUDA(floatingpoint* energyLambda, floatingpoint* currentEnergy, floatingpoint* FDotFA, floatingpoint *fmax, floatingpoint*
lambda, floatingpoint *params, bool* prev_convergence, bool*  current_convergence, bool *safestate, int *status){
#ifdef DETAILEDOUTPUT_LAMBDA
    printf("FCL 1 prev_conv %d curr_conv %d \n", prev_convergence[0],
           current_convergence[0]);
#endif
    if(prev_convergence[0]) return;
//    current_convergence[0] = false;
//    floatingpoint LAMBDAREDUCE = params[1];
//    floatingpoint LAMBDATOL = params[2];
//    floatingpoint MAXDIST = params[4];
    floatingpoint idealEnergyChange = 0.0;
    floatingpoint energyChange = 0.0;

#ifdef DETAILEDOUTPUT_LAMBDA
    printf("CUDA safestate %d\n",safestate[0]);
#endif
    if(safestate[0] == true){
        energyChange = energyLambda[0] - currentEnergy[0];

#ifdef DETAILEDOUTPUT_LAMBDA
        printf("CUDA energyChange %f currentEnergy %f energyLambda %f\n", energyChange,
               currentEnergy[0], energyLambda[0]);
#endif
        if (energyChange <= 0.0) {
            current_convergence[0] = true;
            atomicAdd(&status[0], 1);

#ifdef DETAILEDOUTPUT_LAMBDA
            printf("CUDA energyChange %f lambda_converged %.8f\n", energyChange, lambda[0]);
#endif
            return;
        }
    }
    else {
        floatingpoint BACKTRACKSLOPE = params[0];
        idealEnergyChange = -BACKTRACKSLOPE * lambda[0] * FDotFA[0];
        energyChange = energyLambda[0] - currentEnergy[0];

#ifdef DETAILEDOUTPUT_LAMBDA
        printf("BACKTRACKSLOPE %f lambda %.10e allFDotFA %f \n", BACKTRACKSLOPE, lambda[0],
               FDotFA[0]);
        printf("idealEnergyChange %f energyChange %f \n",
               idealEnergyChange, energyChange);
#endif
        if (energyChange <= idealEnergyChange) {
            current_convergence[0] = true;
//            printf("lambda_statusb %d\n", status[0]);
            atomicAdd(&status[0], 1);
//            printf("lambda_status %d\n", status[0]);
#ifdef DETAILEDOUTPUT_LAMBDA
            printf("idealEnergyChange %f energyChange %f lambda_converged %.8f \n",
                   idealEnergyChange, energyChange,
                   lambda[0]);
#endif
            return;
        }
    }

    //Reduce lambda if not converged.
    floatingpoint LAMBDAREDUCE = params[1];
    floatingpoint LAMBDATOL = params[2];
    floatingpoint MAXDIST = params[4];
    lambda[0] *= LAMBDAREDUCE;
    if(lambda[0] <= 0.0 || lambda[0] <= LAMBDATOL) {
        current_convergence[0] = true;
        atomicAdd(&status[0], 2);
        if(safestate[0] == true)
            lambda[0] = MAXDIST / fmax[0];
        else
            lambda[0] = 0.0;

    }
#ifdef DETAILEDOUTPUT_LAMBDA
    printf("lambda2 %.8f state %d\n", lambda[0], current_convergence[0]);
#endif
}
__global__ void findLambdaCUDA2(floatingpoint *fmax, floatingpoint* lambda, floatingpoint *params, bool* prev_convergence, bool*
current_convergence, bool *safestate, int *status){

#ifdef DETAILEDOUTPUT_LAMBDA
    printf("FCL 1 prev_conv %d curr_conv %d \n", prev_convergence[0],
           current_convergence[0]);
    printf("prev_conv %d \n", prev_convergence[0]);
#endif
    if(prev_convergence[0]) return;
    floatingpoint LAMBDAREDUCE = params[1];
    floatingpoint MAXDIST = params[4];
    floatingpoint LAMBDATOL = params[2];
    lambda[0] *= LAMBDAREDUCE;
    if(lambda[0] <= 0.0 || lambda[0] <= LAMBDATOL) {
        current_convergence[0] = true;
        atomicAdd(&status[0], 2);
        if(safestate[0] == true)
            lambda[0] = MAXDIST / fmax[0];
        else
            lambda[0] = 0.0;

    }
#ifdef DETAILEDOUTPUT_LAMBDA
    printf("lambda2 %.8f state %d\n", lambda[0], current_convergence[0]);
#endif
}
//__global__ void CUDAbacktrackingfindlambda(floatingpoint* energyLambda, floatingpoint* currentEnergy, floatingpoint* FDotFA, floatingpoint*
//lambda, floatingpoint *params, bool* prev_convergence, bool*  current_convergence){
////    printf("%d %d %f %f\n", prev_convergence[0],current_convergence[0],energyLambda[0],currentEnergy[0]);
//    if(prev_convergence[0]) return;
//    current_convergence[0] = false;
//    floatingpoint BACKTRACKSLOPE = params[0];
//    floatingpoint LAMBDAREDUCE = params[1];
//    floatingpoint LAMBDATOL = params[2];
//    floatingpoint idealEnergyChange = -BACKTRACKSLOPE * lambda[0] * FDotFA[0];
//    floatingpoint energyChange =  energyLambda[0] - currentEnergy[0];
////    printf("%f \n", lambda[0]);
//    if(energyChange <= idealEnergyChange) {
//        current_convergence[0] = true;
//        return;}
//    lambda[0] *= LAMBDAREDUCE;
////    printf("lambdareduce %f \n", lambda[0]);
//    if(lambda[0] <= 0.0 || lambda[0] <= LAMBDATOL) {
//        current_convergence[0] = true;
//        lambda[0] = 0.0;
//    }
//}

//__global__ void CUDAsafebacktrackingfindlambda(floatingpoint* energyLambda, floatingpoint* currentenergy, floatingpoint* fmax,
//                                               floatingpoint* lambda, floatingpoint* params,  bool* prev_convergence,
//                                               bool*  current_convergence){
//    if(prev_convergence[0]) return;
//    current_convergence[0] = false;
//    floatingpoint energyChange = energyLambda[0] - currentenergy[0];
//    floatingpoint LAMBDAREDUCE = params[1];
//    floatingpoint LAMBDATOL = params[2];
//    floatingpoint MAXDIST = params[4];
//    //return if ok
//    if(energyChange <= 0.0) {current_convergence[0] = true; return;}
//    //reduce lambda
//    lambda[0] *= LAMBDAREDUCE;
//
//    //just shake if we cant find an energy min,
//    //so we dont get stuck
//    if(lambda[0] <= 0.0 || lambda[0] <= LAMBDATOL) {current_convergence[0] = true; lambda[0] = MAXDIST / fmax[0];  }
//}


__global__ void getsafestateCUDA(floatingpoint* FDotFA, floatingpoint* curGrad, floatingpoint* newGrad, bool* checkout){
    //if(checkin[0] == true) return;
    checkout[0] = false;
    if(FDotFA[0] <= 0.0 || areEqual(curGrad[0], newGrad[0]))
        checkout[0] = true;
#ifdef DETAILEDOUTPUT_LAMBDA
    printf("CUDA FDotFA %f newGrad %f curGrad %f\n", FDotFA[0], newGrad[0], curGrad[0]);
    printf("CUDA safe state %d \n", checkout[0]);
#endif
    curGrad[0] = newGrad[0];
}
__global__ void getminimizestateCUDA(floatingpoint *fmax, floatingpoint *GRADTOL, bool *checkin, bool *checkout) {
    //maxF() > GRADTOL
#ifdef DETAILEDOUTPUT_ENERGY
    printf("maxF CUDA %f \n", fmax[0]);
#endif
    checkout[0] = false;
    if(checkin[0] == false) return;
    if(fmax[0] > GRADTOL[0])
        checkout[0] = true;
//    printf("minstate %f %f %d %d\n", fmax[0],GRADTOL[0],checkin[0],checkout[0]);
}

__global__ void initializePolak(bool* Mcheckin, bool *Mcheckout, bool *Scheckin, bool *Scheckout){
    Mcheckin[0] = true;
    Mcheckout[0] = true;
    Scheckin[0] = false;
    Scheckout[0] = false;
}
#endif
#endif //CUDA_VEC_CGMETHODCUDA_H

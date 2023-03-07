//
// Created by aravind on 11/1/17.
//

#ifndef CUDA_VEC_FORCEFIELDMANAGERCUDA_H
#define CUDA_VEC_FORCEFIELDMANAGERCUDA_H
#ifdef CUDAACCL
#include <cuda.h>
#include <cuda_runtime.h>
#include "ForceField.h"
#include "MathFunctions.h"
using namespace mathfunc;

__global__ void setenergytozero(floatingpoint *U_tot){
    U_tot[0] = 0.0;
}
__global__ void copyForcesCUDA(floatingpoint *f, floatingpoint *fAux, int* n){

    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

    if (thread_idx < n[0]) {
        for (auto i = 0; i < 3; i++) {
            fAux[3 * thread_idx + i] = f[3 * thread_idx + i];
        }
    }
}

__global__ void resetForcesCUDA(floatingpoint *f, int* n){

    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;
//    printf("%d \n", n[0]);
    if (thread_idx < n[0]) {
        for (auto i = 0; i < 3; i++) {

//            atomicAdd(&f[3 * thread_idx + i], -f[3 * thread_idx + i]);
//            atomicExch(&f[3 * thread_idx + i], 0.0);
            f[3 * thread_idx + i] = 0.0;
        }
    }
}

#endif
#endif //CUDA_VEC_FORCEFIELDMANAGERCUDA_H

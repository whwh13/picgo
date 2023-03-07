//
// Created by aravind on 1/25/18.
//

#ifndef CUDA_VEC_BRANCHINGSTRETCHINGHARMONICCUDA_H
#define CUDA_VEC_BRANCHINGSTRETCHINGHARMONICCUDA_H
#ifdef CUDAACCL
#include "BranchingStretchingHarmonic.h"

#include "Bead.h"

#include "MathFunctions.h"
#include "SysParams.h"
#include <limits>
#include <assert.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cmath>
#include <math.h>

using namespace mathfunc;

//#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600
//
//#else
//static __inline__ __device__ floatingpoint atomicAdd(floatingpoint *address, floatingpoint val) {
//    unsigned long long int* address_as_ull = (unsigned long long int*)address;
//    unsigned long long int old = *address_as_ull, assumed;
//    if (val==0.0)
//      return __longlong_as_floatingpoint(old);
//    do {
//      assumed = old;
//      old = atomicCAS(address_as_ull, assumed, __floatingpoint_as_longlong(val +__longlong_as_floatingpoint(assumed)));
//    } while (assumed != old);
//    return __longlong_as_floatingpoint(old);
//  }
//
//
//#endif

//__global__ void addvectorBS(floatingpoint *U, int *params, floatingpoint *U_sum, floatingpoint *U_tot){
//    U_sum[0] = 0.0;
//    floatingpoint sum = 0.0;
//    for(auto i=0;i<params[1];i++){
//        if(U[i] == -1.0 && sum != -1.0){
//            U_sum[0] = -1.0;
//            U_tot[0] = -1.0;
//            sum = -1.0;
//            break;
//        }
//        else
//            sum  += U[i];
//    }
//    U_sum[0] = sum;
//    atomicAdd(&U_tot[0], sum);
//
//}

__global__ void BranchingStretchingHarmonicenergy(floatingpoint *coord, floatingpoint *force, int *beadSet, floatingpoint *kstr,
                                                 floatingpoint *eql, floatingpoint *pos, int *params, floatingpoint *U_i, floatingpoint *z,  int
                                                  *culpritID,
                                                  char* culpritFF, char* culpritinteraction, char* FF, char*
interaction) {
    if(z[0] == 0.0) {
        extern __shared__ floatingpoint s[];
        floatingpoint *c1 = s;
        floatingpoint *c2 = &c1[3 * blockDim.x];
        floatingpoint *c3 = &c2[3 * blockDim.x];
        floatingpoint dist;
        floatingpoint v1[3];
        int nint = params[1];
        int n = params[0];
        const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

        if (thread_idx < nint) {
            for (auto i = 0; i < 3; i++) {
                U_i[thread_idx] = 0.0;
                c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
                c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
                c3[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 2] + i];
            }

        }
        __syncthreads();
        if (thread_idx < nint) {
            midPointCoordinate(v1, c1, c2, pos[thread_idx], 3 * threadIdx.x);
            dist = twoPointDistancemixedID(v1, c3, 0, 3 * threadIdx.x) - eql[thread_idx];

            U_i[thread_idx] = 0.5 * kstr[thread_idx] * dist * dist;

            if (fabs(U_i[thread_idx]) == __longlong_as_floatingpoint(0x7ff0000000000000) //infinity
                || U_i[thread_idx] != U_i[thread_idx] || U_i[thread_idx] < -1.0) {

                U_i[thread_idx] = -1.0;
                culpritID[0] = thread_idx;
                culpritID[1] = -1;
                int j = 0;
                while (FF[j] != 0) {
                    culpritFF[j] = FF[j];
                    j++;
                }
                j = 0;
                while (interaction[j] != 0) {
                    culpritinteraction[j] = interaction[j];
                    j++;
                }
                assert(0);
                __syncthreads();
            }
        }
    }
}

__global__ void BranchingStretchingHarmonicenergyz(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *kstr,
                                                  floatingpoint *eql, floatingpoint *pos, int *params, floatingpoint *U_i, floatingpoint *z,
                                                   int *culpritID, char* culpritFF, char* culpritinteraction, char* FF,
                                                   char* interaction) {
    if(z[0] != 0.0) {
        extern __shared__ floatingpoint s[];
        floatingpoint *c1 = s;
        floatingpoint *c2 = &c1[3 * blockDim.x];
        floatingpoint *c3 = &c2[3 * blockDim.x];
        floatingpoint *f1 = &c3[3 * blockDim.x];
        floatingpoint *f2 = &f1[3 * blockDim.x];
        floatingpoint *f3 = &f2[3 * blockDim.x];

        floatingpoint dist;
        floatingpoint v1[3];
        floatingpoint vzero[3];
        vzero[0] = 0;
        vzero[1] = 0;
        vzero[2] = 0;
        int nint = params[1];
        int n = params[0];
        const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

        if (thread_idx < nint) {
            U_i[thread_idx] = 0.0;
            for (auto i = 0; i < 3; i++) {
                c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
                c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
                c3[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 2] + i];
                f1[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx] + i];
                f2[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx + 1] + i];
                f3[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx + 2] + i];
            }

        }
        __syncthreads();

        if (thread_idx < nint) {

            midPointCoordinateStretched(v1, c1, f1, c2, f2, pos[thread_idx], z[0], 3 * threadIdx.x);
            dist = twoPointDistanceStretchedmixedID(v1, vzero, c3, f3, z[0], 0, 3 * threadIdx.x) - eql[thread_idx];

            U_i[thread_idx] = 0.5 * kstr[thread_idx] * dist * dist;

            if (fabs(U_i[thread_idx]) == __longlong_as_floatingpoint(0x7ff0000000000000) //infinity
                || U_i[thread_idx] != U_i[thread_idx] || U_i[thread_idx] < -1.0) {

                U_i[thread_idx] = -1.0;
                culpritID[0] = thread_idx;
                culpritID[1] = -1;
                int j = 0;
                while (FF[j] != 0) {
                    culpritFF[j] = FF[j];
                    j++;
                }
                j = 0;
                while (interaction[j] != 0) {
                    culpritinteraction[j] = interaction[j];
                    j++;
                }
                assert(0);
                __syncthreads();
            }

        }
    }
}


__global__ void BranchingStretchingHarmonicforces(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                                 floatingpoint *kstr, floatingpoint *eql, floatingpoint *pos, int *params){
    extern __shared__ floatingpoint s[];
    floatingpoint *c1 = s;
    floatingpoint *c2 = &c1[3 * blockDim.x];
    floatingpoint *c3 = &c2[3 * blockDim.x];
    floatingpoint dist, invL, f0;
    floatingpoint v1[3], f1[3], f2[3], f3[3];
    int nint = params[1];
    int n = params[0];
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

    if(thread_idx<nint) {
        for(auto i=0;i<3;i++){
            c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
            c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
            c3[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 2] + i];
        }

    }
    __syncthreads();

    if(thread_idx<nint) {
        midPointCoordinate(v1, c1, c2, pos[thread_idx], 3 * threadIdx.x);
        dist = twoPointDistancemixedID(v1, c3, 0, 3 * threadIdx.x);
        invL = 1 / dist;
        f0 = kstr[thread_idx] * ( dist - eql[thread_idx]) * invL;

        f1[0] =  -f0 * ( c3[3 * threadIdx.x] - v1[0] ) * (pos[thread_idx] - 1);
        f1[1] =  -f0 * ( c3[3 * threadIdx.x + 1] - v1[1] ) * (pos[thread_idx] - 1);
        f1[2] =  -f0 * ( c3[3 * threadIdx.x + 2] - v1[2] ) * (pos[thread_idx] - 1);

        // force i+1
        f2[0] =  f0 * ( c3[3 * threadIdx.x] - v1[0] ) * pos[thread_idx];
        f2[1] =  f0 * ( c3[3 * threadIdx.x + 1] - v1[1] ) * pos[thread_idx];
        f2[2] =  f0 * ( c3[3 * threadIdx.x + 2] - v1[2] ) * pos[thread_idx];

        //force on j
        f3[0] =  -f0 * ( c3[3 * threadIdx.x] - v1[0] );
        f3[1] =  -f0 * ( c3[3 * threadIdx.x + 1] - v1[1] );
        f3[2] =  -f0 * ( c3[3 * threadIdx.x + 2] - v1[2] );

        for (int i = 0; i < 3; i++) {
            atomicAdd(&f[3 * beadSet[n * thread_idx] + i], f1[i]);
            atomicAdd(&f[3 * beadSet[n * thread_idx + 1] + i], f2[i]);
            atomicAdd(&f[3 * beadSet[n * thread_idx + 2] + i], f3[i]);
        }
    }
}

#endif
#endif //CUDA_VEC_BRANCHINGSTRETCHINGHARMONICCUDA_H

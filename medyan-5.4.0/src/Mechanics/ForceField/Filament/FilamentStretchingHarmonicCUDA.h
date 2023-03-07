//
// Created by aravind on 1/22/18.
//

#ifndef CUDA_VEC_FILAMENTSTRETCHINGHARMONICCUDA_H
#define CUDA_VEC_FILAMENTSTRETCHINGHARMONICCUDA_H
#ifdef CUDAACCL
#include "FilamentStretchingHarmonic.h"

#include "FilamentStretching.h"

#include "Bead.h"

#include "MathFunctions.h"
#include "SysParams.h"
#include <limits>
#include <assert.h>
#include <cuda.h>
#include <cuda_runtime.h>

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

//__global__ void addvectorFS(floatingpoint *U, int *params, floatingpoint *U_sum, floatingpoint *U_tot){
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

__global__ void FilamentStretchingHarmonicenergy(floatingpoint *coord, floatingpoint *force, int *beadSet, floatingpoint *kstr,
                                                 floatingpoint *eql, int *params, floatingpoint *U_i, floatingpoint *z, int *culpritID,
                                                 char* culpritFF, char* culpritinteraction, char* FF, char*
                                                 interaction) {
    if(z[0] == 0.0) {
        extern __shared__ floatingpoint s[];
        floatingpoint *c1 = s;
        floatingpoint *c2 = &c1[3 * blockDim.x];
        floatingpoint dist;
        int nint = params[1];
        int n = params[0];
        const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

        if (thread_idx < nint) {
            U_i[thread_idx] = 0.0;
            for (auto i = 0; i < 3; i++) {
                U_i[thread_idx] = 0.0;
                c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
                c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
            }

            dist = twoPointDistance(c1, c2, 3 * threadIdx.x) - eql[thread_idx];
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

__global__ void FilamentStretchingHarmonicenergyz(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *kstr,
                                                  floatingpoint *eql, int *params, floatingpoint *U_i,
                                                  floatingpoint *U_vec, floatingpoint *z, int *culpritID,
                                                  char* culpritFF, char* culpritinteraction, char* FF, char*
                                                  interaction, bool*
                                                  conv_state1, bool* conv_state2){
if(conv_state1[0]||conv_state2[0]) {
//    if(blockIdx.x ==0 && threadIdx.x == 0)
//        printf("FStretch energy returning\n");
    return;
}
    if(z[0] == 0.0) {
//        extern __shared__ floatingpoint s[];
//        floatingpoint *c1 = s;
//        floatingpoint *c2 = &c1[3 * blockDim.x];
        floatingpoint dist;
        int nint = params[1];
        int n = params[0];
        int offset = max(params[2] -1,0);
        floatingpoint *c1, *c2;
        const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

        if (thread_idx < nint) {
//            if(thread_idx == 0){
//                printf("Offset %d \n", offset);
//            }
            U_i[thread_idx] = 0.0;
//            for (auto i = 0; i < 3; i++) {
//                U_i[thread_idx] = 0.0;
//                c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
//                c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
//            }

//            dist = twoPointDistance(c1, c2, 3 * threadIdx.x) - eql[thread_idx];
//
//            dist = twoPointDistancemixedID(coord, coord, 3 * beadSet[n * thread_idx], 3 *
//                                           beadSet[n * thread_idx + 1]) - eql[thread_idx];
            c1 = &coord[3 * beadSet[n * thread_idx]];
            c2 = &coord[3 * beadSet[n * thread_idx + 1]];
            dist = twoPointDistance(c1,c2) - eql[thread_idx];

            U_i[thread_idx] = 0.5 * kstr[thread_idx] * dist * dist;
            U_vec[offset + thread_idx] = 0.5 * kstr[thread_idx] * dist * dist;
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
    else if(z[0] != 0.0) {
//        extern __shared__ floatingpoint s[];
//        floatingpoint *c1 = s;
//        floatingpoint *c2 = &c1[3 * blockDim.x];
//        floatingpoint *f1 = &c2[3 * blockDim.x];
//        floatingpoint *f2 = &f1[3 * blockDim.x];
        floatingpoint *c1, *c2, *f1, *f2;
        floatingpoint dist;
        int nint = params[1];
        int n = params[0];
        int offset = max(params[2] -1,0);
        const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

        if (thread_idx < nint) {
//            if(thread_idx == 0){
//                printf("Offset %d \n", offset);
//            }
            U_i[thread_idx] = 0.0;
//            for (auto i = 0; i < 3; i++) {
//                c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
//                c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
//                f1[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx] + i];
//                f2[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx + 1] + i];
//            }
//           dist = twoPointDistanceStretched(c1, f1, c2, f2, z[0], 3 * threadIdx.x) - eql[thread_idx];

//            dist = twoPointDistanceStretchedmixedID(coord,coord,f,f, z[0], 3 * beadSet[n * thread_idx], 3 *
//                                                    beadSet[n * thread_idx + 1]) - eql[thread_idx];
            c1 = &coord[3 * beadSet[n * thread_idx]];
            c2 = &coord[3 * beadSet[n * thread_idx + 1]];
            f1 = &f[3 * beadSet[n * thread_idx]];
            f2 = &f[3 * beadSet[n * thread_idx + 1]];
            dist = twoPointDistanceStretched(c1,f1,c2,f2, z[0])- eql[thread_idx];
//            printf("dist %f \n", dist);
           U_i[thread_idx] = 0.5 * kstr[thread_idx] * dist * dist;
           U_vec[offset + thread_idx] = 0.5 * kstr[thread_idx] * dist * dist;
/*            if (fabs(U_i[thread_idx]) == __longlong_as_floatingpoint(0x7ff0000000000000) //infinity
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
            }*/

        }
    }
}


__global__ void FilamentStretchingHarmonicforces(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                                   floatingpoint *kstr, floatingpoint *eql, int *params){
    extern __shared__ floatingpoint s[];
    floatingpoint *c1 = s;
    floatingpoint *c2 = &c1[3 * blockDim.x];
//    floatingpoint *c1, *c2;
    floatingpoint dist, invL, f0, f1[3], f2[3];
    int nint = params[1];
    int n = params[0];
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

    if(thread_idx<nint) {
        for(auto i=0;i<3;i++){
            c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
            c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
        }

        dist = twoPointDistance(c1, c2, 3 * threadIdx.x);
//        c1 = &coord[3 * beadSet[n * thread_idx]];
//        c2 = &coord[3 * beadSet[n * thread_idx + 1]];
//        dist = twoPointDistance(c1,c2);
        invL = 1 / dist;
        f0 = kstr[thread_idx] * (dist - eql[thread_idx]) * invL;

        //force on i
        f2[0] = f0 * (c1[3 * threadIdx.x] - c2[3 * threadIdx.x]);
        f2[1] = f0 * (c1[3 * threadIdx.x + 1] - c2[3 * threadIdx.x + 1]);
        f2[2] = f0 * (c1[3 * threadIdx.x + 2] - c2[3 * threadIdx.x + 2]);

        // force i+1
        f1[0] = f0 * (-c1[3 * threadIdx.x] + c2[3 * threadIdx.x]);
        f1[1] = f0 * (-c1[3 * threadIdx.x + 1] + c2[3 * threadIdx.x + 1]);
        f1[2] = f0 * (-c1[3 * threadIdx.x + 2] + c2[3 * threadIdx.x + 2]);

        for (int i = 0; i < 3; i++) {
            if (fabs(f1[i]) == __longlong_as_floatingpoint(0x7ff0000000000000) //infinity
                || f1[i] != f1[i]) {
                printf("Fil. Stret. Force became infinite %f %f %f\n",f1[0], f1[1], f1[2]);
                assert(0);
            }
            if (fabs(f2[i]) == __longlong_as_floatingpoint(0x7ff0000000000000) //infinity
                || f2[i] != f2[i]) {
                printf("Fil. Stret. became infinite %f %f %f\n",f2[0], f2[1], f2[2]);
                assert(0);
            }
            atomicAdd(&f[3 * beadSet[n * thread_idx] + i], f1[i]);
            atomicAdd(&f[3 * beadSet[n * thread_idx + 1] + i], f2[i]);
        }

    }
}

#endif
#endif //CUDA_VEC_FILAMENTSTRETCHINGHARMONICCUDA_H

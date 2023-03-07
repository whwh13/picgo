//
// Created by aravind on 9/29/17.
//

#ifndef CUDA_VEC_MOTORGHOSTSTRETCHINGHARMONICCUDA_H
#define CUDA_VEC_MOTORGHOSTSTRETCHINGHARMONICCUDA_H
#ifdef CUDAACCL
#include "MotorGhostStretchingHarmonic.h"

#include "MotorGhostStretching.h"

#include "Bead.h"

#include "MathFunctions.h"
#include "SysParams.h"
#include <limits>
#include <assert.h>
#include <cuda.h>
#include <cuda_runtime.h>

using namespace mathfunc;


//
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

//__global__ void addvectorM(floatingpoint *U, int *params, floatingpoint *U_sum, floatingpoint *U_tot){
//  U_sum[0] = 0.0;
//    floatingpoint sum = 0.0;
//  for(auto i=0;i<params[1];i++){
//    if(U[i] == -1.0 && sum != -1.0){
//        U_sum[0] = -1.0;
//        U_tot[0] = -1.0;
//        sum = -1.0;
//      break;
//    }
//    else
//      sum  += U[i];
//  }
//    U_sum[0] = sum;
//    atomicAdd(&U_tot[0], sum);
//
//}

__global__ void MotorGhostStretchingHarmonicenergy(floatingpoint *coord, floatingpoint *force, int *beadSet, floatingpoint *kstr,
                                                   floatingpoint *eql, floatingpoint *pos1, floatingpoint *pos2, int *params,
                                                   floatingpoint *U_i, floatingpoint *z, int *culpritID,
                                                   char* culpritFF, char* culpritinteraction, char* FF, char*
                                                   interaction) {
    if(z[0] == 0.0) {
        extern __shared__ floatingpoint s[];
        floatingpoint *c1 = s;
        floatingpoint *c2 = &c1[3 * blockDim.x];
        floatingpoint *c3 = &c2[3 * blockDim.x];
        floatingpoint *c4 = &c3[3 * blockDim.x];
        floatingpoint v1[3], v2[3];
        floatingpoint dist;
        int nint = params[1];
        int n = params[0];
        const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

//    floatingpoint *coord1, *coord2, *coord3, *coord4, dist, U_i;
//    floatingpoint *v1 = new floatingpoint[3];
//    floatingpoint *v2 = new floatingpoint[3];
        if (thread_idx < nint) {
            for (auto i = 0; i < 3; i++) {
                U_i[thread_idx] = 0.0;
                c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
                c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
                c3[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 2] + i];
                c4[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 3] + i];
            }

//        }
//        __syncthreads();
//        if (thread_idx < nint) {

            midPointCoordinate(v1, c1, c2, pos1[thread_idx], 3 * threadIdx.x);
//        printf("%f %f %f \n",v1[0],v1[1],v1[2]);
            midPointCoordinate(v2, c3, c4, pos2[thread_idx], 3 * threadIdx.x);
//        printf("%f %f %f \n",v2[0],v2[1],v2[2]);
            dist = twoPointDistance(v1, v2) - eql[thread_idx];
//        printf("%f \n",dist);
//        U_i[thread_idx] = dist;
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

__global__ void MotorGhostStretchingHarmonicenergyz(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *kstr,
                                                   floatingpoint *eql, floatingpoint *pos1, floatingpoint *pos2, int *params,
                                                   floatingpoint *U_i, floatingpoint *U_vec, floatingpoint *z,
                                                    int *culpritID, char* culpritFF, char* culpritinteraction, char* FF, char*
                                                    interaction, bool*
                                                    conv_state1, bool* conv_state2){
    if(conv_state1[0]||conv_state2[0]) return;
    if(z[0] == 0.0) {
//        extern __shared__ floatingpoint s[];
//        floatingpoint *c1 = s;
//        floatingpoint *c2 = &c1[3 * blockDim.x];
//        floatingpoint *c3 = &c2[3 * blockDim.x];
//        floatingpoint *c4 = &c3[3 * blockDim.x];
        floatingpoint *c1, *c2, *c3, *c4;
        floatingpoint v1[3], v2[3];
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
//                U_i[thread_idx] = 0.0;
//                c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
//                c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
//                c3[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 2] + i];
//                c4[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 3] + i];
//            }

//            midPointCoordinate(v1, c1, c2, pos1[thread_idx], 3 * threadIdx.x);
//            midPointCoordinate(v2, c3, c4, pos2[thread_idx], 3 * threadIdx.x);

            c1 = &coord[3 * beadSet[n * thread_idx]];
            c2 = &coord[3 * beadSet[n * thread_idx + 1]];
            c3 = &coord[3 * beadSet[n * thread_idx + 2]];
            c4 = &coord[3 * beadSet[n * thread_idx + 3]];
            midPointCoordinate(v1, c1, c2, pos1[thread_idx]);
            midPointCoordinate(v2, c3, c4, pos2[thread_idx]);

            dist = twoPointDistance(v1, v2) - eql[thread_idx];
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

    else if(z[0] != 0.0) {

//        extern __shared__ floatingpoint s[];
//        floatingpoint *c1 = s;
//        floatingpoint *c2 = &c1[3 * blockDim.x];
//        floatingpoint *c3 = &c2[3 * blockDim.x];
//        floatingpoint *c4 = &c3[3 * blockDim.x];
//        floatingpoint *f1 = &c4[3 * blockDim.x];
//        floatingpoint *f2 = &f1[3 * blockDim.x];
//        floatingpoint *f3 = &f2[3 * blockDim.x];
//        floatingpoint *f4 = &f3[3 * blockDim.x];
        floatingpoint *c1, *c2, *c3, *c4, *f1, *f2, *f3, *f4;

        floatingpoint v1[3], v2[3];
        floatingpoint dist;
        int nint = params[1];
        int n = params[0];
        int offset = max(params[2] -1 , 0 );
        const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

        if (thread_idx < nint) {
//            if(thread_idx == 0){
//                printf("Offset %d \n", offset);
//            }
            U_i[thread_idx] = 0.0;
//            for (auto i = 0; i < 3; i++) {
//                c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
//                c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
//                c3[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 2] + i];
//                c4[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 3] + i];
//                f1[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx] + i];
//                f2[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx + 1] + i];
//                f3[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx + 2] + i];
//                f4[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx + 3] + i];
//            }

//            midPointCoordinateStretched(v1, c1, f1, c2, f2, pos1[thread_idx], z[0], 3 * threadIdx.x);
//            midPointCoordinateStretched(v2, c3, f3, c4, f4, pos2[thread_idx], z[0], 3 * threadIdx.x);

            c1 = &coord[3 * beadSet[n * thread_idx]];
            c2 = &coord[3 * beadSet[n * thread_idx + 1]];
            c3 = &coord[3 * beadSet[n * thread_idx + 2]];
            c4 = &coord[3 * beadSet[n * thread_idx + 3]];
            f1 = &f[3 * beadSet[n * thread_idx]];
            f2 = &f[3 * beadSet[n * thread_idx + 1]];
            f3 = &f[3 * beadSet[n * thread_idx + 2]];
            f4 = &f[3 * beadSet[n * thread_idx + 3]];
            midPointCoordinateStretched(v1, c1, f1, c2, f2, pos1[thread_idx], z[0]);
            midPointCoordinateStretched(v2, c3, f3, c4, f4, pos2[thread_idx], z[0]);

            dist = twoPointDistance(v1, v2) - eql[thread_idx];
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

}


__global__ void MotorGhostStretchingHarmonicforces(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                          floatingpoint *kstr, floatingpoint *eql, floatingpoint *pos1, floatingpoint
                                                   *pos2, int *params, floatingpoint
                                                   *MStretchingforce){

//    extern __shared__ floatingpoint s[];
//    floatingpoint *c1 = s;
//    floatingpoint *c2 = &c1[3 * blockDim.x];
//    floatingpoint *c3 = &c2[3 * blockDim.x];
//    floatingpoint *c4 = &c3[3 * blockDim.x];
    floatingpoint *c1, *c2, *c3, *c4;

    floatingpoint v1[3], v2[3];
    floatingpoint dist, invL, f0, f1[3], f2[3], f3[3], f4[3];
    int nint = params[1];
    int n = params[0];
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

    if(thread_idx<nint) {
//        for(auto i=0;i<3;i++){
//            c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
//            c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
//            c3[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 2] + i];
//            c4[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 3] + i];
//        }

//    }
//    __syncthreads();
//
//    if(thread_idx<nint) {

        c1 = &coord[3 * beadSet[n * thread_idx]];
        c2 = &coord[3 * beadSet[n * thread_idx + 1]];
        c3 = &coord[3 * beadSet[n * thread_idx + 2]];
        c4 = &coord[3 * beadSet[n * thread_idx + 3]];
        midPointCoordinate(v1, c1, c2, pos1[thread_idx]);
        midPointCoordinate(v2, c3, c4, pos2[thread_idx]);

//        midPointCoordinate(v1, c1, c2, pos1[thread_idx], 3 * threadIdx.x);
//        midPointCoordinate(v2, c3, c4, pos2[thread_idx], 3 * threadIdx.x);
        dist = twoPointDistance(v1, v2);
        invL = 1 / dist;
        f0 = kstr[thread_idx] * (dist - eql[thread_idx]) * invL;
        MStretchingforce[thread_idx] = f0/invL;

        //force on i
        f1[0] = -f0 * (v1[0] - v2[0]) * (1 - pos1[thread_idx]);
        f1[1] = -f0 * (v1[1] - v2[1]) * (1 - pos1[thread_idx]);
        f1[2] = -f0 * (v1[2] - v2[2]) * (1 - pos1[thread_idx]);

        // force i+1
        f2[0] = -f0 * (v1[0] - v2[0]) * (pos1[thread_idx]);
        f2[1] = -f0 * (v1[1] - v2[1]) * (pos1[thread_idx]);
        f2[2] = -f0 * (v1[2] - v2[2]) * (pos1[thread_idx]);

        //force on j
        f3[0] = f0 * (v1[0] - v2[0]) * (1 - pos2[thread_idx]);
        f3[1] = f0 * (v1[1] - v2[1]) * (1 - pos2[thread_idx]);
        f3[2] = f0 * (v1[2] - v2[2]) * (1 - pos2[thread_idx]);

        // force j+1
        f4[0] = f0 * (v1[0] - v2[0]) * (pos2[thread_idx]);
        f4[1] = f0 * (v1[1] - v2[1]) * (pos2[thread_idx]);
        f4[2] = f0 * (v1[2] - v2[2]) * (pos2[thread_idx]);

        for (int i = 0; i < 3; i++) {
            if (fabs(f1[i]) == __longlong_as_floatingpoint(0x7ff0000000000000) //infinity
                || f1[i] != f1[i]) {
                printf("Motor Force became infinite %f %f %f\n",f1[0], f1[1], f1[2]);
                assert(0);
            }
            if (fabs(f2[i]) == __longlong_as_floatingpoint(0x7ff0000000000000) //infinity
                || f2[i] != f2[i]) {
                printf("Motor Force became infinite %f %f %f\n",f2[0], f2[1], f2[2]);
                assert(0);
            }
            if (fabs(f3[i]) == __longlong_as_floatingpoint(0x7ff0000000000000) //infinity
                || f3[i] != f3[i]) {
                printf("Motor Force became infinite %f %f %f\n",f3[0], f3[1], f3[2]);
                assert(0);
            }
            if (fabs(f4[i]) == __longlong_as_floatingpoint(0x7ff0000000000000) //infinity
                || f4[i] != f4[i]) {
                printf("Motor Force became infinite %f %f %f\n",f4[0], f4[1], f4[2]);
                assert(0);
            }
            atomicAdd(&f[3 * beadSet[n * thread_idx] + i], f1[i]);
            atomicAdd(&f[3 * beadSet[n * thread_idx + 1] + i], f2[i]);
            atomicAdd(&f[3 * beadSet[n * thread_idx + 2] + i], f3[i]);
            atomicAdd(&f[3 * beadSet[n * thread_idx + 3] + i], f4[i]);
        }
    }
    }

#endif
#endif //CUDA_VEC_MOTORGHOSTSTRETCHINGHARMONICCUDA_H

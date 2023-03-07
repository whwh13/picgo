//
// Created by aravind on 1/25/18.
//

#ifndef CUDA_VEC_BRANCHINGPOSITIONCOSINECUDA_H
#define CUDA_VEC_BRANCHINGPOSITIONCOSINECUDA_H
#ifdef CUDAACCL
#include "BranchingPositionCosine.h"

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

//__global__ void addvectorBPC(floatingpoint *U, int *params, floatingpoint *U_sum, floatingpoint *U_tot){
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

__global__ void BranchingPositionCosineenergy(floatingpoint *coord, floatingpoint *force, int *beadSet, floatingpoint *kpos,
                                            floatingpoint *pos, int *params, floatingpoint *U_i, floatingpoint *z,  int *culpritID,
                                              char* culpritFF, char* culpritinteraction, char* FF, char*
                                            interaction) {
    if(z[0] == 0.0) {
        extern __shared__ floatingpoint s[];
        floatingpoint *c1 = s;
        floatingpoint *c2 = &c1[3 * blockDim.x];
        floatingpoint *c3 = &c2[3 * blockDim.x];
        floatingpoint X, D, XD, xd, theta, posheta, dTheta;
        floatingpoint mp[3];
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
            midPointCoordinate(mp, c1, c2, pos[thread_idx], 3 * threadIdx.x);
            X = sqrt(scalarProductmixedID(mp, c2, mp, c2, 0, 3 * threadIdx.x, 0, 3 * threadIdx.x));
            D = sqrt(scalarProductmixedID(mp, c3, mp, c3, 0, 3 * threadIdx.x, 0, 3 * threadIdx.x));

            XD = X * D;

            xd = (scalarProductmixedID(mp, c2, mp, c3, 0, 3 * threadIdx.x, 0, 3 * threadIdx.x));

            theta = safeacos(xd / XD);
            posheta = 0.5 * M_PI;
            dTheta = theta - posheta;

            U_i[thread_idx] = kpos[thread_idx] * (1 - cos(dTheta));

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
//    __syncthreads();
}

__global__ void BranchingPositionCosineenergyz(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *kpos,
                                             floatingpoint *pos, int *params, floatingpoint *U_i, floatingpoint *z,  int *culpritID,
                                               char* culpritFF, char* culpritinteraction, char* FF, char*
                                             interaction) {
    if(z[0] != 0.0) {
        extern __shared__ floatingpoint s[];
        floatingpoint *c1 = s;
        floatingpoint *c2 = &c1[3 * blockDim.x];
        floatingpoint *c3 = &c2[3 * blockDim.x];
        floatingpoint *f1 = &c3[3 * blockDim.x];
        floatingpoint *f2 = &f1[3 * blockDim.x];
        floatingpoint *f3 = &f2[3 * blockDim.x];
        floatingpoint X, D, XD, xd, theta, posheta, dTheta;
        floatingpoint mp[3];
        floatingpoint vzero[3];
        vzero[0] = 0.0;
        vzero[1] = 0.0;
        vzero[2] = 0.0;

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
            midPointCoordinateStretched(mp, c1, f1, c2, f2, pos[thread_idx], z[0], 3 * threadIdx.x);
            X = sqrt(scalarProductStretchedmixedID(mp, vzero, c2, f2, mp, vzero, c2, f2, z[0], 0, 3 * threadIdx
                    .x, 0, 3 * threadIdx.x));
            D = sqrt(scalarProductStretchedmixedID(mp, vzero, c3, f3, mp, vzero, c3, f3, z[0], 0, 3 * threadIdx
                    .x, 0, 3 * threadIdx.x));
            XD = X * D;
            xd = (scalarProductStretchedmixedID(mp, vzero, c2, f2, mp, vzero, c3, f3, z[0], 0, 3 * threadIdx
                    .x, 0, 3 * threadIdx.x));

            theta = safeacos(xd / XD);
            posheta = 0.5 * M_PI;
            dTheta = theta - posheta;

            U_i[thread_idx] = kpos[thread_idx] * (1 - cos(dTheta));

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

__global__ void BranchingPositionCosineforces(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                            floatingpoint *kpos, floatingpoint *pos, int *params){
    extern __shared__ floatingpoint s[];
    floatingpoint *c1 = s;
    floatingpoint *c2 = &c1[3 * blockDim.x];
    floatingpoint *c3 = &c2[3 * blockDim.x];
    floatingpoint X, D, XD, xd, invX, invD, position, A, B, C, k, theta, posheta, dTheta;
    floatingpoint mp[3], f1[3], f2[3], f3[3];

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
        midPointCoordinate(mp, c1, c2, pos[thread_idx], 3 * threadIdx.x);
        X = sqrt(scalarProductmixedID(mp, c2, mp, c2, 0, 3 * threadIdx.x, 0 , 3 * threadIdx.x));
        D = sqrt(scalarProductmixedID(mp, c3, mp, c3, 0, 3 * threadIdx.x, 0 , 3 * threadIdx.x));

        XD = X * D;

        xd = scalarProductmixedID(mp, c2, mp, c3, 0, 3 * threadIdx.x, 0 , 3 * threadIdx.x);

        invX = 1/X;
        invD = 1/D;
        A = invX*invD;
        B = invX*invX;
        C = invD*invD;

        theta = safeacos(xd / XD);
        posheta = 0.5*M_PI;
        dTheta = theta-posheta;

        position = pos[thread_idx];

        k =  kpos[thread_idx] * A * sin(dTheta)/sin(theta);

        //bead 1
        f1[0] =  k * (1-position)* (- (1-position)*(c2[3 * threadIdx.x] - c1[3 * threadIdx.x]) - (c3[3 * threadIdx.x] - (1-position)*c1[3 * threadIdx.x] - position*c2[3 * threadIdx.x])
                                     + xd *(B*(1-position)*(c2[3 * threadIdx.x] - c1[3 * threadIdx.x]) + C*(c3[3 * threadIdx.x] - (1-position)*c1[3 * threadIdx.x] - position*c2[3 * threadIdx.x])) );

        f1[1] =  k * (1-position)* (- (1-position)*(c2[3 * threadIdx.x + 1] - c1[3 * threadIdx.x + 1]) - (c3[3 * threadIdx.x + 1] - (1-position)*c1[3 * threadIdx.x + 1] - position*c2[3 * threadIdx.x + 1])
                                     + xd *(B*(1-position)*(c2[3 * threadIdx.x + 1] - c1[3 * threadIdx.x + 1]) + C*(c3[3 * threadIdx.x + 1] - (1-position)*c1[3 * threadIdx.x + 1] - position*c2[3 * threadIdx.x + 1])) );

        f1[2] =  k * (1-position)* (- (1-position)*(c2[3 * threadIdx.x + 2] - c1[3 * threadIdx.x + 2]) - (c3[3 * threadIdx.x + 2] - (1-position)*c1[3 * threadIdx.x + 2] - position*c2[3 * threadIdx.x + 2])
                                     + xd *(B*(1-position)*(c2[3 * threadIdx.x + 2] - c1[3 * threadIdx.x + 2]) + C*(c3[3 * threadIdx.x + 2] - (1-position)*c1[3 * threadIdx.x + 2] - position*c2[3 * threadIdx.x + 2])) );

        //bead 2

        f2[0] =  k * (- position*(1-position)*(c2[3 * threadIdx.x] - c1[3 * threadIdx.x]) + (1-position)*(c3[3 * threadIdx.x]- (1-position)*c1[3 * threadIdx.x] -
                position*c2[3 * threadIdx.x]) + xd *( (1-position)*B*(1-position)*(c2[3 * threadIdx.x] - c1[3 * threadIdx.x]) - position*C*(c3[3 * threadIdx.x] -
                (1-position)*c1[3 * threadIdx.x] - position*c2[3 * threadIdx.x])) );

        f2[1] =  k * (- position*(1-position)*(c2[3 * threadIdx.x + 1] - c1[3 * threadIdx.x + 1]) + (1-position)*(c3[3 * threadIdx.x + 1]- (1-position)*c1[3 * threadIdx.x + 1] - position*c2[3 * threadIdx.x + 1])
                       + xd *( (1-position)*B*(1-position)*(c2[3 * threadIdx.x + 1] - c1[3 * threadIdx.x + 1]) - position*C*(c3[3 * threadIdx.x + 1]
                       - (1-position)*c1[3 * threadIdx.x + 1] - position*c2[3 * threadIdx.x + 1])) );

        f2[2] =  k * (- position*(1-position)*(c2[3 * threadIdx.x + 2] - c1[3 * threadIdx.x + 2]) + (1-position)*(c3[3 * threadIdx.x + 2]- (1-position)*c1[3 * threadIdx.x + 2] - position*c2[3 * threadIdx.x + 2])
                       + xd *( (1-position)*B*(1-position)*(c2[3 * threadIdx.x + 2] - c1[3 * threadIdx.x + 2]) - position*C*(c3[3 * threadIdx.x + 2] - (1-position)*c1[3 * threadIdx.x + 2] -
                position*c2[3 * threadIdx.x + 2])) );

        //bead3

        f3[0] =  k * ( (1-position)*(c2[3 * threadIdx.x] - c1[3 * threadIdx.x]) - xd * C*(c3[3 * threadIdx.x] - (1-position)*c1[3 * threadIdx.x] - position*c2[3 * threadIdx.x]) );
        f3[1] =  k * ( (1-position)*(c2[3 * threadIdx.x + 1] - c1[3 * threadIdx.x + 1]) - xd * C*(c3[3 * threadIdx.x + 1] - (1-position)*c1[3 * threadIdx.x + 1] - position*c2[3 * threadIdx.x + 1]) );
        f3[2] =  k * ( (1-position)*(c2[3 * threadIdx.x + 2] - c1[3 * threadIdx.x + 2]) - xd * C*(c3[3 * threadIdx.x + 2] - (1-position)*c1[3 * threadIdx.x + 2] - position*c2[3 * threadIdx.x + 2]) );
        for (int i = 0; i < 3; i++) {
            atomicAdd(&f[3 * beadSet[n * thread_idx] + i], f1[i]);
            atomicAdd(&f[3 * beadSet[n * thread_idx + 1] + i], f2[i]);
            atomicAdd(&f[3 * beadSet[n * thread_idx + 2] + i], f3[i]);
        }
    }
}
#endif
#endif //CUDA_VEC_BRANCHINGPOSITIONCOSINECUDA_H

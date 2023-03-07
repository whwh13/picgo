//
// Created by aravind on 9/26/17.
//

#ifndef CUDA_VEC_CYLINDEREXCLVOLREPULSIONCUDA_H
#define CUDA_VEC_CYLINDEREXCLVOLREPULSIONCUDA_H
#ifdef CUDAACCL
#include "CylinderExclVolRepulsion.h"

#include "CylinderExclVolume.h"

#include "Bead.h"

#include "MathFunctions.h"
#include "SysParams.h"
#include <limits>

#include <cuda.h>
#include <cuda_runtime.h>
#include <assert.h>
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

__global__
void saxpy(int n, float a, float *x, float *y)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i < n) y[i] = a*x[i] + y[i];
}
__global__ void CUDAExclVolRepulsionenergy(floatingpoint *coord, floatingpoint *force, int *beadSet, floatingpoint *krep,
                                           int *params, floatingpoint *U_i, floatingpoint *z, int *culpritID,
                                           char* culpritFF, char* culpritinteraction, char* FField, char*
                                           interaction) {
//memory needed: 34*THREADSPERBLOCK*sizeof(floatingpoint)+2*THREADSPERBLOCK*sizeof(int);
    if(z[0] == 0.0) {
        extern __shared__ floatingpoint s[];
//    floatingpoint *coord_image=s;
//    floatingpoint *c1 = &coord_image[3 * blockDim.x];
        floatingpoint *c1 = s;
        floatingpoint *c2 = &c1[3 * blockDim.x];
        floatingpoint *c3 = &c2[3 * blockDim.x];
        floatingpoint *c4 = &c3[3 * blockDim.x];
        floatingpoint d, invDSquare;
        floatingpoint a, b, c, e, F, AA, BB, CC, DD, EE, FF, GG, HH, JJ;
        floatingpoint ATG1, ATG2, ATG3, ATG4;
        int nint = params[1];
        int n = params[0];
        const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

        if (thread_idx < nint) {
            U_i[thread_idx] = 0.0;
            for (auto i = 0; i < 3; i++) {
                c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
                c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
                c3[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 2] + i];
                c4[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 3] + i];
            }

//        }
//        __syncthreads();
//
//        if (thread_idx < nint) {
            //check if parallel
            if (areParallel(c1, c2, c3, c4, 3 * threadIdx.x)) {

                d = twoPointDistance(c1, c3, 3 * threadIdx.x);
                invDSquare = 1 / (d * d);
                U_i[thread_idx] = krep[thread_idx] * invDSquare * invDSquare;

                if (U_i[thread_idx] == __longlong_as_floatingpoint(0x7ff0000000000000) //infinity
                    || U_i[thread_idx] != U_i[thread_idx] || U_i[thread_idx] < -1.0) {
                    U_i[thread_idx] = -1.0;
                    culpritID[0] = beadSet[n * thread_idx];
                    culpritID[1] = beadSet[n * thread_idx + 1];
                    culpritID[2] = beadSet[n * thread_idx + 2];
                    culpritID[3] = beadSet[n * thread_idx + 3];
                    int j = 0;
                    while (FField[j] != 0) {
                        culpritFF[j] = FField[j];
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

            } else {

                //check if in same plane
                if (areInPlane(c1, c2, c3, c4, 3 * threadIdx.x)) {

                    //slightly move point
                    movePointOutOfPlane(c1, c2, c3, c4, 2, 0.01f, 3 * threadIdx.x);
                }

                a = scalarProduct(c1, c2, c1, c2, 3 * threadIdx.x);
                b = scalarProduct(c3, c4, c3, c4, 3 * threadIdx.x);
                c = scalarProduct(c3, c1, c3, c1, 3 * threadIdx.x);
                d = scalarProduct(c1, c2, c3, c4, 3 * threadIdx.x);
                e = scalarProduct(c1, c2, c3, c1, 3 * threadIdx.x);
                F = scalarProduct(c3, c4, c3, c1, 3 * threadIdx.x);

                AA = sqrt(a * c - e * e);
                BB = sqrt(b * c - F * F);

                CC = d * e - a * F;
                DD = b * e - d * F;

                EE = sqrt(a * (b + c - 2 * F) - (d - e) * (d - e));
                FF = sqrt(b * (a + c + 2 * e) - (d + F) * (d + F));

                GG = d * d - a * b - CC;
                HH = CC + GG - DD;
                JJ = c * (GG + CC) + e * DD - F * CC;


                ATG1 = atan((a + e) / AA) - atan(e / AA);
                ATG2 = atan((a + e - d) / EE) - atan((e - d) / EE);
                ATG3 = atan((F) / BB) - atan((F - b) / BB);
                ATG4 = atan((d + F) / FF) - atan((d + F - b) / FF);

                U_i[thread_idx] = 0.5 * krep[thread_idx] / JJ *
                                  (CC / AA * ATG1 + GG / EE * ATG2 + DD / BB * ATG3 + HH / FF * ATG4);

                if (U_i[thread_idx] == __longlong_as_floatingpoint(0x7ff0000000000000)
                    || U_i[thread_idx] != U_i[thread_idx] || U_i[thread_idx] < -1.0) {
                    U_i[thread_idx] = -1.0;
                    culpritID[0] = beadSet[n * thread_idx];
                    culpritID[1] = beadSet[n * thread_idx + 1];
                    culpritID[2] = beadSet[n * thread_idx + 2];
                    culpritID[3] = beadSet[n * thread_idx + 3];
                    int j = 0;
                    while (FField[j] != 0) {
                        culpritFF[j] = FField[j];
                        j++;
                    }
                    j = 0;
                    while (interaction[j] != 0) {
                        culpritinteraction[j] = interaction[j];
                        j++;
                    }
                    printf("Coordiantes \n %f %f %f \n %f %f %f \n %f %f %f \n %f %f %f \n", c1[3 * threadIdx.x ],
                           c1[3 * threadIdx
                    .x + 1], c1[3 * threadIdx.x + 2], c2[3 * threadIdx.x ], c2[3 * threadIdx
                            .x + 1], c2[3 * threadIdx.x + 2], c3[3 * threadIdx.x ], c3[3 * threadIdx
                            .x + 1], c3[3 * threadIdx.x + 2], c4[3 * threadIdx.x ], c4[3 * threadIdx
                            .x + 1], c4[3 * threadIdx.x + 2]);
                    assert(0);
                    __syncthreads();
                }
            }
        }
    }
}


__global__ void CUDAExclVolRepulsionenergyz(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                            floatingpoint *krep, int *params, floatingpoint *U_i,
                                            floatingpoint *U_vec, floatingpoint *z, int *culpritID,
                                            char* culpritFF, char* culpritinteraction, char* FField, char*
                                            interaction, bool*
                                            conv_state1, bool* conv_state2){
    if(conv_state1[0]||conv_state2[0]) return;
    if(z[0] == 0.0) {
        extern __shared__ floatingpoint s[];
//    floatingpoint *coord_image=s;
//    floatingpoint *c1 = &coord_image[3 * blockDim.x];
        floatingpoint *c1 = s;
        floatingpoint *c2 = &c1[3 * blockDim.x];
        floatingpoint *c3 = &c2[3 * blockDim.x];
        floatingpoint *c4 = &c3[3 * blockDim.x];
        floatingpoint d, invDSquare;
        floatingpoint a, b, c, e, F, AA, BB, CC, DD, EE, FF, GG, HH, JJ;
        floatingpoint ATG1, ATG2, ATG3, ATG4;
        int nint = params[1];
        int n = params[0];
        int offset = max(params[2] -1 , 0 );
        const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

        if (thread_idx < nint) {
//            if(thread_idx == 0){
//                printf("Offset %d, nint %d \n", offset, nint);
//            }
            U_i[thread_idx] = 0.0;
            for (auto i = 0; i < 3; i++) {
                c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
                c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
                c3[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 2] + i];
                c4[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 3] + i];
            }

            //check if parallel
            if (areParallel(c1, c2, c3, c4, 3 * threadIdx.x)) {

                d = twoPointDistance(c1, c3, 3 * threadIdx.x);
                invDSquare = 1 / (d * d);
                U_i[thread_idx] = krep[thread_idx] * invDSquare * invDSquare;
                U_vec[offset + thread_idx] = krep[thread_idx] * invDSquare * invDSquare;
//                printf("%d %f P\n", thread_idx, U_i[thread_idx]);
                if (U_i[thread_idx] == __longlong_as_floatingpoint(0x7ff0000000000000) //infinity
                    || U_i[thread_idx] != U_i[thread_idx] || U_i[thread_idx] < -1.0) {
                    U_i[thread_idx] = -1.0;
                    culpritID[0] = beadSet[n * thread_idx];
                    culpritID[1] = beadSet[n * thread_idx + 1];
                    culpritID[2] = beadSet[n * thread_idx + 2];
                    culpritID[3] = beadSet[n * thread_idx + 3];
                    int j = 0;
                    while (FField[j] != 0) {
                        culpritFF[j] = FField[j];
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

            } else {

                //check if in same plane
                if (areInPlane(c1, c2, c3, c4, 3 * threadIdx.x)) {

                    //slightly move point
                    movePointOutOfPlane(c1, c2, c3, c4, 2, 0.01f, 3 * threadIdx.x);
                }

                a = scalarProduct(c1, c2, c1, c2, 3 * threadIdx.x);
                b = scalarProduct(c3, c4, c3, c4, 3 * threadIdx.x);
                c = scalarProduct(c3, c1, c3, c1, 3 * threadIdx.x);
                d = scalarProduct(c1, c2, c3, c4, 3 * threadIdx.x);
                e = scalarProduct(c1, c2, c3, c1, 3 * threadIdx.x);
                F = scalarProduct(c3, c4, c3, c1, 3 * threadIdx.x);

                AA = sqrt(a * c - e * e);
                BB = sqrt(b * c - F * F);

                CC = d * e - a * F;
                DD = b * e - d * F;

                EE = sqrt(a * (b + c - 2 * F) - (d - e) * (d - e));
                FF = sqrt(b * (a + c + 2 * e) - (d + F) * (d + F));

                GG = d * d - a * b - CC;
                HH = CC + GG - DD;
                JJ = c * (GG + CC) + e * DD - F * CC;


                ATG1 = atan((a + e) / AA) - atan(e / AA);
                ATG2 = atan((a + e - d) / EE) - atan((e - d) / EE);
                ATG3 = atan((F) / BB) - atan((F - b) / BB);
                ATG4 = atan((d + F) / FF) - atan((d + F - b) / FF);

                U_i[thread_idx] = 0.5 * krep[thread_idx] / JJ *
                                  (CC / AA * ATG1 + GG / EE * ATG2 + DD / BB * ATG3 + HH / FF * ATG4);
                U_vec[offset + thread_idx] = 0.5 * krep[thread_idx] / JJ *
                                  (CC / AA * ATG1 + GG / EE * ATG2 + DD / BB * ATG3 + HH / FF * ATG4);
//                printf("%d %f NP\n", thread_idx, U_i[thread_idx]);
                if (U_i[thread_idx] == __longlong_as_floatingpoint(0x7ff0000000000000)
                    || U_i[thread_idx] != U_i[thread_idx] || U_i[thread_idx] < -1.0) {
                    U_i[thread_idx] = -1.0;
                    culpritID[0] = beadSet[n * thread_idx];
                    culpritID[1] = beadSet[n * thread_idx + 1];
                    culpritID[2] = beadSet[n * thread_idx + 2];
                    culpritID[3] = beadSet[n * thread_idx + 3];
                    int j = 0;
                    while (FField[j] != 0) {
                        culpritFF[j] = FField[j];
                        j++;
                    }
                    j = 0;
                    while (interaction[j] != 0) {
                        culpritinteraction[j] = interaction[j];
                        j++;
                    }
                    printf("Coordiantes \n %f %f %f \n %f %f %f \n %f %f %f \n %f %f %f \n", c1[3 * threadIdx.x ],
                           c1[3 * threadIdx
                                   .x + 1], c1[3 * threadIdx.x + 2], c2[3 * threadIdx.x ], c2[3 * threadIdx
                                    .x + 1], c2[3 * threadIdx.x + 2], c3[3 * threadIdx.x ], c3[3 * threadIdx
                                    .x + 1], c3[3 * threadIdx.x + 2], c4[3 * threadIdx.x ], c4[3 * threadIdx
                                    .x + 1], c4[3 * threadIdx.x + 2]);
                    assert(0);
                    __syncthreads();
                }
            }
        }
    }

    else if(z[0] != 0.0f) {

        extern __shared__ floatingpoint s[];

        floatingpoint *c1 = s;
        floatingpoint *c2 = &c1[3 * blockDim.x];
        floatingpoint *c3 = &c2[3 * blockDim.x];
        floatingpoint *c4 = &c3[3 * blockDim.x];

        floatingpoint d, invDSquare;
        floatingpoint a, b, c, e, F, AA, BB, CC, DD, EE, FF, GG, HH, JJ;
        floatingpoint ATG1, ATG2, ATG3, ATG4;
        int nint = params[1];
        int n = params[0];
        int offset = max(params[2] - 1 , 0);
        const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

        if (thread_idx < nint) {
//            if(thread_idx == 0){
//                printf("Offset %d \n", offset);
//            }
            U_i[thread_idx] = 0.0f;
            for (auto i = 0; i < 3; i++) {
                c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i]
                                          + z[0] * f[3 * beadSet[n * thread_idx] + i];
                c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i]
                                          + z[0] * f[3 * beadSet[n * thread_idx + 1] + i];
                c3[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 2] + i]
                                          + z[0] * f[3 * beadSet[n * thread_idx + 2] + i];
                c4[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 3] + i]
                                          + z[0] * f[3 * beadSet[n * thread_idx + 3] + i];
            }

            //check if parallel
            if (areParallel(c1, c2, c3, c4, 3 * threadIdx.x)) {

                d = twoPointDistance(c1, c3, 3 * threadIdx.x);
                invDSquare = 1 / (d * d);
                U_i[thread_idx] = krep[thread_idx] * invDSquare * invDSquare;
                U_vec[offset + thread_idx] = krep[thread_idx] * invDSquare * invDSquare;

                if (U_i[thread_idx] == __longlong_as_floatingpoint(0x7ff0000000000000) //infinity
                    || U_i[thread_idx] != U_i[thread_idx] || U_i[thread_idx] < (floatingpoint)-1.0) {
                    U_i[thread_idx] = -1.0;
                    culpritID[0] = beadSet[n * thread_idx];
                    culpritID[1] = beadSet[n * thread_idx + 1];
                    culpritID[2] = beadSet[n * thread_idx + 2];
                    culpritID[3] = beadSet[n * thread_idx + 3];
                    int j = 0;
                    while (FField[j] != 0) {
                        culpritFF[j] = FField[j];
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
            } else {
                //check if in same plane
                if (areInPlane(c1, c2, c3, c4, 3 * threadIdx.x)) {
                    //slightly move point
                    movePointOutOfPlane(c1, c2, c3, c4, 2, 0.01f, 3 * threadIdx.x);
                }
                // else
                a = scalarProduct(c1, c2, c1, c2, 3 * threadIdx.x);
                b = scalarProduct(c3, c4, c3, c4, 3 * threadIdx.x);
                c = scalarProduct(c3, c1, c3, c1, 3 * threadIdx.x);
                d = scalarProduct(c1, c2, c3, c4, 3 * threadIdx.x);
                e = scalarProduct(c1, c2, c3, c1, 3 * threadIdx.x);
                F = scalarProduct(c3, c4, c3, c1, 3 * threadIdx.x);

                AA = sqrt(a * c - e * e);
                BB = sqrt(b * c - F * F);

                CC = d * e - a * F;
                DD = b * e - d * F;

                EE = sqrt(a * (b + c - 2 * F) - (d - e) * (d - e));
                FF = sqrt(b * (a + c + 2 * e) - (d + F) * (d + F));

                GG = d * d - a * b - CC;
                HH = CC + GG - DD;
                JJ = c * (GG + CC) + e * DD - F * CC;


                ATG1 = atan((a + e) / AA) - atan(e / AA);
                ATG2 = atan((a + e - d) / EE) - atan((e - d) / EE);
                ATG3 = atan((F) / BB) - atan((F - b) / BB);
                ATG4 = atan((d + F) / FF) - atan((d + F - b) / FF);

                U_i[thread_idx] = 0.5 * krep[thread_idx] / JJ *
                                  (CC / AA * ATG1 + GG / EE * ATG2 + DD / BB * ATG3 + HH / FF * ATG4);
                U_vec[offset + thread_idx] = 0.5 * krep[thread_idx] / JJ *
                                  (CC / AA * ATG1 + GG / EE * ATG2 + DD / BB * ATG3 + HH / FF * ATG4);

                /*if (U_i[thread_idx] == __longlong_as_floatingpoint(0x7ff0000000000000)
                    || U_i[thread_idx] != U_i[thread_idx] || U_i[thread_idx] < -1.0) {
                    U_i[thread_idx] = -1.0;
                    culpritID[0] = beadSet[n * thread_idx];
                    culpritID[1] = beadSet[n * thread_idx + 1];
                    culpritID[2] = beadSet[n * thread_idx + 2];
                    culpritID[3] = beadSet[n * thread_idx + 3];
                    int j = 0;
                    while (FField[j] != 0) {
                        culpritFF[j] = FField[j];
                        j++;
                    }
                    j = 0;
                    while (interaction[j] != 0) {
                        culpritinteraction[j] = interaction[j];
                        j++;
                    }
                    printf("Coordiantes \n %f %f %f \n %f %f %f \n %f %f %f \n %f %f %f \n", c1[3 * threadIdx.x ],
                           c1[3 * threadIdx
                                   .x + 1], c1[3 * threadIdx.x + 2], c2[3 * threadIdx.x ], c2[3 * threadIdx
                                    .x + 1], c2[3 * threadIdx.x + 2], c3[3 * threadIdx.x ], c3[3 * threadIdx
                                    .x + 1], c3[3 * threadIdx.x + 2], c4[3 * threadIdx.x ], c4[3 * threadIdx
                                    .x + 1], c4[3 * threadIdx.x + 2]);
                    printf("Forces \n %f %f %f \n %f %f %f \n %f %f %f \n %f %f %f \n", c1[3 * threadIdx.x ],
                           f[3 * beadSet[n * thread_idx]], f[3 * beadSet[n * thread_idx] + 1], f[3 * beadSet[n *
                           thread_idx] + 2], f[3 * beadSet[n * thread_idx +1]], f[3 * beadSet[n * thread_idx +1] + 1],
                           f[3 * beadSet[n * thread_idx +1] + 2], f[3 * beadSet[n * thread_idx +2]], f[3 * beadSet[n *
                           thread_idx +2] + 1], f[3 * beadSet[n * thread_idx +2] + 2], f[3 * beadSet[n * thread_idx +3]],
                           f[3 * beadSet[n * thread_idx+3] + 1], f[3 * beadSet[n * thread_idx+3] + 2]);
                    assert(0);
                    __syncthreads();
                }*/
            }
        }
    }
}


__global__ void CUDAExclVolRepulsionforce(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *krep, int *params){


    int nint = params[1];
    int n = params[0];
    floatingpoint d, invDSquare, U;
    floatingpoint a, b, c, e, F, AA, BB, CC, DD, EE, FF, GG, HH, JJ, invJJ;
    floatingpoint ATG1, ATG2, ATG3, ATG4;
    floatingpoint A1, A2, E1, E2, B1, B2, F1, F2, A11, A12, A13, A14;
    floatingpoint E11, E12, E13, E14, B11, B12, B13, B14, F11, F12, F13, F14;
    floatingpoint f1[3], f2[3], f3[3], f4[3];

    extern __shared__ floatingpoint s[];
    floatingpoint *c1 = s;
    floatingpoint *c2 = &c1[3 * blockDim.x];
    floatingpoint *c3 = &c2[3 * blockDim.x];
    floatingpoint *c4 = &c3[3 * blockDim.x];
//    floatingpoint *f1 = &c4[3 * blockDim.x];
//    floatingpoint *f2 = &f1[3 * blockDim.x];
//    floatingpoint *f3 = &f2[3 * blockDim.x];
//    floatingpoint *f4 = &f3[3 * blockDim.x];

    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

    if(thread_idx<nint) {
        for(auto i=0;i<3;i++){
            c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
            c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
            c3[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 2] + i];
            c4[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 3] + i];
//            f1[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx] + i];
//            f2[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx + 1] + i];
//            f3[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx + 2] + i];
//            f4[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx + 3] + i];
        }

//    }
//    __syncthreads();
//    if(thread_idx<nint) {

        //check if parallel
        if(areParallel(c1, c2, c3, c4, 3 * threadIdx.x))  {

            d = twoPointDistance(c1, c3, 3 * threadIdx.x);
            invDSquare =  1 / (d * d);
            U = krep[thread_idx] * invDSquare * invDSquare;

            floatingpoint f0 = 4 * krep[thread_idx] * invDSquare * invDSquare * invDSquare;

            f1[0] = - f0 * (c3[0] - c1[0]);
            f1[1] = - f0 * (c3[1] - c1[1]);
            f1[2] = - f0 * (c3[2] - c1[2]);

            f2[0] = - f0 * (c4[0] - c2[0]);
            f2[1] = - f0 * (c4[1] - c2[1]);
            f2[2] = - f0 * (c4[2] - c2[2]);

            f3[0] = f0 * (c3[0] - c1[0]);
            f3[1] = f0 * (c3[1] - c1[1]);
            f3[2] = f0 * (c3[2] - c1[2]);

            f4[0] = f0 * (c4[0] - c2[0]);
            f4[1] = f0 * (c4[1] - c2[1]);
            f4[2] = f0 * (c4[2] - c2[2]);

            for(int i=0;i<3;i++) {
                atomicAdd(&f[3 * beadSet[n * thread_idx]  +i], f1[i]);
                atomicAdd(&f[3 * beadSet[n * thread_idx+1]  +i], f2[i]);
                atomicAdd(&f[3 * beadSet[n * thread_idx+2]  +i], f3[i]);
                atomicAdd(&f[3 * beadSet[n * thread_idx+3]  +i], f4[i]);
            }
        }

        else {

            //check if in same plane
            if(areInPlane(c1, c2, c3, c4, 3 * threadIdx.x)) {

                //slightly move point
                movePointOutOfPlane(c1, c2, c3, c4, 2, 0.01f, 3 * threadIdx.x);
            }

            a = scalarProduct(c1, c2, c1, c2, 3 * threadIdx.x);
            b = scalarProduct(c3, c4, c3, c4, 3 * threadIdx.x);
            c = scalarProduct(c3, c1, c3, c1, 3 * threadIdx.x);
            d = scalarProduct(c1, c2, c3, c4, 3 * threadIdx.x);
            e = scalarProduct(c1, c2, c3, c1, 3 * threadIdx.x);
            F = scalarProduct(c3, c4, c3, c1, 3 * threadIdx.x);

            AA = sqrt(a*c - e*e);
            BB = sqrt(b*c - F*F);

            CC = d*e - a*F;
            DD = b*e - d*F;

            EE = sqrt( a*(b + c - 2*F) - (d - e)*(d - e) );
            FF = sqrt( b*(a + c + 2*e) - (d + F)*(d + F) );

            GG = d*d - a*b - CC;
            HH = CC + GG - DD;
            JJ = c*(GG + CC) + e*DD - F*CC;

            invJJ = 1/JJ;

            ATG1 = atan( (a + e)/AA) - atan(e/AA);
            ATG2 = atan((a + e - d)/EE) - atan((e - d)/EE);
            ATG3 = atan((F)/BB) - atan((F - b)/BB);
            ATG4 = atan((d + F)/FF) - atan((d + F - b)/FF);

            U = 0.5 * krep[thread_idx]/ JJ * ( CC/AA*ATG1 + GG/EE*ATG2 + DD/BB*ATG3 + HH/FF*ATG4);


            A1 = AA * AA / (AA * AA + e * e);
            A2 = AA * AA / (AA * AA + (a + e) * (a + e));

            E1 = EE * EE / (EE * EE + (a + e - d) * (a + e - d));
            E2 = EE * EE / (EE * EE + (e - d) * (e - d));

            B1 = BB * BB / (BB * BB + (F - b) * (F - b));;
            B2 = BB * BB / (BB * BB + F * F);

            F1 = FF * FF / (FF * FF + (d + F - b) * (d + F - b));
            F2 = FF * FF / (FF * FF + (d + F) * (d + F));

            A11 = ATG1 / AA;
            A12 = -((ATG1 * CC) / (AA * AA)) + (A1 * CC * e) / (AA * AA * AA) -
                  (A2 * CC * (a + e)) / (AA * AA * AA);
            A13 = -((A1 * CC) / (AA * AA)) + (A2 * CC) / (AA * AA);
            A14 = (A2 * CC) / (AA * AA);

            E11 = ATG2 / EE;
            E12 = (E2 * (-a + d - e) * GG) / (EE * EE * EE) + (E1 * (-d + e) * GG) / (EE * EE * EE) -
                  (ATG2 * GG) / (EE * EE);
            E13 = -((E1 * GG) / (EE * EE)) + (E2 * GG) / (EE * EE);
            E14 = (E2 * GG) / (EE * EE);

            B11 = ATG3 / BB;
            B12 = -((ATG3 * DD) / (BB * BB)) - (B2 * DD * F) / (BB * BB * BB) +
                  (B1 * DD * (-b + F)) / (BB * BB * BB);
            B13 = -((B1 * DD) / (BB * BB)) + (B2 * DD) / (BB * BB);
            B14 = (B1 * DD) / (BB * BB);

            F11 = ATG4 / FF;
            F12 = (F2 * (-d - F) * HH) / (FF * FF * FF) + (F1 * (-b + d + F) * HH) / (FF * FF * FF) -
                  (ATG4 * HH) / (FF * FF);
            F13 = -((F1 * HH) / (FF * FF)) + (F2 * HH) / (FF * FF);
            F14 = (F1 * HH) / (FF * FF);


            f1[0] =  - 0.5*invJJ*( (c2[3 * threadIdx.x] - c1[3 * threadIdx.x] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*F))/(2*EE) - A11*F + E11*F - 2*U*F*F + (F12*b)/(2*FF)) ) + (c4[3 * threadIdx.x] - c3[3 * threadIdx.x] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*F - F11*F - 2*U*a*F - 4*U*e*F + 2*U*(d*e - a*F) - (B12*F)/BB) +  (c1[3 * threadIdx.x] - c3[3 * threadIdx.x] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*F + 2*U*(b*e - d*F) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );

            f1[1] =  - 0.5*invJJ*( (c2[3 * threadIdx.x + 1] - c1[3 * threadIdx.x + 1] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*F))/(2*EE) - A11*F + E11*F - 2*U*F*F + (F12*b)/(2*FF)) ) + (c4[3 * threadIdx.x + 1] - c3[3 * threadIdx.x + 1] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*F - F11*F - 2*U*a*F - 4*U*e*F + 2*U*(d*e - a*F) - (B12*F)/BB) +  (c1[3 * threadIdx.x + 1] - c3[3 * threadIdx.x + 1] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*F + 2*U*(b*e - d*F) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );

            f1[2] =  - 0.5*invJJ*( (c2[3 * threadIdx.x + 2] - c1[3 * threadIdx.x + 2] ) *( A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF - 2*(A14 + E14 - E11*b - F11*b + 2*U*b*c + (A12*c)/(2*AA) + (E12*(b + c - 2*F))/(2*EE) - A11*F + E11*F - 2*U*F*F + (F12*b)/(2*FF)) ) + (c4[3 * threadIdx.x + 2] - c3[3 * threadIdx.x + 2] ) *(B13 + E13 - A11*a + E11*a - B11*d - 2*E11*d - F11*d + 4*U*c*d - A11*e + E11*e + 2*U*d*e - (E12*a)/EE + (E12*(d - e))/EE + B11*F - F11*F - 2*U*a*F - 4*U*e*F + 2*U*(d*e - a*F) - (B12*F)/BB) +  (c1[3 * threadIdx.x + 2] - c3[3 * threadIdx.x + 2] )* (-A13 - E13 - B11*b + F11*b - A11*d + E11*d + 2*U*b*e + (A12*e)/AA - (E12*(d - e))/EE - 2*U*d*F + 2*U*(b*e - d*F) - (F12*b)/FF + 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12*a)/(2*EE) +(B12*b)/(2*BB) + (F12*b)/(2*FF))) );


            f2[0] =  - invJJ*( (c2[3 * threadIdx.x] - c1[3 * threadIdx.x] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*F))/(2*EE)-A11*F+E11*F-2*U*F*F+(F12*b)/(2*FF) ) + 0.5*(c4[3 * threadIdx.x] - c3[3 * threadIdx.x])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F + 4*U*e*F - (F12*(d + F))/FF)  + 0.5*(c1[3 * threadIdx.x] - c3[3 * threadIdx.x] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF) );

            f2[1] = - invJJ*( (c2[3 * threadIdx.x + 1] - c1[3 * threadIdx.x + 1] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*F))/(2*EE)-A11*F+E11*F-2*U*F*F+(F12*b)/(2*FF) ) + 0.5*(c4[3 * threadIdx.x + 1] - c3[3 * threadIdx.x + 1])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F + 4*U*e*F - (F12*(d + F))/FF)  + 0.5*(c1[3 * threadIdx.x + 1] - c3[3 * threadIdx.x + 1] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF) );

            f2[2] = - invJJ*( (c2[3 * threadIdx.x + 2] - c1[3 * threadIdx.x + 2] )*( A14+E14-E11*b-F11*b+2*U*b*c+(A12*c)/(2*AA)+(E12*(b+c-2*F))/(2*EE)-A11*F+E11*F-2*U*F*F+(F12*b)/(2*FF) ) + 0.5*(c4[3 * threadIdx.x + 2] - c3[3 * threadIdx.x + 2])*(-E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F + 4*U*e*F - (F12*(d + F))/FF)  + 0.5*(c1[3 * threadIdx.x + 2] - c3[3 * threadIdx.x + 2] )* (A13 + E13 + B11*b - F11*b + A11*d - E11*d - 2*U*b*e - (A12*e)/AA + (E12*(d - e))/EE + 2*U*d*F - 2*U*(b*e - d*F) + (F12*b)/FF) );

            f3[0] =  - 0.5*invJJ*( (c2[3 * threadIdx.x] - c1[3 * threadIdx.x] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*F - F11*F - 2*U*d*F - 4*U*e*F + 2*U*(b*e - d*F) - (F12*b)/FF + (F12*(d + F))/FF) + (c4[3 * threadIdx.x] - c3[3 * threadIdx.x] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (c1[3 * threadIdx.x] - c3[3 * threadIdx.x] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) );

            f3[1] =  - 0.5*invJJ*( (c2[3 * threadIdx.x + 1] - c1[3 * threadIdx.x + 1] )*(-A13 - F13 - B11*b + F11*b -
                    A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*F - F11*F - 2*U*d*F - 4*U*e*F + 2*U*(b*e - d*F) - (F12*b)/FF + (F12*(d + F))/FF) + (c4[3 * threadIdx.x + 1] - c3[3 * threadIdx.x + 1] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (c1[3 * threadIdx.x + 1] - c3[3 * threadIdx.x + 1] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) ) ;

            f3[2] =  - 0.5*invJJ*( (c2[3 * threadIdx.x + 2] - c1[3 * threadIdx.x + 2] )*(-A13 - F13 - B11*b + F11*b - A11*d - E11*d - 2*F11*d + 4*U*c*d - A11*e + E11*e + 2*U*b*e + (A12*e)/AA + B11*F - F11*F - 2*U*d*F - 4*U*e*F + 2*U*(b*e - d*F) - (F12*b)/FF + (F12*(d + F))/FF) + (c4[3 * threadIdx.x + 2] - c3[3 * threadIdx.x + 2] )*(-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))) + (c1[3 * threadIdx.x + 2] - c3[3 * threadIdx.x + 2] ) * (-B13 - F13 + A11*a - E11*a + B11*d - F11*d - 2*U*d*e + (E12*a)/EE + 2*U*a*F - 2*U*(d*e - a*F) + (B12*F)/BB + (F12*(d + F))/FF - 2*(-2*U*((-a)*b + d*d) + (A12*a)/(2*AA) + (E12* a)/(2*EE) + (B12*b)/(2*BB) + (F12*b)/(2*FF))) );


            f4[0] =  - invJJ*( 0.5*(c2[3 * threadIdx.x] - c1[3 * threadIdx.x] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F +4*U*e*F - (F12*(d + F))/FF ) + (c4[3 * threadIdx.x] - c3[3 * threadIdx.x])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(c1[3 * threadIdx.x] - c3[3 * threadIdx.x] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*F + 2*U*(d*e - a*F) - (B12*F)/BB - (F12*(d + F))/FF) )  ;

            f4[1] =  - invJJ*( 0.5*(c2[3 * threadIdx.x + 1] - c1[3 * threadIdx.x + 1] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F +4*U*e*F - (F12*(d + F))/FF ) + (c4[3 * threadIdx.x + 1] - c3[3 * threadIdx.x + 1])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(c1[3 * threadIdx.x + 1] - c3[3 * threadIdx.x + 1] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*F + 2*U*(d*e - a*F) - (B12*F)/BB - (F12*(d + F))/FF) ) ;

            f4[2] =  - invJJ*( 0.5*(c2[3 * threadIdx.x + 2] - c1[3 * threadIdx.x + 2] )*( -E13 + F13 + 2*E11*d + 2*F11*d - 4*U*c*d + A11*e - E11*e - (E12*(d - e))/EE - B11*F + F11*F +4*U*e*F - (F12*(d + F))/FF ) + (c4[3 * threadIdx.x + 2] - c3[3 * threadIdx.x + 2])*(B14 + F14 - E11*a - F11*a + 2*U*a*c + B11*e - F11*e - 2*U*e*e + (E12*a)/(2*EE) + (B12*c)/(2*BB) + (F12*(a + c + 2*e))/(2*FF))  + 0.5*(c1[3 * threadIdx.x + 2] - c3[3 * threadIdx.x + 2] )* (B13 + F13 - A11*a + E11*a - B11*d + F11*d + 2*U*d*e - (E12*a)/EE - 2*U*a*F + 2*U*(d*e - a*F) - (B12*F)/BB - (F12*(d + F))/FF) ) ;


            for(int i=0;i<3;i++) {
                atomicAdd(&f[3 * beadSet[n * thread_idx]  +i], f1[i]);
                atomicAdd(&f[3 * beadSet[n * thread_idx+1]  +i], f2[i]);
                atomicAdd(&f[3 * beadSet[n * thread_idx+2]  +i], f3[i]);
                atomicAdd(&f[3 * beadSet[n * thread_idx+3]  +i], f4[i]);
            }

        }
//        for(int i=0;i<3;i++) {
//            f1c[3 * thread_idx + i] =  f1[i];
//            f2c[3 * thread_idx + i] = f2[i];
//            f3c[3 * thread_idx + i] =  f3[i];
//            f4c[3 * thread_idx + i] =  f4[i];
//        }

//            f1c[3 * thread_idx]  =  f1[0];
//        f1c[3 * thread_idx + 1]  = f1[1];
//        f1c[3 * thread_idx + 2]  =  f1[2];
//        f2c[3 * thread_idx]  =  d;
//        f2c[3 * thread_idx + 1]  = e;
//        f2c[3 * thread_idx + 2]  =  F;
//        f3c[3 * thread_idx]  =  AA;
//        f3c[3 * thread_idx + 1]  =  BB;
//        f3c[3 * thread_idx + 2]  =  A1;
//        f4c[3 * thread_idx]  =  A2;
//        f4c[3 * thread_idx + 1]  = 0.0;
//        f4c[3 * thread_idx + 2]  =  0.0;

//        f5c[44 * thread_idx]  =  a;
//        f5c[44 * thread_idx+1]  =  b;
//        f5c[44 * thread_idx+2]  =  c;
//        f5c[44 * thread_idx+3]  =  d;
//        f5c[44 * thread_idx+4]  =  e;
//        f5c[44 * thread_idx+5]  =  F;
//        f5c[44 * thread_idx+6]  =  AA;
//        f5c[44 * thread_idx+7]  =  BB;
//        f5c[44 * thread_idx+8]  =  CC;
//        f5c[44 * thread_idx+9]  =  DD;
//        f5c[44 * thread_idx+10]  =  EE;
//        f5c[44 * thread_idx+11]  =  FF;
//        f5c[44 * thread_idx+12]  =  GG;
//        f5c[44 * thread_idx+13]  =  HH;
//        f5c[44 * thread_idx+14]  =  JJ;
//        f5c[44 * thread_idx+15]  =  ATG1;
//        f5c[44 * thread_idx+16]  =  ATG2;
//        f5c[44 * thread_idx+17]  =  ATG3;
//        f5c[44 * thread_idx+18]  =  ATG4;
//        f5c[44 * thread_idx+19]  =  U;
//        f5c[44 * thread_idx+20]  =  A1;
//        f5c[44 * thread_idx+21]  =  A2;
//        f5c[44 * thread_idx+22]  =  E1;
//        f5c[44 * thread_idx+23]  =  E2;
//        f5c[44 * thread_idx+24]  =  B1;
//        f5c[44 * thread_idx+25]  =  B2;
//        f5c[44 * thread_idx+26]  =  F1;
//        f5c[44 * thread_idx+27]  =  F2;
//        f5c[44 * thread_idx+28]  =  A11;
//        f5c[44 * thread_idx+29]  =  A12;
//        f5c[44 * thread_idx+30]  =  A13;
//        f5c[44 * thread_idx+31]  =  A14;
//        f5c[44 * thread_idx+32]  =  E11;
//        f5c[44 * thread_idx+33]  =  E12;
//        f5c[44 * thread_idx+34]  =  E13;
//        f5c[44 * thread_idx+35]  =  E14;
//        f5c[44 * thread_idx+36]  =  B11;
//        f5c[44 * thread_idx+37]  =  B12;
//        f5c[44 * thread_idx+38]  =  B13;
//        f5c[44 * thread_idx+39]  =  B14;
//        f5c[44 * thread_idx+40]  =  F11;
//        f5c[44 * thread_idx+41]  =  F12;
//        f5c[44 * thread_idx+42]  =  F13;
//        f5c[44 * thread_idx+43]  =  F14;


    }

}



#endif
#endif //CUDA_VEC_CYLINDEREXCLVOLREPULSIONCUDA_H

//
// Created by aravind on 1/25/18.
//

#ifndef CUDA_VEC_BRANCHINGDIHEDRALCOSINECUDA_H
#define CUDA_VEC_BRANCHINGDIHEDRALCOSINECUDA_H
#ifdef CUDAACCL
#include "BranchingDihedralCosine.h"

#include "BranchingDihedral.h"

#include "Bead.h"

#include "MathFunctions.h"
#include "SysParams.h"
#include <limits>
#include <assert.h>
#ifdef CUDAACCL
#include <cuda.h>
#include <cuda_runtime.h>
#endif
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

//__global__ void addvectorBDC(floatingpoint *U, int *params, floatingpoint *U_sum, floatingpoint *U_tot){
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

__global__ void BranchingDihedralCosineenergy(floatingpoint *coord, floatingpoint *force, int *beadSet, floatingpoint *kdih,
                                               floatingpoint *pos, int *params,
                                               floatingpoint *U_i, floatingpoint *z,  int *culpritID,
                                              char* culpritFF, char* culpritinteraction, char* FF, char*
                                                interaction) {
    if(z[0] == 0.0) {
        extern __shared__ floatingpoint s[];
        floatingpoint *c1 = s;
        floatingpoint *c2 = &c1[3 * blockDim.x];
        floatingpoint *c3 = &c2[3 * blockDim.x];
        floatingpoint *c4 = &c3[3 * blockDim.x];
        floatingpoint n1n2;
        floatingpoint mp[3], n1[3], n2[3];

        int nint = params[1];
        int n = params[0];
        const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

        if (thread_idx < nint) {
            for (auto i = 0; i < 3; i++) {
                U_i[thread_idx] = 0.0;
                c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
                c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
                c3[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 2] + i];
                c4[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 3] + i];
            }

        }
        __syncthreads();
        if (thread_idx < nint) {
            midPointCoordinate(mp, c1, c2, pos[thread_idx], 3 * threadIdx.x);

            vectorProductmixedID(n1, mp, c2, mp, c3, 0, 3 * threadIdx.x, 0, 3 * threadIdx.x);
            vectorProductmixedID(n2, c3, c4, mp, c3, 3 * threadIdx.x, 3 * threadIdx.x, 0, 3 * threadIdx.x);

            normalizeVector(n1);
            normalizeVector(n2);
            n1n2 = dotProduct(n1, n2);

            U_i[thread_idx] = kdih[thread_idx] * (1 - n1n2);

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

__global__ void BranchingDihedralCosineenergyz(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *kdih,
                                                floatingpoint *pos, int *params,
                                                floatingpoint *U_i, floatingpoint *z,  int *culpritID,
                                               char* culpritFF, char* culpritinteraction, char* FF, char*
                                                interaction) {
    if(z[0] != 0.0) {
        extern __shared__ floatingpoint s[];
        floatingpoint *c1 = s;
        floatingpoint *c2 = &c1[3 * blockDim.x];
        floatingpoint *c3 = &c2[3 * blockDim.x];
        floatingpoint *c4 = &c3[3 * blockDim.x];
        floatingpoint *f1 = &c4[3 * blockDim.x];
        floatingpoint *f2 = &f1[3 * blockDim.x];
        floatingpoint *f3 = &f2[3 * blockDim.x];
        floatingpoint *f4 = &f3[3 * blockDim.x];

        floatingpoint n1n2;
        floatingpoint mp[3], n1[3], n2[3], zero[3];
        zero[0] = 0;
        zero[1] = 0;
        zero[2] = 0;
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
                f1[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx] + i];
                f2[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx + 1] + i];
                f3[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx + 2] + i];
                f4[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx + 3] + i];
            }

        }
        __syncthreads();

        if (thread_idx < nint) {
            midPointCoordinateStretched(mp, c1, f1, c2, f2, pos[thread_idx], z[0], 3 * threadIdx.x);

            vectorProductStretchedmixedID(n1, mp, zero, c2, f2, mp, zero, c3, f3, z[0], 0, 3 * threadIdx.x, 0,
                                          3 * threadIdx.x);
            vectorProductStretchedmixedID(n2, c3, f3, c4, f4, mp, zero, c3, f3, z[0], 3 * threadIdx.x, 3 * threadIdx.x,
                                          0, 3 * threadIdx.x);

            normalizeVector(n1);
            normalizeVector(n2);
            n1n2 = dotProduct(n1, n2);

            U_i[thread_idx] = kdih[thread_idx] * (1 - n1n2);

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


__global__ void BranchingDihedralCosineforces(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                               floatingpoint *kdih, floatingpoint *pos, int *params){

    extern __shared__ floatingpoint s[];
    floatingpoint *c1 = s;
    floatingpoint *c2 = &c1[3 * blockDim.x];
    floatingpoint *c3 = &c2[3 * blockDim.x];
    floatingpoint *c4 = &c3[3 * blockDim.x];
    floatingpoint N1, N2, n1n2, f0, NN1, NN2, X, D, Y, position;
    floatingpoint n2x, n1y, xd, yd, xx, xy, yy, XD, X1, X2, Y1, Y2, D1, D2, YD;
    floatingpoint mp[3], n1[3], n2[3], zero[3], f1[3], f2[3], f3[3], f4[3]; zero[0] = 0; zero[1] = 0; zero[2] = 0;
    int nint = params[1];
    int n = params[0];
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

    if(thread_idx<nint) {
        for(auto i=0;i<3;i++){
            c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
            c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
            c3[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 2] + i];
            c4[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 3] + i];
        }

    }
    __syncthreads();

    if(thread_idx<nint) {

        midPointCoordinate(mp, c1, c2, pos[thread_idx], 3 * threadIdx.x);

        vectorProductmixedID(n1, mp, c2, mp, c3, 0, 3 * threadIdx.x, 0, 3 * threadIdx.x);
        vectorProductmixedID(n2, c3, c4, mp, c3, 3 * threadIdx.x, 3 * threadIdx.x, 0, 3 * threadIdx.x);

        N1 = sqrt(dotProduct(n1, n1));
        N2 = sqrt(dotProduct(n2, n2));
        n1n2 = dotProduct(n1, n2);

        f0 = kdih[thread_idx]/N1/N2;

        NN1 = n1n2/N1/N1;
        NN2 = n1n2/N2/N2;

        X = sqrt(scalarProductmixedID(mp, c2, mp, c2, 0, 3 * threadIdx.x, 0, 3 * threadIdx.x));
        D = sqrt(scalarProductmixedID(mp, c3, mp, c3, 0, 3 * threadIdx.x, 0, 3 * threadIdx.x));
        Y = sqrt(scalarProductmixedID(c3, c4, c3, c4, 3 * threadIdx.x, 3 * threadIdx.x, 3 * threadIdx
            .x, 3 * threadIdx.x));

        n2x = scalarProductmixedID(zero, n2, mp, c2, 0, 0, 0, 3 * threadIdx.x);
        n1y = scalarProductmixedID(zero, n1, c3, c4, 0, 0, 3 * threadIdx.x, 3 * threadIdx.x);
        xd = scalarProductmixedID(mp, c2, mp, c3, 0, 3 * threadIdx.x, 0, 3 * threadIdx.x);
        yd = scalarProductmixedID(c3, c4, mp, c3, 3 * threadIdx.x, 3 * threadIdx.x, 0, 3 * threadIdx.x);

        xx = scalarProductmixedID(mp, c2, mp, c2, 0, 3 * threadIdx.x, 0, 3 * threadIdx.x);
        xy = scalarProductmixedID(mp, c2, c3, c4, 0, 3 * threadIdx.x, 3 * threadIdx.x, 3 * threadIdx.x);
        yy = scalarProductmixedID(c3, c4, c3, c4, 3 * threadIdx.x, 3 * threadIdx.x, 3 * threadIdx.x, 3 * threadIdx.x);

        XD = n2x/D/X/X/X;
        X1 = -NN2*xd/D/X + yd/D/Y + yd/D/D/X/Y;
        X2 = xd*yd/D/D/X/X/X/Y;
        Y1 = -xd/D/X - xd/D/D/X/Y + NN1*yd/D/Y;
        Y2 = xd*yd/D/D/X/Y/Y/Y;
        D1 = NN2*xx/D/X - xy/D/X-xy/D/Y - 2*xy/D/D/X/Y + NN1*yy/D/Y;
        D2 = xd*xy/D/D/X/X/X/Y;
        YD = n1y/D/Y/Y/Y;

        position = pos[thread_idx];

        //force on b1:
        f1[0] = f0*(- (1 - position)*XD*(1-position)*( (c2[3 * threadIdx.x + 1] - c1[3 * threadIdx.x + 1])*(c3[3 * threadIdx.x + 2] - (1-position)*c1[3 * threadIdx.x + 2] - position*c2[3 * threadIdx.x + 2]) - (c2[3 * threadIdx.x + 2] - c1[3 * threadIdx.x + 2])*(c3[3 * threadIdx.x + 1] - (1-position)*c1[3 * threadIdx.x + 1] - position*c2[3 * threadIdx.x + 1]) ) + (1 - position)*(X1 - X2)*(1-position)*(c2[3 * threadIdx.x] - c1[3 * threadIdx.x]) - (1 - position)*Y1*(c4[3 * threadIdx.x] - c3[3 * threadIdx.x]) + (1 - position)*(D1 + D2)*(c3[3 * threadIdx.x] - (1-position)*c1[3 * threadIdx.x] - position*c2[3 * threadIdx.x]));


        f1[1] = f0*(- (1 - position)*XD*(1-position)*( (c2[3 * threadIdx.x + 2] - c1[3 * threadIdx.x + 2])*(c3[3 * threadIdx.x] - (1-position)*c1[3 * threadIdx.x] -
                position*c2[3 * threadIdx.x]) - (c2[3 * threadIdx.x] - c1[3 * threadIdx.x])*(c3[3 * threadIdx.x + 2] - (1-position)*c1[3 * threadIdx.x + 2] - position*c2[3 * threadIdx.x + 2]) ) + (1 - position)*(X1 - X2)*(1-position)*(c2[3 * threadIdx.x + 1] - c1[3 * threadIdx.x + 1]) - (1 - position)*Y1*(c4[3 * threadIdx.x + 1] - c3[3 * threadIdx.x + 1]) + (1 - position)*(D1 + D2)*(c3[3 * threadIdx.x + 1] - (1-position)*c1[3 * threadIdx.x + 1] - position*c2[3 * threadIdx.x + 1]));

        f1[2] = f0*(- (1 - position)*XD*(1-position)*( (c2[3 * threadIdx.x] - c1[3 * threadIdx.x])*(c3[3 * threadIdx.x + 1] - (1-position)*c1[3 * threadIdx.x + 1] - position*c2[3 * threadIdx.x + 1]) - (c2[3 * threadIdx.x + 1] - c1[3 * threadIdx.x + 1])*(c3[3 * threadIdx.x] - (1-position)*c1[3 * threadIdx.x] - position*c2[3 * threadIdx.x]) ) + (1 - position)*(X1 - X2)*(1-position)*(c2[3 * threadIdx.x + 2] - c1[3 * threadIdx.x + 2]) - (1 - position)*Y1*(c4[3 * threadIdx.x + 2] - c3[3 * threadIdx.x + 2]) + (1 - position)*(D1 + D2)*(c3[3 * threadIdx.x + 2] - (1-position)*c1[3 * threadIdx.x + 2] - position*c2[3 * threadIdx.x + 2]));


        //force on b2:
        f2[0] = f0*( (1 - position)*XD*(1-position)*( (c2[3 * threadIdx.x + 1] - c1[3 * threadIdx.x + 1])*(c3[3 * threadIdx.x + 2] - (1-position)*c1[3 * threadIdx.x + 2] - position*c2[3 * threadIdx.x + 2]) - (c2[3 * threadIdx.x + 2] - c1[3 * threadIdx.x + 2])*(c3[3 * threadIdx.x + 1] - (1-position)*c1[3 * threadIdx.x + 1] - position*c2[3 * threadIdx.x + 1]) ) + (X2 + position*(X1 - X2))*(1-position)*(c2[3 * threadIdx.x] - c1[3 * threadIdx.x]) - position*Y1*(c4[3 * threadIdx.x] - c3[3 * threadIdx.x]) + (position*(D1 + D2) - D2)*(c3[3 * threadIdx.x] - (1-position)*c1[3 * threadIdx.x] - position*c2[3 * threadIdx.x]) );

        f2[1] = f0*( (1 - position)*XD*(1-position)*( (c2[3 * threadIdx.x + 2] - c1[3 * threadIdx.x + 2])*(c3[3 * threadIdx.x] - (1-position)*c1[3 * threadIdx.x] -
                position*c2[3 * threadIdx.x]) - (c2[3 * threadIdx.x] - c1[3 * threadIdx.x])*(c3[3 * threadIdx.x + 2] - (1-position)*c1[3 * threadIdx.x + 2] - position*c2[3 * threadIdx.x + 2]) ) + (X2 + position*(X1 - X2))*(1-position)*(c2[3 * threadIdx.x + 1] - c1[3 * threadIdx.x + 1]) - position*Y1*(c4[3 * threadIdx.x + 1] - c3[3 * threadIdx.x + 1]) + (position*(D1 + D2) - D2)*(c3[3 * threadIdx.x + 1] - (1-position)*c1[3 * threadIdx.x + 1] - position*c2[3 * threadIdx.x + 1]) );

        f2[2] = f0*( (1 - position)*XD*(1-position)*( (c2[3 * threadIdx.x] - c1[3 * threadIdx.x])*(c3[3 * threadIdx.x + 1] - (1-position)*c1[3 * threadIdx.x + 1] - position*c2[3 * threadIdx.x + 1]) - (c2[3 * threadIdx.x + 1] - c1[3 * threadIdx.x + 1])*(c3[3 * threadIdx.x] - (1-position)*c1[3 * threadIdx.x] - position*c2[3 * threadIdx.x]) ) + (X2 + position*(X1 - X2))*(1-position)*(c2[3 * threadIdx.x + 2] - c1[3 * threadIdx.x + 2]) - position*Y1*(c4[3 * threadIdx.x + 2] - c3[3 * threadIdx.x + 2]) + (position*(D1 + D2) - D2)*(c3[3 * threadIdx.x + 2] - (1-position)*c1[3 * threadIdx.x + 2] - position*c2[3 * threadIdx.x + 2]) );

        //force on b3:
        f3[0] = f0*(-YD*( (c4[3 * threadIdx.x + 1] - c3[3 * threadIdx.x + 1])*(c3[3 * threadIdx.x + 2] - (1-position)*c1[3 * threadIdx.x + 2] - position*c2[3 * threadIdx.x + 2]) - (c4[3 * threadIdx.x + 2] - c3[3 * threadIdx.x + 2])*(c3[3 * threadIdx.x + 1] - (1-position)*c1[3 * threadIdx.x + 1] - position*c2[3 * threadIdx.x + 1]) ) - X1*(1-position)*(c2[3 * threadIdx.x] - c1[3 * threadIdx.x]) + (Y1 - Y2)*(c4[3 * threadIdx.x] - c3[3 * threadIdx.x]) + (D2 - D1)*(c3[3 * threadIdx.x] - (1-position)*c1[3 * threadIdx.x] - position*c2[3 * threadIdx.x]));

        f3[1] = f0*(-YD*( (c4[3 * threadIdx.x + 2] - c3[3 * threadIdx.x + 2])*(c3[3 * threadIdx.x] - (1-position)*c1[3 * threadIdx.x] - position*c2[3 * threadIdx.x]) - (c4[3 * threadIdx.x] - c3[3 * threadIdx.x])*(c3[3 * threadIdx.x + 2] - (1-position)*c1[3 * threadIdx.x + 2] - position*c2[3 * threadIdx.x + 2]) ) - X1*(1-position)*(c2[3 * threadIdx.x + 1] - c1[3 * threadIdx.x + 1]) + (Y1 - Y2)*(c4[3 * threadIdx.x + 1] - c3[3 * threadIdx.x + 1]) + (D2 - D1)*(c3[3 * threadIdx.x + 1] - (1-position)*c1[3 * threadIdx.x + 1] - position*c2[3 * threadIdx.x + 1]));

        f3[2] = f0*(-YD*( (c4[3 * threadIdx.x] - c3[3 * threadIdx.x])*(c3[3 * threadIdx.x + 1] - (1-position)*c1[3 * threadIdx.x + 1] - position*c2[3 * threadIdx.x + 1]) - (c4[3 * threadIdx.x + 1] - c3[3 * threadIdx.x + 1])*(c3[3 * threadIdx.x] - (1-position)*c1[3 * threadIdx.x] - position*c2[3 * threadIdx.x]) ) - X1*(1-position)*(c2[3 * threadIdx.x + 2] - c1[3 * threadIdx.x + 2]) + (Y1 - Y2)*(c4[3 * threadIdx.x + 2] - c3[3 * threadIdx.x + 2]) + (D2 - D1)*(c3[3 * threadIdx.x + 2] - (1-position)*c1[3 * threadIdx.x + 2] - position*c2[3 * threadIdx.x + 2]));


        //force on b4:
        f4[0] =f0*( YD*( (c4[3 * threadIdx.x + 1] - c3[3 * threadIdx.x + 1])*(c3[3 * threadIdx.x + 2] - (1-position)*c1[3 * threadIdx.x + 2] - position*c2[3 * threadIdx.x + 2]) - (c4[3 * threadIdx.x + 2] - c3[3 * threadIdx.x + 2])*(c3[3 * threadIdx.x + 1] - (1-position)*c1[3 * threadIdx.x + 1] - position*c2[3 * threadIdx.x + 1]) ) + Y2*(c4[3 * threadIdx.x] - c3[3 * threadIdx.x]) - D2*(c3[3 * threadIdx.x] - (1-position)*c1[3 * threadIdx.x] - position*c2[3 * threadIdx.x]) );

        f4[1] =f0*( YD*( (c4[3 * threadIdx.x + 2] - c3[3 * threadIdx.x + 2])*(c3[3 * threadIdx.x] - (1-position)*c1[3 * threadIdx.x] - position*c2[3 * threadIdx.x]) - (c4[3 * threadIdx.x] - c3[3 * threadIdx.x])*(c3[3 * threadIdx.x + 2] - (1-position)*c1[3 * threadIdx.x + 2] - position*c2[3 * threadIdx.x + 2]) ) + Y2*(c4[3 * threadIdx.x + 1] - c3[3 * threadIdx.x + 1]) - D2*(c3[3 * threadIdx.x + 1] - (1-position)*c1[3 * threadIdx.x + 1] - position*c2[3 * threadIdx.x + 1]) );

        f4[2] =f0*( YD*( (c4[3 * threadIdx.x] - c3[3 * threadIdx.x])*(c3[3 * threadIdx.x + 1] - (1-position)*c1[3 * threadIdx.x + 1] - position*c2[3 * threadIdx.x + 1]) - (c4[3 * threadIdx.x + 1] - c3[3 * threadIdx.x + 1])*(c3[3 * threadIdx.x] - (1-position)*c1[3 * threadIdx.x] - position*c2[3 * threadIdx.x]) ) + Y2*(c4[3 * threadIdx.x + 2] - c3[3 * threadIdx.x + 2]) - D2*(c3[3 * threadIdx.x + 2] - (1-position)*c1[3 * threadIdx.x + 2] - position*c2[3 * threadIdx.x + 2]) );
        for (int i = 0; i < 3; i++) {
            atomicAdd(&f[3 * beadSet[n * thread_idx] + i], f1[i]);
            atomicAdd(&f[3 * beadSet[n * thread_idx + 1] + i], f2[i]);
            atomicAdd(&f[3 * beadSet[n * thread_idx + 2] + i], f3[i]);
            atomicAdd(&f[3 * beadSet[n * thread_idx + 3] + i], f4[i]);
        }


    }
}

#endif
#endif //CUDA_VEC_BRANCHINGDIHEDRALCOSINECUDA_H

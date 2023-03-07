//
// Created by aravind on 1/24/18.
//

#ifndef CUDA_VEC_FILAMENTBENDINGHARMONICCUDA_H
#define CUDA_VEC_FILAMENTBENDINGHARMONICCUDA_H
#ifdef CUDAACCL
#include "FilamentBendingHarmonic.h"

#include "FilamentBending.h"

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

//__global__ void addvectorFBH(floatingpoint *U, int *params, floatingpoint *U_sum, floatingpoint *U_tot){
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

__global__ void FilamentBendingHarmonicenergy(floatingpoint *coord, floatingpoint *force, int *beadSet, floatingpoint *kbend,
                                            floatingpoint *eqt, int *params, floatingpoint *U_i, int *culpritID,
                                              char* culpritFF, char* culpritinteraction, char* FF, char*
                                            interaction) {

    extern __shared__ floatingpoint s[];
    floatingpoint *c1 = s;
    floatingpoint *c2 = &c1[3 * blockDim.x];
    floatingpoint *c3 = &c2[3 * blockDim.x];
    floatingpoint L1, L2, L1L2, l1l2;

    int nint = params[1];
    int n = params[0];
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

    if(thread_idx<nint) {
        for(auto i=0;i<3;i++){
            U_i[thread_idx] =0.0;
            c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
            c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
            c3[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 2] + i];
        }
//    }
//    __syncthreads();
//    if(thread_idx<nint) {
        L1 = sqrt(scalarProduct(c1, c2,
                                c1, c2, 3 * threadIdx.x));
        L2 = sqrt(scalarProduct(c2, c3,
                                c2, c3, 3 * threadIdx.x));

        L1L2 = L1*L2;
        l1l2 = scalarProduct(c1, c2,
                             c2, c3, 3 * threadIdx.x);

        U_i[thread_idx] = kbend[thread_idx] * ( 1 - l1l2 / L1L2 );

        if (fabs(U_i[thread_idx]) == __longlong_as_floatingpoint(0x7ff0000000000000) //infinity
            || U_i[thread_idx] != U_i[thread_idx] || U_i[thread_idx] < -1.0) {
            U_i[thread_idx]=-1.0;
            culpritID[0] = thread_idx;
            culpritID[1] = -1;
            int j = 0;
            while(FF[j]!=0){
                culpritFF[j] = FF[j];
                j++;
            }
            j = 0;
            while(interaction[j]!=0){
                culpritinteraction[j] = interaction[j];
                j++;
            }
            assert(0);
            __syncthreads();
        }
    }
}

__global__ void FilamentBendingHarmonicenergyz(floatingpoint *coord, floatingpoint *f, int *beadSet, floatingpoint *kbend,
                                             floatingpoint *eqt, int *params, floatingpoint *U_i, floatingpoint *z, int *culpritID,
                                             char* culpritFF, char* culpritinteraction, char* FF, char*
                                             interaction) {

    extern __shared__ floatingpoint s[];
    floatingpoint *c1 = s;
    floatingpoint *c2 = &c1[3 * blockDim.x];
    floatingpoint *c3 = &c2[3 * blockDim.x];
    floatingpoint *f1 = &c3[3 * blockDim.x];
    floatingpoint *f2 = &f1[3 * blockDim.x];
    floatingpoint *f3 = &f2[3 * blockDim.x];
    floatingpoint L1, L2, L1L2, l1l2;

    int nint = params[1];
    int n = params[0];
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

    if(thread_idx<nint) {
        U_i[thread_idx] = 0.0;
        for(auto i=0;i<3;i++){
            c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
            c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
            c3[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 2] + i];
            f1[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx] + i];
            f2[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx + 1] + i];
            f3[3 * threadIdx.x + i] = f[3 * beadSet[n * thread_idx + 2] + i];
        }

//    }
//    __syncthreads();
//
//    if(thread_idx<nint) {
        L1 = sqrt(scalarProductStretched(c1, f1, c2, f2,
                                         c1, f1, c2, f2, z[0], 3 * threadIdx.x));
        L2 = sqrt(scalarProductStretched(c2, f2, c3, f3,
                                         c2, f2, c3, f3, z[0], 3 * threadIdx.x));

        L1L2 = L1*L2;
        l1l2 = scalarProductStretched(c1, f1, c2, f2,
                                      c2, f2, c3, f3, z[0], 3 * threadIdx.x);

        U_i[thread_idx] = kbend[thread_idx] * ( 1 - l1l2 / L1L2 );
        if (fabs(U_i[thread_idx]) == __longlong_as_floatingpoint(0x7ff0000000000000) //infinity
            || U_i[thread_idx] != U_i[thread_idx] || U_i[thread_idx] < -1.0) {
            U_i[thread_idx]=-1.0;
            culpritID[0] = thread_idx;
            culpritID[1] = -1;
            int j = 0;
            while(FF[j]!=0){
                culpritFF[j] = FF[j];
                j++;
            }
            j = 0;
            while(interaction[j]!=0){
                culpritinteraction[j] = interaction[j];
                j++;
            }
            assert(0);
            __syncthreads();
        }

    }

}


__global__ void FilamentBendingHarmonicforces(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                            floatingpoint *kbend, floatingpoint *eqt, int *params){
    extern __shared__ floatingpoint s[];
    floatingpoint *c1 = s;
    floatingpoint *c2 = &c1[3 * blockDim.x];
    floatingpoint *c3 = &c2[3 * blockDim.x];
    floatingpoint  f1[3], f2[3], f3[3], L1, L2, l1l2, invL1, invL2, A,B,C, k;

    int nint = params[1];
    int n = params[0];
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

    if(thread_idx<nint) {
        for(auto i=0;i<3;i++){
            c1[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx] + i];
            c2[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 1] + i];
            c3[3 * threadIdx.x + i] = coord[3 * beadSet[n * thread_idx + 2] + i];
        }

//    }
//    __syncthreads();
//
//    if(thread_idx<nint) {
        L1 = sqrt(scalarProduct(c1, c2,
                                c1, c2, 3 * threadIdx.x));
        L2 = sqrt(scalarProduct(c2, c3,
                                c2, c3, 3 * threadIdx.x));

        l1l2 = scalarProduct(c1, c2,
                             c2, c3, 3 * threadIdx.x);
        invL1 = 1/L1;
        invL2 = 1/L2;
        A = invL1*invL2;
        B = l1l2*invL1*A*A*L2;
        C = l1l2*invL2*A*A*L1;

        k = kbend[thread_idx];

        //force on i-1, f = k*(-A*l2 + B*l1):

        f1[0] =  k * ((-c3[3 * threadIdx.x] + c2[3 * threadIdx.x])*A +
                       (c2[3 * threadIdx.x] - c1[3 * threadIdx.x])*B );
        f1[1] =  k * ((-c3[3 * threadIdx.x + 1] + c2[3 * threadIdx.x + 1])*A +
                       (c2[3 * threadIdx.x + 1] - c1[3 * threadIdx.x + 1])*B );
        f1[2] =  k * ((-c3[3 * threadIdx.x + 2] + c2[3 * threadIdx.x + 2])*A +
                       (c2[3 * threadIdx.x + 2] - c1[3 * threadIdx.x + 2])*B );

        //force on i, f = k*(A*(l1-l2) - B*l1 + C*l2):
        f2[0] =  k *( (c3[3 * threadIdx.x] - 2*c2[threadIdx.x] + c1[3 * threadIdx.x])*A -
                       (c2[3 * threadIdx.x] - c1[3 * threadIdx.x])*B +
                       (c3[3 * threadIdx.x] - c2[3 * threadIdx.x])*C );

        f2[1] =  k *( (c3[3 * threadIdx.x + 1] - 2*c2[3 * threadIdx.x + 1] + c1[3 * threadIdx.x + 1])*A -
                       (c2[3 * threadIdx.x + 1] - c1[3 * threadIdx.x + 1])*B +
                       (c3[3 * threadIdx.x + 1] - c2[3 * threadIdx.x + 1])*C );

        f2[2] =  k *( (c3[3 * threadIdx.x + 2] - 2*c2[3 * threadIdx.x + 2] + c1[3 * threadIdx.x + 2])*A -
                       (c2[3 * threadIdx.x + 2] - c1[3 * threadIdx.x + 2])*B +
                       (c3[3 * threadIdx.x + 2] - c2[3 * threadIdx.x + 2])*C );

        //force on i-1, f = k*(A*l - B*l2):
        f3[0] =  k *( (c2[3 * threadIdx.x] - c1[3 * threadIdx.x])*A -
                       (c3[3 * threadIdx.x] - c2[3 * threadIdx.x])*C );

        f3[1] =  k *( (c2[3 * threadIdx.x + 1] - c1[3 * threadIdx.x + 1])*A -
                       (c3[3 * threadIdx.x + 1] - c2[3 * threadIdx.x + 1])*C );

        f3[2] =  k *( (c2[3 * threadIdx.x + 2] - c1[3 * threadIdx.x + 2])*A -
                       (c3[3 * threadIdx.x + 2] - c2[3 * threadIdx.x + 2])*C );

        for (int i = 0; i < 3; i++) {
            if (fabs(f1[i]) == __longlong_as_floatingpoint(0x7ff0000000000000) //infinity
                || f1[i] != f1[i]) {
                printf("Fil. Bend. Force became infinite %f %f %f\n",f1[0], f1[1], f1[2]);
                assert(0);
            }
            if (fabs(f2[i]) == __longlong_as_floatingpoint(0x7ff0000000000000) //infinity
                || f2[i] != f2[i]) {
                printf("Fil. Bend. Force became infinite %f %f %f\n",f2[0], f2[1], f2[2]);
                assert(0);
            }
            if (fabs(f3[i]) == __longlong_as_floatingpoint(0x7ff0000000000000) //infinity
                || f3[i] != f3[i]) {
                printf("Fil. Bend. Force became infinite %f %f %f\n",f3[0], f3[1], f3[2]);
                assert(0);
            }
            atomicAdd(&f[3 * beadSet[n * thread_idx] + i], f1[i]);
            atomicAdd(&f[3 * beadSet[n * thread_idx + 1] + i], f2[i]);
            atomicAdd(&f[3 * beadSet[n * thread_idx + 2] + i], f3[i]);
        }
    }
}

#endif
#endif //CUDA_VEC_FILAMENTBENDINGHARMONICCUDA_H

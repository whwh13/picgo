//
// Created by aravind on 3/3/18.
//

#ifndef CUDA_VEC_BINDINGMANAGERCUDA_H
#define CUDA_VEC_BINDINGMANAGERCUDA_H
#ifdef CUDAACCL
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cmath>
#include "math.h"
#include "utility.h"
#include "assert.h"
#include "MathFunctions.h"
using namespace mathfunc;

__global__ void updateAllPossibleBindingsCUDA(floatingpoint *coord, int *beadSet, int *cylID, int *filID, int
*filType, unsigned int *cmpID,  int *NL, int *numNLpairs, int *numpairs, int *params,
                                              floatingpoint *params2, int *possibleBindings,
int *relevant_spboundvec, int *bindingSites){
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;
    int numBindingSites = params[0];
    int filamenttype = params[1];
    int numMonomerperCyl = params[2];
    floatingpoint rmin = params2[0];
    floatingpoint rmax = params2[1];
    floatingpoint cylf[3], cyls[3], cylnf[3], cylns[3], leg1[3], leg2[3] ;
    bool checkstate = true;
//    printf("Num NL pairs %d %f %f %d\n", numNLpairs[0], rmin, rmax, filamenttype);
    if(thread_idx < numNLpairs[0]){
        int cIndex = NL[2 * thread_idx];
        int cnIndex = NL[2 * thread_idx +1];
//        printf("%d %d %d %d %d\n", filamenttype, filType[cIndex], filType[cnIndex], filID[cIndex], filID[cnIndex]);
        if (filamenttype != filType[cIndex] || filamenttype != filType[cnIndex] ||
            filID[cIndex] == filID[cnIndex]|| cylID[cIndex] <= cylID[cnIndex])
            checkstate = false;
//        printf("c %d %d %d %d \n",  numNLpairs[0], checkstate, cIndex, cnIndex);
//        __syncthreads();
        if(checkstate) {
            //get coordiante
            for (auto i = 0; i < 3; i++) {
                cylf[i] = coord[3 * beadSet[2 * cIndex] + i];
                cyls[i] = coord[3 * beadSet[2 * cIndex + 1] + i];
                cylnf[i] = coord[3 * beadSet[2 * cnIndex] + i];
                cylns[i] = coord[3 * beadSet[2 * cnIndex + 1] + i];
            }
            //loop through binding sites
//            int count = 0;
            for (int bs = 0; bs < numBindingSites; bs++) {
                floatingpoint mp1 = floatingpoint (bindingSites[bs])/floatingpoint(numMonomerperCyl);
                midPointCoordinate(leg1, cylf, cyls, mp1, 0);
                for (int bsn = 0; bsn < numBindingSites; bsn++) {
//                    count++;
                    //make sure binding sites are not occupied
                    if (int(relevant_spboundvec[numBindingSites * cIndex + bs] - 1.0) != 0 ||
                            int(relevant_spboundvec[numBindingSites * cnIndex + bsn] - 1.0) != 0)
                        checkstate = false;
                    else
                        checkstate = true;
//                    printf("2 checkstate %d %d %d %d %d %d %d\n",  numNLpairs[0], checkstate, cIndex, cnIndex,
//                           numBindingSites, bs, bsn);
//                    printf("%d %d %d %d %d\n", checkstate, relevant_spboundvec[numBindingSites * cIndex + bs],
//                           relevant_spboundvec[numBindingSites * cnIndex + bsn], int(relevant_spboundvec[numBindingSites
//                           * cIndex + bs] - 1.0) , int(relevant_spboundvec[numBindingSites * cnIndex + bsn] - 1.0));
                    if(checkstate){
                        floatingpoint mp2 = (floatingpoint) (bindingSites[bsn])/floatingpoint(numMonomerperCyl);
                        midPointCoordinate(leg2, cylnf, cylns, mp2, 0);
//                        printf("%f %f %d %d %d\n", mp1, mp2, bindingSites[bs],bindingSites[bsn], numMonomerperCyl );
//                        printf("%d %d %d %d\n",bindingSites[0],bindingSites[1],bindingSites[2], bindingSites[3]);
//                        printf("c3 %d %d %d %d %f %f %f %f %f %f %f %f\n", cIndex, cnIndex, bs, bsn, leg1[0], leg1[1],
//                               leg1[2], leg2[0], leg2[1], leg2[2],
//                               mp1, mp2);
                        floatingpoint dist = twoPointDistance(leg1, leg2, 0);
//                        printf("%f %f %f\n", dist, rmin, rmax);
                        if(dist >= rmin && dist <= rmax){
//                            printf("c %d %d %d %d %f %f %f %f %f %f %f %f\n", cIndex, cnIndex, bs, bsn, leg1[0], leg1[1],
//                                   leg1[2], leg2[0], leg2[1], leg2[2],
//                                   mp1, mp2);
                            int numpair_prev = atomicAdd(&numpairs[0], 1);
                            possibleBindings[5 * numpair_prev] = cmpID[cIndex];
                            possibleBindings[5 * numpair_prev + 1] = cIndex;
                            possibleBindings[5 * numpair_prev + 2] = bindingSites[bs];
                            possibleBindings[5 * numpair_prev + 3] = cnIndex;
                            possibleBindings[5 * numpair_prev + 4] = bindingSites[bsn];
//                            printf("C %f %d %d %d %d\n", dist, cIndex, bindingSites[bs], cnIndex,
//                                   bindingSites[bsn]);
                            //Make sure comopartment is also added
                        }
                    }
                }
            }
//            printf("checkstate %d %d \n", thread_idx, count);
//            __syncthreads(); //TODO remove later
        }//checkstate
    }
}

__global__ void updateAllPossibleBindingsBrancherCUDA (floatingpoint *coord, int *beadSet, int *cylID, int *filID, int *
                                                      filType, unsigned int *cmpID, int
*numpairs, int *params,
                                                      floatingpoint *distances, int *zone, int *possibleBindings,
                                                      int *bindingSites, int
                                                       *relevant_spboundvec, floatingpoint *beListplane){
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;
    int numBindingSites = params[0];
    int filamenttype = params[1];
    int numMonomerperCyl = params[2];
    int numcyl = params[3];
    floatingpoint leg1[3];
    bool checkstate = true;
    if(thread_idx < numcyl){
        int cIndex = thread_idx;
        if (filamenttype != filType[cIndex])
            checkstate = false;
        if(checkstate){
            floatingpoint cylf[3], cyls[3];
            for(int i = 0; i < 3; i++) {
                cylf[i] = coord[3 * beadSet[2 * cIndex] + i];
                cyls[i] = coord[3 * beadSet[2 * cIndex + 1] + i];
            }
            for (int bs = 0; bs < numBindingSites; bs++) {
                floatingpoint dist = -1.0;
                if (zone[0] > 0) {// Type 0 = ALL
                    floatingpoint mp1 = floatingpoint(bindingSites[bs]) / floatingpoint(numMonomerperCyl);
                    midPointCoordinate(leg1, cylf, cyls, mp1, 0);
                //Get distance from the closest boundary.
                    for (auto b = 0; b < 6; b++) {
                        floatingpoint plane[4];
                        for (auto i = 0; i < 4; i++)
                            plane[i] = beListplane[4*b + i];
                        floatingpoint dist_temp = getdistancefromplane(leg1, plane, 0);
                        if(dist == -1.0)
                            dist = dist_temp;
                        else if (dist_temp < dist)
                            dist = dist_temp;
                    }
                    //Check if within nucleationzone

                    if(dist < distances[0]) {
                        if (zone[0] == 2) {
                            if (leg1[2] >= distances[1])
                                checkstate = true;
                            else
                                checkstate = false;
                        }
                        else
                            checkstate = true;
                    }
                    else
                        checkstate = false;
                }

                //make sure binding sites are not occupied
                if (int(relevant_spboundvec[numBindingSites * cIndex + bs] - 1.0) != 0)
                    checkstate = false;
//                printf("%d %d %d %f %f %f\n",thread_idx, bindingSites[bs], zone[0],dist, distances[0],
//                       distances[1]);
//                else
//                    checkstate = true;
                if(checkstate){
//                    printf("%d %d %d %f %f %f\n",thread_idx, bindingSites[bs], zone[0],dist, distances[0],
//                       distances[1]);
                    int numpair_prev = atomicAdd(&numpairs[0], 1);
                    possibleBindings[3 * numpair_prev] = cmpID[cIndex];
                    possibleBindings[3 * numpair_prev + 1] = cIndex;
                    possibleBindings[3 * numpair_prev + 2] = bindingSites[bs];
                }
            }
        }
    }
}
#endif
#endif //CUDA_VEC_BINDINGMANAGERCUDA_H

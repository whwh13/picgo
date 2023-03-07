
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include <cmath>
#include <vector>
#include <math.h>

#include "MathFunctions.h"
#include "Rand.h"

namespace medyan {
namespace mathfunc {

    tuple<vector<floatingpoint>, vector<floatingpoint>> branchProjection(const vector<floatingpoint>& n,
                                                           const vector<floatingpoint>& p,
                                                           floatingpoint l, floatingpoint m, floatingpoint theta){
        //get random permutation from p
        vector<floatingpoint> r = {p[0] + Rand::randfloatingpoint(-1, 1),
                            p[1] + Rand::randfloatingpoint(-1, 1),
                            p[2] + Rand::randfloatingpoint(-1, 1)};

        //construct vector z which is r-p
        auto z = twoPointDirection(p, r);

        //construct u and v, which creates an orthogonal set n, u, v
        auto u = crossProduct(n, z);
        auto v = crossProduct(n, u);

        normalize(u); normalize(v);

        //find random point on circle defining the branching point
        floatingpoint thetaRandom = Rand::randfloatingpoint(0, 2*M_PI);
        vector<floatingpoint> bp1;
        bp1.push_back(p[0] + l * (u[0] * cos(thetaRandom) + v[0] * sin(thetaRandom)));
        bp1.push_back(p[1] + l * (u[1] * cos(thetaRandom) + v[1] * sin(thetaRandom)));
        bp1.push_back(p[2] + l * (u[2] * cos(thetaRandom) + v[2] * sin(thetaRandom)));

        //now find the second point
        vector<floatingpoint> newP;
        floatingpoint dist = m * cos(theta);
        newP.push_back(p[0] + n[0] * dist);
        newP.push_back(p[1] + n[1] * dist);
        newP.push_back(p[2] + n[2] * dist);
        floatingpoint newL = (l + m * sin(theta));

        vector<floatingpoint> bp2;
        bp2.push_back(newP[0] + newL * (u[0] * cos(thetaRandom) + v[0] * sin(thetaRandom)));
        bp2.push_back(newP[1] + newL * (u[1] * cos(thetaRandom) + v[1] * sin(thetaRandom)));
        bp2.push_back(newP[2] + newL * (u[2] * cos(thetaRandom) + v[2] * sin(thetaRandom)));

        //get direction
        auto direction = twoPointDirection(bp1, bp2);

        return tuple<vector<floatingpoint>, vector<floatingpoint>>(direction, bp1);
    }

#ifdef CUDAACCL
//     __global__ void addvector(floatingpoint *U, int *params, floatingpoint *U_sum, floatingpoint *U_tot){
//        U_sum[0] = 0.0;
//        floatingpoint sum = 0.0;
//        for(auto i=0;i<params[1];i++){
//            if(U[i] == -1.0 && sum != -1.0){
//                U_sum[0] = -1.0;
//                U_tot[0] = -1.0;
//                sum = -1.0;
//                break;
//            }
//            else
//                sum  += U[i];
//        }
//        U_sum[0] = sum;
//        atomicAdd(&U_tot[0], sum);
//        printf("ADD1 %f\n", sum);
//
//    }

//    __global__ void addvectorred(floatingpoint *U, int *params, floatingpoint *U_sum, floatingpoint *U_tot){
//        extern __shared__ floatingpoint s[];
//        floatingpoint *c1 = s;
//        int start = 0;
//        int end = params[1];
//        int factor = params[1]/blockDim.x;
//        if(factor > 0) {
//            if (threadIdx.x > 0)
//                start = threadIdx.x * factor;
//            if (threadIdx.x < blockDim.x - 1)
//                end = (threadIdx.x + 1) * factor;
//            c1[threadIdx.x] = 0.0;
//            for (auto i = start; i < end; i++) {
//                if (c1[threadIdx.x] != -1.0 || U[i] != -1.0)
//                    c1[threadIdx.x] += U[i];
//                else
//                    c1[threadIdx.x] = -1.0;
//            }
////    printf("%d \n", params[0]);
//            __syncthreads();
//
//            if (threadIdx.x == 0) {
//                U_sum[0] = 0.0;
//                for (auto i = 0; i < blockDim.x; i++) {
//                    if (c1[i] != -1.0 && U_sum[0] != -1.0)
//                        U_sum[0] += c1[i];
//                    else
//                        U_sum[0] = -1.0;
//                }
//            if(U_sum[0] == -1.0){
//                auto val = -U_tot[0]-1.0;
//                atomicAdd(&U_tot[0], val);
//            }
//            else
//                atomicAdd(&U_tot[0], U_sum[0]);
////                printf("ADD2 %f %f \n", U_tot[0], U_sum[0]);
//            }
//        }
//        else if(threadIdx.x == 0) {
//            U_sum[0] = 0.0;
//            floatingpoint sum = 0.0;
//            for (auto i = 0; i < params[1]; i++) {
//                if (U[i] == -1.0 && sum != -1.0) {
//                    U_sum[0] = -1.0;
//                    U_tot[0] = -1.0;
//                    sum = -1.0;
//                    break;
//                } else
//                    sum += U[i];
//            }
//            U_sum[0] = sum;
//            atomicAdd(&U_tot[0], sum);
////            printf("ADD3 %f %f \n", U_tot[0], U_sum[0]);
//        }
//    }

    __global__ void addvectorred2(floatingpoint *g_idata, int *num, floatingpoint *U_sum, floatingpoint *g_odata)
    {
        extern __shared__ floatingpoint sdata[];
        unsigned int tid = threadIdx.x;
        auto blockSize = blockDim.x;
        unsigned int i = blockIdx.x*(blockSize*2) + tid;
        unsigned int gridSize = blockSize*2*gridDim.x;
        int n = num[1];
        sdata[tid] = 0;
        while (i < n) {
//            printf("%d %d %d\n", tid, i, blockSize);
            if(g_idata[i] == -1.0 || g_idata[i+blockSize] == -1.0)
            {printf("CUDA addition of energies. Energy is infinite\n");assert(0);}
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
            if (blockSize >= 64) {sdata[tid] += sdata[tid + 32]; __syncthreads(); }
            if (blockSize >= 32) {sdata[tid] += sdata[tid + 16]; __syncthreads(); }
            if (blockSize >= 16) {sdata[tid] += sdata[tid + 8]; __syncthreads(); }
            if (blockSize >= 8) {sdata[tid] += sdata[tid + 4]; __syncthreads(); }
            if (blockSize >= 4) {sdata[tid] += sdata[tid + 2]; __syncthreads(); }
            if (blockSize >= 2) {sdata[tid] += sdata[tid + 1]; __syncthreads(); }
        }
        if (tid == 0) {
            atomicAdd(&g_odata[0], sdata[0]);
            atomicAdd(&U_sum[0], sdata[0]);
//        printf("addv2 %f %f %f\n", sdata[0], sdata[1], sdata[2]);
        }
    }

    __global__ void addvectorred3(floatingpoint *g_idata, int *num, floatingpoint *U_sum)
    {
        extern __shared__ floatingpoint sdata[];
        unsigned int tid = threadIdx.x;
        auto blockSize = blockDim.x;
        unsigned int i = blockIdx.x*(blockSize*2) + tid;
        unsigned int gridSize = blockSize*2*gridDim.x;
        int n = num[0];
        sdata[tid] = 0.0;
        while (i < n) {
//            printf("%d %d %d %d\n", tid, i, i+blockSize,num[0]);
            if(g_idata[i] == -1.0 || g_idata[i+blockSize] == -1.0)
            {printf("CUDA addition of energies. Energy is infinite\n");assert(0);}
            sdata[tid] += g_idata[i] + g_idata[i+blockSize];
            i += gridSize;
        }
        /*if(tid ==0) printf("%d %f %d\n", blockIdx.x, sdata[0],n);*/
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
        if (tid == 0) {
            /*printf("%d %f \n", blockIdx.x, sdata[0]);*/
            atomicAdd(&U_sum[0], sdata[0]);
        }
    }

    __global__ void resetintvariableCUDA(int *variable){
        variable[0] = 0;
    }
    __global__ void resetfloatingpointvariableCUDA(floatingpoint *variable){
        variable[0] = 0.0;
    }
#endif
//    __global__ void addvector(floatingpoint *U, int *params, floatingpoint *U_sum, floatingpoint *U_tot, int *culpritID, char* culpritFF,
//                              char* culpritinteraction, char* FF, char* interaction){
//        U_sum[0] = 0.0;
//        floatingpoint sum = 0.0;
//        for(auto i=0;i<params[1];i++){
//            if(U[i] == -1.0 && sum != -1.0){
//                U_sum[0] = -1.0;
//                U_tot[0] = -1.0;
//                sum = -1.0;
//            }
//            else if(sum != -1.0)
//                sum  += U[i];
//        }
//        U_sum[0] = sum;
////        printf("sum %f\n",sum);
//        if(sum != -1.0)
//            atomicAdd(&U_tot[0], sum);
//        else{
//            assert(0);
//        }
//    }
    
    float delGGenChem(float delGZero, vector<species_copy_t> reacN, vector<int> reacNu, vector<species_copy_t> prodN, vector<int> prodNu){
        using std::log;
        if(reacN.size() != reacNu.size()){
            cout << "Reactant copy number vector and stoichiometry vector mismatch" << endl;
            return 0;
        }
        
        if(prodN.size() != prodNu.size()){
            cout << "Product copy number vector and stoichiometry vector mismatch" << endl;
            return 0;
        }
        
        vector<int> reacNminNu;
        
        for(int i=0; (i<reacN.size()); i++){
            reacNminNu.push_back(reacN[i]-reacNu[i]);
        }
        
        vector<int> prodNplusNu;
        
        for(int i=0; (i<prodN.size()); i++){
            prodNplusNu.push_back(prodN[i]+prodNu[i]);
        }
        
        float sumreacs = 0;
        
        for (int i=0; (i<reacN.size()); i++){
            sumreacs += reacNminNu[i]*log(reacNminNu[i]) - reacN[i]*log(reacN[i]) ;
        }
        
        float sumprods = 0;
        
        for (int i=0; (i<prodN.size()); i++){
            sumprods += prodNplusNu[i]*log(prodNplusNu[i]) - prodN[i]*log(prodN[i]) ;
        }
        
        float delG = delGZero +  sumreacs +  sumprods;
        
        return delG;
    }
    
    float delGGenChemI(float delGZero, vector<species_copy_t> reacN, vector<int> reacNu, vector<species_copy_t> prodN, vector<int> prodNu){
        using std::log;
        if(reacN.size() != reacNu.size()){
            cout << "Reactant copy number vector and stoichiometry vector mismatch" << endl;
            return 0;
        }
        
        if(prodN.size() != prodNu.size()){
            cout << "Product copy number vector and stoichiometry vector mismatch" << endl;
            return 0;
        }
        
        vector<int> reacNminNu;
        
        for(int i=0; (i<reacN.size()); i++){
            reacNminNu.push_back(-reacNu[i]);
        }
        
        vector<int> prodNplusNu;
        
        for(int i=0; (i<prodN.size()); i++){
            prodNplusNu.push_back(+prodNu[i]);
        }
        
        float sumreacs = 0;
        
        for (int i=0; (i<reacN.size()); i++){
            sumreacs += reacNminNu[i]*log(reacN[i]) ;
        }
        
        float sumprods = 0;
        
        for (int i=0; (i<prodN.size()); i++){
            sumprods += prodNplusNu[i]*log(prodN[i]) ;
        }
        
        float delG = delGZero +  sumreacs +  sumprods;
        
        return delG;
    }
    
    float delGDifChem(species_copy_t reacN, species_copy_t prodN){
        using std::log;
        if(prodN==0){
            prodN=1;
        }
        
        // float delG =  - log(reacN) + log(prodN);
        float delG;
        delG = (reacN-1)*log(reacN-1) - reacN*log(reacN) + (prodN+1)*log(prodN+1) - prodN*log(prodN);
        //delG = log(prodN) -log(reacN);
        
        return delG;
    }
    
    
    
    float delGPolyChem(float delGzero, species_copy_t reacN, string whichWay){
        
        using std::log;

        float sumreacs=0;
        
        if(whichWay=="P"){
            sumreacs = (reacN-1)*log(reacN-1) - reacN*log(reacN);
        } else if (whichWay=="D"){
            sumreacs = ((reacN+1)*log(reacN+1) - reacN*log(reacN));
        }
        
        float delG =  delGzero +   sumreacs ;
        
        return delG;
        
    }
    
    float delGMyoChem(float nh, float rn){
        float delG;
        delG = - std::log(pow(rn,nh)-1);
        
        return delG;
    }

} // namespace mathfunc

} // namespace medyan

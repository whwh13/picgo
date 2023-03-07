
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

#include "BranchingBendingCosine.h"
#include "BranchingBending.h"
#include "BranchingBendingCosineCUDA.h"

#include "BranchingPoint.h"
#include "Bead.h"

#include "MathFunctions.h"
#include "Cylinder.h"
#ifdef CUDAACCL
#include "nvToolsExt.h"
#endif

namespace medyan {
using namespace mathfunc;
#ifdef CUDAACCL
void BranchingBendingCosine::deallocate(){
    if(!(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamDestroy(stream));
    CUDAcommon::handleerror(cudaFree(gU_i));
    CUDAcommon::handleerror(cudaFree(gU_sum));
    CUDAcommon::handleerror(cudaFree(gFF));
    CUDAcommon::handleerror(cudaFree(ginteraction));
}
void BranchingBendingCosine::optimalblocksnthreads( int nint){
    //CUDA stream create
    if(stream == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamCreate(&stream));
    blocksnthreadse.clear();
    blocksnthreadsez.clear();
    blocksnthreadsf.clear();
    int blockSize;   // The launch configurator returned block size
    int minGridSize; // The minimum grid size needed to achieve the
    // maximum occupancy for a full device launch
    if(nint>0) {
        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       BranchingBendingCosineenergy, blockToSmem, 0);
        blocksnthreadse.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadse.push_back(blockSize);
        blockSize = 0;

        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       BranchingBendingCosineenergyz, blockToSmemez, 0);
        blocksnthreadsez.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadsez.push_back(blockSize);
        blockSize = 0;

        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       BranchingBendingCosineforces, blockToSmem, 0);
        blocksnthreadsf.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadsf.push_back(blockSize);

//        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, nint*sizeof(floatingpoint)));
//        CUDAcommon::handleerror(cudaMalloc((void **) &gU_sum, sizeof(floatingpoint)));
        //get addition vars
        bntaddvec2.clear();
        bntaddvec2 = getaddred2bnt(nint);
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, bntaddvec2.at(0)*sizeof(floatingpoint)));
        CUDAcommon::handleerror(cudaMemset(gU_i, 0, bntaddvec2.at(0) * sizeof(floatingpoint)));
//        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, nint*sizeof(floatingpoint)));
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_sum, sizeof(floatingpoint)));

        char a[] = "BranchingFF";
        char b[] = "Branching Bending Cosine";
        CUDAcommon::handleerror(cudaMalloc((void **) &gFF, 100 * sizeof(char)));
        CUDAcommon::handleerror(cudaMalloc((void **) &ginteraction, 100 * sizeof(char)));
        CUDAcommon::handleerror(cudaMemcpy(gFF, a, 100 * sizeof(char), cudaMemcpyHostToDevice));
        CUDAcommon::handleerror(cudaMemcpy(ginteraction, b, 100 * sizeof(char), cudaMemcpyHostToDevice));
    }
    else{
        blocksnthreadse.push_back(0);
        blocksnthreadse.push_back(0);
        blocksnthreadsez.push_back(0);
        blocksnthreadsez.push_back(0);
        blocksnthreadsf.push_back(0);
        blocksnthreadsf.push_back(0);
    }

}
floatingpoint* BranchingBendingCosine::energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                       floatingpoint *kbend, floatingpoint *eqt, int *params) {
//    if(blocksnthreadse[1]>0) {
//        BranchingBendingCosineenergy<<<blocksnthreadse[0], blocksnthreadse[1], (12 * blocksnthreadse[1]) * sizeof
//                (floatingpoint), stream>>> (coord, f, beadSet, kbend, eqt, params, gU_i, CUDAcommon::getCUDAvars().gculpritID,
//                CUDAcommon::getCUDAvars().gculpritFF,
//                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//        CUDAcommon::handleerror( cudaGetLastError() ,"BranchingBendingCosineenergy", "BranchingBendingCosine.cu");
//        floatingpoint* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror( cudaGetLastError(),"BranchingBendingCosineenergy", "BranchingBendingCosine.cu");
//        return gU_sum;}
//    else
//        return NULL;
}


floatingpoint* BranchingBendingCosine::energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                       floatingpoint *kbend, floatingpoint *eqt, floatingpoint *z, int *params) {
        if(blocksnthreadse[1]>0) {
        BranchingBendingCosineenergy<<<blocksnthreadse[0], blocksnthreadse[1], (12 * blocksnthreadse[1]) * sizeof
                (floatingpoint), stream>>> (coord, f, beadSet, kbend, eqt, params, gU_i, z, CUDAcommon::getCUDAvars()
                .gculpritID,
                CUDAcommon::getCUDAvars().gculpritFF,
                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
            CUDAcommon::handleerror( cudaGetLastError() ,"BranchingBendingCosineenergy", "BranchingBendingCosine.cu");
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//        floatingpoint* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror( cudaGetLastError(),"BranchingBendingCosineenergy", "BranchingBendingCosine.cu");
//        return gU_sum;
        }
    if(blocksnthreadsez[1]>0) {

        BranchingBendingCosineenergyz << < blocksnthreadsez[0], blocksnthreadsez[1], (24 * blocksnthreadsez[1]) *
                                      sizeof(floatingpoint), stream>> > (coord, f, beadSet, kbend, eqt, params, gU_i, z,
                CUDAcommon::getCUDAvars().gculpritID,
                CUDAcommon::getCUDAvars().gculpritFF,
                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction );
        CUDAcommon::handleerror(cudaGetLastError(),"BranchingBendingCosineenergyz", "BranchingBendingCosine.cu");
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//        floatingpoint* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror(cudaGetLastError(),"BranchingBendingCosineenergyz", "BranchingBendingCosine.cu");

//        BranchingBendingCosineenergyz2 << < blocksnthreadsf[0], blocksnthreadsf[1], (12 * blocksnthreadsf[1]) *
//                                                                                     sizeof(floatingpoint), stream>> > (coord, f, beadSet, kbend, eqt, params, gU_i, z );
//        CUDAcommon::handleerror(cudaGetLastError(),"BranchingBendingCosineenergyz2", "BranchingBendingCosine.cu");
//        return gU_sum;
    }
    if(blocksnthreadse[1]<=0 && blocksnthreadsez[1]<=0)
        return NULL;
    else {
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
        floatingpoint *gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;

//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        cudaStreamSynchronize(stream);
//        addvectorred<<<1,200,200*sizeof(floatingpoint),stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        cudaStreamSynchronize(stream);
//        std::cout<<"bntaddvec "<<bntaddvec2.at(0)<<" "<<bntaddvec2.at(1)<<" "<<bntaddvec2.at(0)<<" "
//                ""<<bntaddvec2.at(2)<<" "<<bntaddvec2.at(3)<<endl;
        resetfloatingpointvariableCUDA << < 1, 1, 0, stream >> > (gU_sum);
        addvectorred2 << < bntaddvec2.at(2), bntaddvec2.at(3), bntaddvec2.at(3) * sizeof(floatingpoint), stream >> > (gU_i,
                params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror(cudaDeviceSynchronize(),"FilamentBendingCosineenergyz", "FilamentBendingCosine.cu");
        CUDAcommon::handleerror(cudaGetLastError(), "FilamentBendingCosineenergyz", "FilamentBendingCosine.cu");
        return gU_sum;
    }
}

void BranchingBendingCosine::forces(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                    floatingpoint *kbend, floatingpoint *eqt, int *params){
    if(blocksnthreadsf[1]>0) {
        BranchingBendingCosineforces << < blocksnthreadsf[0], blocksnthreadsf[1], (12 * blocksnthreadsf[1]) *
                                                                                  sizeof(floatingpoint), stream >> > (coord, f, beadSet, kbend, eqt, params);
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
        CUDAcommon::handleerror(cudaGetLastError(),"BranchingBendingCosineforces", "BranchingBendingCosine.cu");
    }
}
void BranchingBendingCosine::checkforculprit() {
    CUDAcommon::printculprit("BranchingBending","BranchingBendingCosine");
    BranchingPoint* br;
    br = (BranchingPoint::getBranchingPoints()[CUDAcommon::getCUDAvars().culpritID[0]]);
    cout<<"Printing culprit branching point information."<<endl;
    br->printSelf();
    exit(EXIT_FAILURE);
}
#endif

floatingpoint BranchingBendingCosine::energy(floatingpoint *coord, int *beadSet,
                                      floatingpoint *kbend, floatingpoint *eqt){

    int n = BranchingBending<BranchingBendingCosine>::n;
    int nint = BranchingPoint::getBranchingPoints().size();

    floatingpoint *coord1, *coord2, *coord3, *coord4, U_i, L1, L2, L1L2, l1l2;

    floatingpoint U = 0.0;

    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[beadSet[n * i]];
        coord2 = &coord[beadSet[n * i + 1]];
        coord3 = &coord[beadSet[n * i + 2]];
        coord4 = &coord[beadSet[n * i + 3]];

        L1 = sqrt(scalarProduct(coord1, coord2,
                                coord1, coord2));
        L2 = sqrt(scalarProduct(coord3, coord4,
                                coord3, coord4));

        L1L2 = L1*L2;
        l1l2 = scalarProduct(coord1, coord2,
                             coord3, coord4);

        floatingpoint x=l1l2 / L1L2;
        if (x < -1.0) x = -1.0;
        else if (x > 1.0) x = 1.0;
        //Option 1 ignore eqt as it is always 0.
        if(areEqual(eqt[i],0.0))
            U_i = kbend[i] * (1 - x);
            //Option 2 Need to calculate Cos(A-B).
        else{
            floatingpoint cosA = x;
            floatingpoint sinA = max<floatingpoint>(sqrt(1-cosA*cosA),(floatingpoint)0.0);
            floatingpoint cosAminusB = cosA*cos(eqt[i]) + sinA*sin(eqt[i]);
            U_i = kbend[i] *(1-cosAminusB);

        }

        U += U_i;
    }

    return U;
}

void BranchingBendingCosine::forces(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                    floatingpoint *kbend, floatingpoint *eqt, floatingpoint *stretchforce){


    int n = BranchingBending<BranchingBendingCosine>::n;
    int nint = BranchingPoint::getBranchingPoints().size();

    floatingpoint *coord1, *coord2, *coord3, *coord4;
	floatingpoint  *force1, *force2, *force3, *force4;
    floatingpoint l1l2, A, B, C, invL1, invL2, k;
    //    floatingpoint L1, L2;
    floatingpoint L1sq, L2sq, x;
    floatingpoint r1x, r1y, r1z, r2x, r2y, r2z;
    floatingpoint Fr1x, Fr1y, Fr1z, Fr2x, Fr2y, Fr2z;

    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[beadSet[n * i]];
        coord2 = &coord[beadSet[n * i + 1]];
        coord3 = &coord[beadSet[n * i + 2]];
        coord4 = &coord[beadSet[n * i + 3]];

        force1 = &f[beadSet[n * i]];
        force2 = &f[beadSet[n * i + 1]];
        force3 = &f[beadSet[n * i + 2]];
        force4 = &f[beadSet[n * i + 3]];

        //Method1
//        L1 = sqrt(scalarProduct(coord1, coord2,
//                                coord1, coord2));
//        L2 = sqrt(scalarProduct(coord3, coord4,
//                                coord3, coord4));
//
////        L1L2 = L1*L2;
//        l1l2 = scalarProduct(coord1, coord2,
//                             coord3, coord4);
//
//        invL1 = 1/L1;
//        invL2 = 1/L2;
//        A = invL1*invL2;
//        B = l1l2*invL1*A*A*L2;
//        C = l1l2*invL2*A*A*L1;

        //Method 2
        L1sq = (scalarProduct(coord1, coord2, coord1, coord2));
        L2sq = (scalarProduct(coord3, coord4, coord3, coord4));

        l1l2 = scalarProduct(coord1, coord2, coord3, coord4);

        A = 1/sqrt(L1sq*L2sq);
        x = l1l2*A;
        B = x/L1sq;
        C = x/L2sq;

        if (areEqual(eqt[i], 0.0)) k = kbend[i];

        else{
            if(abs(abs(x) - 1.0)<0.001)
                x = 0.999*x;

            if (x < -1.0) x = -1.0;
            else if (x > 1.0) x = 1.0;

            floatingpoint cosp =  x;
            floatingpoint sinp = sqrt(max<floatingpoint>((1-cosp*cosp),(floatingpoint)0.0));
            floatingpoint sinpminusq = sinp * cos(eqt[i]) - cosp * sin(eqt[i]);

            k = kbend[i] * sinpminusq/sinp;
        }

        /*	    if(abs(l1l2*A - 1.0)<0.001){
		    l1l2 = 0.999*l1l2;
	    }

        phi = safeacos(l1l2 * A);
        dPhi = phi-eqt[i];

        k =  kbend[i] * sin(dPhi)/sin(phi);*/

        //Bond vectors connecting the 3 beads
        r1x = coord2[0] - coord1[0];
        r1y = coord2[1] - coord1[1];
        r1z = coord2[2] - coord1[2];
        r2x = coord4[0] - coord3[0];
        r2y = coord4[1] - coord3[1];
        r2z = coord4[2] - coord3[2];
        //Forces acting along the bond vectors, r1 and r2.
        Fr1x =  k * ( r2x*A - r1x*B );
        Fr1y =  k * ( r2y*A - r1y*B );
        Fr1z =  k * ( r2z*A - r1z*B );
        Fr2x =  k * ( r1x*A - r2x*C );
        Fr2y =  k * ( r1y*A - r2y*C );
        Fr2z =  k * ( r1z*A - r2z*C );

        //force on i, f = k*(-A*l2 + 2*B*l1):
        force1[0] += -Fr1x;
        force1[1] += -Fr1y;
        force1[2] += -Fr1z;


        //force on i+1, f = k*(A*l2 - 2*B*l1):
        force2[0] += Fr1x;
        force2[1] += Fr1y;
        force2[2] += Fr1z;

        //force on j, k*(-A*l1 + 2*C*l2):
        floatingpoint f3x = -Fr2x;
	    floatingpoint f3y = -Fr2y;
	    floatingpoint f3z = -Fr2z;

	    force3[0] += f3x;
	    force3[1] += f3y;
	    force3[2] += f3z;

	    stretchforce[3*i] = f3x;
	    stretchforce[3*i + 1] = f3y;
	    stretchforce[3*i + 2] = f3z;

        //force on j+1, k*(A*l1 - 2*C*l2):
        force4[0] += Fr2x;
        force4[1] += Fr2y;
        force4[2] += Fr2z;


        #ifdef CHECKFORCES_INF_NAN
        if(checkNaN_INF<floatingpoint>(force1, 0, 2)||checkNaN_INF<floatingpoint>(force2,0,2)||checkNaN_INF<floatingpoint>(force3,0,2)
        ||checkNaN_INF<floatingpoint>(force4,0,2)){
            cout<<"Branching Bending Force becomes infinite. Printing data "<<endl;

            auto b = BranchingPoint::getBranchingPoints()[i];
            auto cyl1 = b->getFirstCylinder();
            auto cyl2 = b->getSecondCylinder();
            cout<<"Cylinder IDs "<<cyl1->getId()<<" "<<cyl2->getId()<<" with cIndex "
                <<cyl1->getStableIndex()<<" "<<cyl2->getStableIndex()<<" and bIndex "
                <<cyl1->getFirstBead()->getStableIndex()<<" "
                <<cyl1->getSecondBead()->getStableIndex()<<" "
                <<cyl2->getFirstBead()->getStableIndex()<<" "
                <<cyl2->getSecondBead()->getStableIndex()<<endl;

            cout<<"Printing intermediate variables"<<endl;
            cout<<"l1l2="<<l1l2<<", A="<<A<<", B="<<B
                <<", C="<<C<<endl;

            cout<<"Printing coords"<<endl;
            cout<<coord1[0]<<" "<<coord1[1]<<" "<<coord1[2]<<endl;
            cout<<coord2[0]<<" "<<coord2[1]<<" "<<coord2[2]<<endl;
            cout<<coord3[0]<<" "<<coord3[1]<<" "<<coord3[2]<<endl;
            cout<<coord4[0]<<" "<<coord4[1]<<" "<<coord4[2]<<endl;
            cout<<"Printing force"<<endl;
            cout<<force1[0]<<" "<<force1[1]<<" "<<force1[2]<<endl;
            cout<<force2[0]<<" "<<force2[1]<<" "<<force2[2]<<endl;
            cout<<force3[0]<<" "<<force3[1]<<" "<<force3[2]<<endl;
            cout<<force4[0]<<" "<<force4[1]<<" "<<force4[2]<<endl;
            cout<<"Printing binary Coords"<<endl;
            printvariablebinary(coord1,0,2);
            printvariablebinary(coord2,0,2);
            printvariablebinary(coord3,0,2);
            printvariablebinary(coord4,0,2);
            cout<<"Printing binary Force"<<endl;
            printvariablebinary(force1,0,2);
            printvariablebinary(force2,0,2);
            printvariablebinary(force3,0,2);
            printvariablebinary(force4,0,2);
            exit(EXIT_FAILURE);
        }
        #endif

    }
}

} // namespace medyan

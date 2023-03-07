
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

#include "FilamentBendingCosine.h"
#include "FilamentBending.h"
#include "FilamentBendingCosineCUDA.h"

#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"
#include "MathFunctions.h"

#ifdef CUDAACCL
#include <cuda.h>
#include <cuda_runtime.h>
#include "nvToolsExt.h"
#endif

namespace medyan {
using namespace mathfunc;
#ifdef CUDAACCL
void FilamentBendingCosine::deallocate(){
//    if(!(CUDAcommon::getCUDAvars().conservestreams))
//        CUDAcommon::handleerror(cudaStreamDestroy(stream));
    CUDAcommon::handleerror(cudaFree(gU_i));
    CUDAcommon::handleerror(cudaFree(gU_sum));
    CUDAcommon::handleerror(cudaFree(gFF));
    CUDAcommon::handleerror(cudaFree(ginteraction));
}
void FilamentBendingCosine::optimalblocksnthreads( int nint, cudaStream_t stream_pass){
//    //CUDA stream create
//    if(stream == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
//        CUDAcommon::handleerror(cudaStreamCreate(&stream));
    stream = stream_pass;
    blocksnthreadse.clear();
    blocksnthreadsez.clear();
    blocksnthreadsf.clear();
    int blockSize;   // The launch configurator returned block size
    int minGridSize; // The minimum grid size needed to achieve the
    // maximum occupancy for a full device launch
    if(nint>0) {
        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       FilamentBendingCosineenergy, blockToSmemFB, 0);
        blocksnthreadse.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadse.push_back(blockSize);
        blockSize = 0;

        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
//                                                       FilamentBendingCosineenergyz, blockToSmemZero, 0);
                                                       FilamentBendingCosineenergyz, blockToSmemFB2, 0);
        blocksnthreadsez.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadsez.push_back(blockSize);
        blockSize = 0;

        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       FilamentBendingCosineforces, blockToSmemFB, 0);
        blocksnthreadsf.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadsf.push_back(blockSize);
        //get addition vars
        bntaddvec2.clear();
        bntaddvec2 = getaddred2bnt(nint);
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, bntaddvec2.at(0)*sizeof(floatingpoint)));
        CUDAcommon::handleerror(cudaMemset(gU_i, 0, bntaddvec2.at(0) * sizeof(floatingpoint)));
//        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, nint*sizeof(floatingpoint)));
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_sum, sizeof(floatingpoint)));

        char a[] = "FilamentFF";
        char b[] = "Filament Bending Cosine";
        CUDAcommon::handleerror(cudaMalloc((void **) &gFF, 100 * sizeof(char)));
        CUDAcommon::handleerror(cudaMalloc((void **) &ginteraction, 100 * sizeof(char)));
        CUDAcommon::handleerror(cudaMemcpyAsync(gFF, a, 100 * sizeof(char),
                                            cudaMemcpyHostToDevice, stream));
        CUDAcommon::handleerror(cudaMemcpyAsync(ginteraction, b, 100 * sizeof(char),
                                            cudaMemcpyHostToDevice, stream));

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
floatingpoint* FilamentBendingCosine::energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                      floatingpoint *kbend, floatingpoint *eqt, int *params) {
//    if(blocksnthreadse[1]>0) {
//        FilamentBendingCosineenergy<<<blocksnthreadse[0], blocksnthreadse[1], (9 * blocksnthreadse[1]) * sizeof
//                (floatingpoint), stream>>> (coord, f, beadSet, kbend, eqt, params, gU_i, CUDAcommon::getCUDAvars()
//                .gculpritID,
//                CUDAcommon::getCUDAvars().gculpritFF,
//                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//        CUDAcommon::handleerror( cudaGetLastError() ,"FilamentBendingCosineenergy", "FilamentBendingCosine.cu");
//        floatingpoint* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror( cudaGetLastError() ,"FilamentBendingCosineenergy", "FilamentBendingCosine.cu");
//        return gU_sum;}
//    else
//        return NULL;
}


floatingpoint* FilamentBendingCosine::energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                      floatingpoint *kbend, floatingpoint *eqt, floatingpoint *z, int *params) {

//    if(blocksnthreadse[1]>0) {
//        FilamentBendingCosineenergy<<<blocksnthreadse[0], blocksnthreadse[1], (9 * blocksnthreadse[1]) * sizeof
//                (floatingpoint), stream>>> (coord, f, beadSet, kbend, eqt, params, gU_i, z, CUDAcommon::getCUDAvars()
//                .gculpritID,
//                CUDAcommon::getCUDAvars().gculpritFF,
//                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
//        CUDAcommon::handleerror(cudaGetLastError(),"FilamentBendingCosineenergy", "FilamentBendingCosine.cu");
//    }

    if(blocksnthreadsez[1]>0) {
        auto boolvarvec = CUDAcommon::cudavars.backtrackbools;
        FilamentBendingCosineenergyz << < blocksnthreadsez[0], blocksnthreadsez[1], (18 * blocksnthreadsez[1]) *
//        FilamentBendingCosineenergyz << < blocksnthreadsez[0], blocksnthreadsez[1], (0) *
                sizeof(floatingpoint), stream>> > (coord, f, beadSet, kbend, eqt, params, gU_i,
                                        CUDAcommon::cudavars.gpu_energyvec, z,
                                        CUDAcommon::getCUDAvars().gculpritID,
                                          CUDAcommon::getCUDAvars().gculpritFF,
                                          CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction, boolvarvec.at(0),
                boolvarvec.at(1) );
        CUDAcommon::handleerror(cudaGetLastError(),"FilamentBendingCosineenergyz", "FilamentBendingCosine.cu");
    }
    if(blocksnthreadse[1]<=0 && blocksnthreadsez[1]<=0)
        return NULL;
    else{
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
#ifdef CUDA_INDIVIDUAL_ESUM
        floatingpoint* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
        resetfloatingpointvariableCUDA<<<1,1,0,stream>>>(gU_sum);
        addvectorred2<<<bntaddvec2.at(2),bntaddvec2.at(3), bntaddvec2.at(3) * sizeof(floatingpoint),stream>>>(gU_i,
                params, gU_sum, gpu_Utot);
#endif
        CUDAcommon::handleerror(cudaGetLastError(),"FilamentBendingCosineenergyz", "FilamentBendingCosine.cu");
        return gU_sum;
    }
}

void FilamentBendingCosine::forces(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                   floatingpoint *kbend, floatingpoint *eqt, int *params){
    if(blocksnthreadsf[1]>0) {
        FilamentBendingCosineforces << < blocksnthreadsf[0], blocksnthreadsf[1], (9 * blocksnthreadsf[1]) *
                                    sizeof(floatingpoint), stream >> > (coord, f, beadSet, kbend, eqt, params);
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
        CUDAcommon::handleerror(cudaGetLastError(),"FilamentBendingCosineforces", "FilamentBendingCosine.cu");
//        CUDAcommon::handleerror(cudaDeviceSynchronize());
//        floatingpoint U_i[Bead::getBeads().size() - 2 * Filament::getFilaments().size()];
//        CUDAcommon::handleerror(cudaMemcpy(U_i, eqt, (Bead::getBeads().size() - 2 * Filament::getFilaments().size()) *
//                                                                               sizeof(floatingpoint), cudaMemcpyDeviceToHost));
//        for(auto i = 0; i< Bead::getBeads().size() - 2 * Filament::getFilaments().size();i ++){
//            std::cout<<U_i[i]<<endl;
//        }
//        std::cout<<endl;
    }
}

void FilamentBendingCosine::checkforculprit() {
    CUDAcommon::printculprit("FilamentBending","FilamentBendingCosine");
    Filament* fil;
    int i = 0;
    bool found = false;
    for (auto f: Filament::getFilaments()) {

        if (f->getCylinderVector().size() > 1){
            i = i + 2 * f->getCylinderVector().size() - 2;
            if(i > CUDAcommon::getCUDAvars().culpritID[0] && !found){
                found = true;
                fil = (Filament*)(Cylinder::getCylinders()[i]->getParent());
            }
        }
    }
    cout<<"Printing culprit Filament information."<<endl;
    fil->printSelf();
    exit(EXIT_FAILURE);
}
#endif
floatingpoint FilamentBendingCosine::energy(floatingpoint *coord, size_t nint, int *beadSet,
                                     floatingpoint *kbend, floatingpoint *eqt) {

    int n = FilamentBending<FilamentBendingCosine>::n;

    floatingpoint *coord1, *coord2, *coord3, L1, L2, L1L2, l1l2;

    floatingpoint U = 0.0, U_i;

    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[beadSet[n * i]];
        coord2 = &coord[beadSet[n * i + 1]];
        coord3 = &coord[beadSet[n * i + 2]];

        L1 = sqrt(scalarProduct(coord1, coord2,
                                coord1, coord2));
        L2 = sqrt(scalarProduct(coord2, coord3,
                                coord2, coord3));

        L1L2 = L1*L2;
        l1l2 = scalarProduct(coord1, coord2,
                             coord2, coord3);

        floatingpoint x = l1l2/L1L2;

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


void FilamentBendingCosine::forces(floatingpoint *coord, floatingpoint *f, size_t nint, int *beadSet,
                                   floatingpoint *kbend, floatingpoint *eqt){

    int n = FilamentBending<FilamentBendingCosine>::n;

    floatingpoint *coord1, *coord2, *coord3, l1l2, invL1, invL2, A,B,C, k;
    floatingpoint *force1, *force2, *force3;
//    floatingpoint L1, L2;
    floatingpoint L1sq, L2sq,x;
    floatingpoint r1x, r1y, r1z, r2x, r2y, r2z;
    floatingpoint Fr1x, Fr1y, Fr1z, Fr2x, Fr2y, Fr2z;

    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[beadSet[n * i]];
        coord2 = &coord[beadSet[n * i + 1]];
        coord3 = &coord[beadSet[n * i + 2]];

        force1 = &f[beadSet[n * i]];
        force2 = &f[beadSet[n * i + 1]];
        force3 = &f[beadSet[n * i + 2]];
		//Method 1
//        L1 = sqrt(scalarProduct(coord1, coord2, coord1, coord2));//|x1|
//        L2 = sqrt(scalarProduct(coord2, coord3, coord2, coord3));
//
//        l1l2 = scalarProduct(coord1, coord2, coord2, coord3);
//
//        invL1 = 1/L1;
//        invL2 = 1/L2;
//        A = invL1*invL2;//1/|x1||x2|
//        B = l1l2*invL1*A*A*L2;//vec(x1).vec(x2)*(1/|x1|)*(1/|x1||x2|)^2*(|x2|)
//        // = vec(x1).vec(x2)/|x1|^3.|x2|
//        C = l1l2*invL2*A*A*L1;//vec(x1).vec(x2)*(1/|x2|)*(1/|x1||x2|)^2(|x1|)
//		// = vec(x1).vec(x2)/|x2|^3.|x1|

		//Method 2
		L1sq = (scalarProduct(coord1, coord2, coord1, coord2));
		L2sq = (scalarProduct(coord2, coord3, coord2, coord3));

	    l1l2 = scalarProduct(coord1, coord2, coord2, coord3);

	    A = 1/sqrt(L1sq*L2sq);
	    x = l1l2*A;
	    B = x/L1sq;
	    C = x/L2sq;

        if (areEqual(eqt[i], 0.0)) k = kbend[i];

        else{
            if(abs(abs(x) - 1.0)<0.001) {
                x = 0.999*x;
            }
	        if (x < -1.0) x = -1.0;
	        else if (x > 1.0) x = 1.0;

            floatingpoint cosp =  x;
            floatingpoint sinp = sqrt(max<floatingpoint>((1-cosp*cosp),(floatingpoint)0.0));
            floatingpoint sinpminusq = sinp * cos(eqt[i]) - cosp * sin(eqt[i]);

/*            phi = safeacos(l1l2 *A);
            dPhi = phi-eqt[i];
			k = kbend[i] * sin(dPhi)/sin(phi);*/
	        k = kbend[i] * sinpminusq/sinp;
        }

        r1x = coord2[0] - coord1[0];
        r1y = coord2[1] - coord1[1];
        r1z = coord2[2] - coord1[2];
        r2x = coord3[0] - coord2[0];
        r2y = coord3[1] - coord2[1];
        r2z = coord3[2] - coord2[2];

        //Force acting along vectors r1 and r2
        Fr1x =  k * ( r2x*A - r1x*B );
	    Fr1y =  k * ( r2y*A - r1y*B );
	    Fr1z =  k * ( r2z*A - r1z*B );
	    Fr2x =  k * ( r1x*A - r2x*C );
	    Fr2y =  k * ( r1y*A - r2y*C );
	    Fr2z =  k * ( r1z*A - r2z*C );

        //force on i-1, f = k*(-A*l2 + B*l1):
        force1[0] += -Fr1x;

        force1[1] += -Fr1y;

        force1[2] += -Fr1z;

        //force on i, f = k*(A*(l1-l2) - B*l1 + C*l2):
        force2[0] +=  Fr1x - Fr2x;

        force2[1] +=  Fr1y - Fr2y;

        force2[2] +=  Fr1z - Fr2z;

        //force on i-1, f = k*(A*l1 - B*l2):
        force3[0] +=  Fr2x;

        force3[1] +=  Fr2y;

        force3[2] +=  Fr2z;

	    #ifdef CHECKFORCES_INF_NAN
	    if(checkNaN_INF<floatingpoint>(force1, 0, 2)||checkNaN_INF<floatingpoint>(force2,0,2)
	       ||checkNaN_INF<floatingpoint>(force3,0,2)){
		    cout<<"Filament Bending Force becomes infinite. Printing data "<<endl;

		    short found = 0;
		    Cylinder *cyl1, *cyl2;
		    for(auto cyl:Cylinder::getCylinders()){
			    auto dbIndex1 = cyl->getFirstBead()->getIndex() * 3;
			    auto dbIndex2 = cyl->getSecondBead()->getIndex() * 3;
			    if(dbIndex1 == beadSet[n * i] && dbIndex2 == beadSet[n * i + 1]) { // FIXME this is unsafe
				    cyl1 = cyl;
				    found++;
				    if(found>=2)
					    break;
			    }
			    else if(dbIndex1 == beadSet[n * i + 1] && dbIndex2 == beadSet[n * i + 2]){ // FIXME this is unsafe
				    cyl2 = cyl;
				    found++;
				    if(found>=2)
					    break;
			    }
		    }
		    cout<<"Cylinder IDs "<<cyl1->getId()<<" "<<cyl2->getId()<<" with cIndex "
		        <<cyl1->getStableIndex()<<" "<<cyl2->getStableIndex()<<" and bIndex "
		        <<cyl1->getFirstBead()->getStableIndex()<<" "
		        <<cyl1->getSecondBead()->getStableIndex()<<" "
		        <<cyl2->getFirstBead()->getStableIndex()<<" "
		        <<cyl2->getSecondBead()->getStableIndex()<<endl;

		    cout<<"Printing coords"<<endl;
		    cout<<coord1[0]<<" "<<coord1[1]<<" "<<coord1[2]<<endl;
		    cout<<coord2[0]<<" "<<coord2[1]<<" "<<coord2[2]<<endl;
		    cout<<coord3[0]<<" "<<coord3[1]<<" "<<coord3[2]<<endl;
		    cout<<"Printing force"<<endl;
		    cout<<force1[0]<<" "<<force1[1]<<" "<<force1[2]<<endl;
		    cout<<force2[0]<<" "<<force2[1]<<" "<<force2[2]<<endl;
		    cout<<force3[0]<<" "<<force3[1]<<" "<<force3[2]<<endl;
		    cout<<"Printing binary Coords"<<endl;
		    printvariablebinary(coord1,0,2);
		    printvariablebinary(coord2,0,2);
		    printvariablebinary(coord3,0,2);
		    cout<<"Printing binary Force"<<endl;
		    printvariablebinary(force1,0,2);
		    printvariablebinary(force2,0,2);
		    printvariablebinary(force3,0,2);
		    exit(EXIT_FAILURE);
	    }
	    #endif

    }
}

} // namespace medyan

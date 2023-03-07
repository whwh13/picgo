
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
#include <math.h>

#include "BranchingPositionCosine.h"
#include "BranchingPositionCosineCUDA.h"
#include "BranchingPosition.h"

#include "BranchingPoint.h"
#include "MathFunctions.h"
#include "Cylinder.h"

namespace medyan {
using namespace mathfunc;
#ifdef CUDAACCL
void BranchingPositionCosine::deallocate(){
    if(!(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamDestroy(stream));
    CUDAcommon::handleerror(cudaFree(gU_i));
    CUDAcommon::handleerror(cudaFree(gU_sum));
    CUDAcommon::handleerror(cudaFree(gFF));
    CUDAcommon::handleerror(cudaFree(ginteraction));
}
void BranchingPositionCosine::optimalblocksnthreads( int nint){
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
                                                       BranchingPositionCosineenergy, blockToSmemFB, 0);
        blocksnthreadse.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadse.push_back(blockSize);
        blockSize = 0;

        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       BranchingPositionCosineenergyz, blockToSmemFB2, 0);
        blocksnthreadsez.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadsez.push_back(blockSize);
        blockSize = 0;

        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       BranchingPositionCosineforces, blockToSmemFB, 0);
        blocksnthreadsf.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadsf.push_back(blockSize);
        //get addition vars
        bntaddvec2.clear();
        bntaddvec2 = getaddred2bnt(nint);
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, bntaddvec2.at(0)*sizeof(floatingpoint)));
        CUDAcommon::handleerror(cudaMemset(gU_i, 0, bntaddvec2.at(0) * sizeof(floatingpoint)));
//        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, nint*sizeof(floatingpoint)));
        CUDAcommon::handleerror(cudaMalloc((void **) &gU_sum, sizeof(floatingpoint)));

//        CUDAcommon::handleerror(cudaMalloc((void **) &gU_i, nint*sizeof(floatingpoint)));
//        CUDAcommon::handleerror(cudaMalloc((void **) &gU_sum, sizeof(floatingpoint)));

        char a[] = "BranchingFF";
        char b[] = "Branching Position Cosine";
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
floatingpoint* BranchingPositionCosine::energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                        floatingpoint *kpos, floatingpoint *pos, int *params) {
//    if(blocksnthreadse[1]>0) {
//        BranchingPositionCosineenergy<<<blocksnthreadse[0], blocksnthreadse[1], (9 * blocksnthreadse[1]) * sizeof
//                (floatingpoint), stream>>> (coord, f, beadSet, kpos, pos, params, gU_i, CUDAcommon::getCUDAvars().gculpritID,
//                CUDAcommon::getCUDAvars().gculpritFF,
//                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//        CUDAcommon::handleerror( cudaGetLastError() ,"BranchingPositionCosineenergy", "BranchingPositionCosine.cu");
//        floatingpoint* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror( cudaGetLastError() ,"BranchingPositionCosineenergy", "BranchingPositionCosine.cu");
//        return gU_sum;}
//    else
//        return NULL;
}


floatingpoint* BranchingPositionCosine::energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                        floatingpoint *kpos, floatingpoint *pos, floatingpoint *z, int *params) {
    if(blocksnthreadse[1]>0) {
        BranchingPositionCosineenergy<<<blocksnthreadse[0], blocksnthreadse[1], (9 * blocksnthreadse[1]) * sizeof
                (floatingpoint), stream>>> (coord, f, beadSet, kpos, pos, params, gU_i, z, CUDAcommon::getCUDAvars()
                .gculpritID,
                CUDAcommon::getCUDAvars().gculpritFF,
                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//        CUDAcommon::handleerror( cudaGetLastError() ,"BranchingPositionCosineenergy", "BranchingPositionCosine.cu");
//        floatingpoint* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror( cudaGetLastError() ,"BranchingPositionCosineenergy", "BranchingPositionCosine.cu");
//        return gU_sum;
    }

    if(blocksnthreadsez[1]>0) {
        BranchingPositionCosineenergyz << < blocksnthreadsez[0], blocksnthreadsez[1], (18 * blocksnthreadsez[1]) *
                                          sizeof(floatingpoint), stream>> > (coord, f, beadSet, kpos, pos, params, gU_i, z,
                CUDAcommon::getCUDAvars().gculpritID,
                CUDAcommon::getCUDAvars().gculpritFF,
                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction );
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//        CUDAcommon::handleerror(cudaGetLastError(),"BranchingPositionCosineenergyz", "BranchingPositionCosine.cu");
//        floatingpoint* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//
//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror(cudaGetLastError(),"BranchingPositionCosineenergyz", "BranchingPositionCosine.cu");
//        return gU_sum;
    }
    if(blocksnthreadse[1]<=0 && blocksnthreadsez[1]<=0)
        return NULL;
    else{
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
        floatingpoint* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;

//        addvector<<<1,1,0,stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        cudaStreamSynchronize(stream);
//        addvectorred<<<1,200,200*sizeof(floatingpoint),stream>>>(gU_i,params, gU_sum, gpu_Utot);
//        cudaStreamSynchronize(stream);
//        std::cout<<"bntaddvec "<<bntaddvec2.at(0)<<" "<<bntaddvec2.at(1)<<" "<<bntaddvec2.at(0)<<" "
//                ""<<bntaddvec2.at(2)<<" "<<bntaddvec2.at(3)<<endl;
        resetfloatingpointvariableCUDA<<<1,1,0,stream>>>(gU_sum);
        addvectorred2<<<bntaddvec2.at(2),bntaddvec2.at(3), bntaddvec2.at(3) * sizeof(floatingpoint),stream>>>(gU_i,
                params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror(cudaDeviceSynchronize(),"FilamentBendingCosineenergyz", "FilamentBendingCosine.cu");
        CUDAcommon::handleerror(cudaGetLastError(),"FilamentBendingCosineenergyz", "FilamentBendingCosine.cu");
        return gU_sum;
    }

}

void BranchingPositionCosine::forces(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                     floatingpoint *kpos, floatingpoint *pos, int *params){
    if(blocksnthreadsf[1]>0) {
        BranchingPositionCosineforces << < blocksnthreadsf[0], blocksnthreadsf[1], (9 * blocksnthreadsf[1]) *
                                                                                   sizeof(floatingpoint), stream >> > (coord, f, beadSet, kpos, pos, params);
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
        CUDAcommon::handleerror(cudaGetLastError(),"BranchingPositionCosineforces", "BranchingPositionCosine.cu");
    }
}
void BranchingPositionCosine::checkforculprit() {
    CUDAcommon::printculprit("BranchingPosition","BranchingPositionCosine");
    BranchingPoint* br;
    br = (BranchingPoint::getBranchingPoints()[CUDAcommon::getCUDAvars().culpritID[0]]);
    cout<<"Printing culprit branching point information."<<endl;
    br->printSelf();
    exit(EXIT_FAILURE);
}
#endif

floatingpoint BranchingPositionCosine::energy(const floatingpoint *coord,
                                                unsigned int *beadSet,
                                                floatingpoint *kpos, floatingpoint *pos) const {


    int n = BranchingPosition<BranchingPositionCosine>::n;
    int nint = BranchingPoint::getBranchingPoints().size();

    const floatingpoint *coord1, *coord2, *coord3;
    floatingpoint X, D, XD, xd, posheta,
    position, U_i;
    floatingpoint *mp = new floatingpoint[3];
    floatingpoint *coord2prime = new floatingpoint[3];
    floatingpoint U = 0.0;

    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[beadSet[n * i]];
        coord2 = &coord[beadSet[n * i + 1]];
        //If the branching point is at the plus end of the cylinder, we end up with
        // singularities in energy and force expressions. To avoid it, we will create a
        // virtual plus end and use that to define vectors.
        position = pos[i];
        if(areEqual(position, 1.0)){
            position = 0.5;//assign a dummy position and extend plus end to new coordinates.
            for(int dim = 0; dim < 3; dim++) {
                coord2prime[dim] = (1 / position) * (coord2[dim] - (1 - position) *
                                    coord1[dim]);//extended plus end coordinate;
            }
            coord2 = &coord2prime[0];
        }
        coord3 = &coord[beadSet[n * i + 2]];
        midPointCoordinate(mp, coord1, coord2, position);
        X = sqrt(scalarProduct(mp, coord2, mp, coord2));
        D = sqrt(scalarProduct(mp, coord3, mp, coord3));

        XD = X * D;

        xd = scalarProduct(mp, coord2, mp, coord3);

        floatingpoint x = xd/XD;

        if(abs(abs(x) - 1.0)<0.001) {
            xd = 0.999 * XD;
            x = xd / XD;
        }

        if (x < -1.0) x = -1.0;
        else if (x > 1.0) x = 1.0;

        floatingpoint cosp =  x;
        posheta = 0.5*M_PI;
        floatingpoint sinp = sqrt(max<floatingpoint>((1-cosp*cosp),(floatingpoint)0.0));
        floatingpoint cospminusq = cosp * cos(posheta) + sinp * sin(posheta);
        U_i = kpos[i] * ( 1 - cospminusq );

        U += U_i;
    }
    delete[] mp;
    delete[] coord2prime;
    return U;
}

floatingpoint BranchingPositionCosine::energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                       floatingpoint *kpos, floatingpoint *pos, floatingpoint d){

    int n = BranchingPosition<BranchingPositionCosine>::n;
    int nint = BranchingPoint::getBranchingPoints().size();

    floatingpoint *coord1, *coord2, *coord3, X, D, XD, xd, posheta, U_i;
    floatingpoint *f1, *f2, *f3;
    floatingpoint *mp = new floatingpoint[3];
    floatingpoint *vzero = new floatingpoint[3]; vzero[0] = 0.0; vzero[1] = 0.0; vzero[2] = 0.0;

    floatingpoint U = 0.0;

    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[beadSet[n * i]];
        coord2 = &coord[beadSet[n * i + 1]];
        coord3 = &coord[beadSet[n * i + 2]];
        f1 = &f[beadSet[n * i]];
        f2 = &f[beadSet[n * i + 1]];
        f3 = &f[beadSet[n * i + 2]];

        midPointCoordinateStretched(mp, coord1, f1, coord2, f2, pos[i], d);
        X = sqrt(scalarProductStretched(mp, vzero, coord2, f2, mp, vzero, coord2, f2, d));
        D = sqrt(scalarProductStretched(mp, vzero, coord3, f3, mp, vzero, coord3, f3, d));

        XD = X * D;

        xd = scalarProductStretched(mp, vzero, coord2, f2, mp, vzero, coord3, f3, d);

        floatingpoint x = xd/XD;

        if(abs(abs(x) - 1.0)<0.001) {
            xd = 0.999 * XD;
            x = xd / XD;
        }

        if (x < -1.0) x = -1.0;
        else if (x > 1.0) x = 1.0;

        floatingpoint cosp =  x;
        posheta = 0.5*M_PI;
        floatingpoint sinp = sqrt(max<floatingpoint>((1-cosp*cosp),(floatingpoint)0.0));
        floatingpoint cospminusq = cosp * cos(posheta) + sinp * sin(posheta);
        U_i = kpos[i] * ( 1 - cospminusq );

        U += U_i;
    }
    delete[] mp;
    delete[] vzero;
    return U;
}

void BranchingPositionCosine::forces(const floatingpoint *coord, floatingpoint *f, unsigned
                                        int *beadSet, floatingpoint *kpos, floatingpoint *pos,
                                     floatingpoint *stretchforce) const {

    int n = BranchingPosition<BranchingPositionCosine>::n;
    int nint = BranchingPoint::getBranchingPoints().size();

    const floatingpoint *coord1, *coord2, *coord3;
    floatingpoint Xsq, Dsq, x, xd, position, A, B, C, k, posheta;
	floatingpoint  *f1, *f2, *f3;
    floatingpoint *mp = new floatingpoint[3];
    floatingpoint *coord2prime = new floatingpoint[3];
    floatingpoint f1tempx, f1tempy, f1tempz, f2tempx, f2tempy, f2tempz, f3tempx, f3tempy,
    f3tempz;
    floatingpoint r1x, r1y, r1z, r2x, r2y, r2z;
    floatingpoint Fr1x, Fr1y, Fr1z, Fr2x, Fr2y, Fr2z;
    floatingpoint fmpx, fmpy, fmpz, fc2x, fc2y, fc2z, fc3x, fc3y, fc3z;

    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[beadSet[n * i]];
        coord2 = &coord[beadSet[n * i + 1]];
        //If the branching point is at the plus end of the cylinder, we end up with
        // singularities in energy and force expressions. To avoid it, we will create a
        // virtual plus end and use that to define vectors.
        position = pos[i];
        if(areEqual(position, (floatingpoint)1.0)){
            position = 0.5;//assign a dummy position and extend plus end to new coordinates.
            for(int dim = 0; dim < 3; dim++) {
                coord2prime[dim] = (1 / position) * (coord2[dim] - (1 - position) *
                                                                   coord1[dim]);//extended plus end coordinate;
            }
            coord2 = &coord2prime[0];
        }
        coord3 = &coord[beadSet[n * i + 2]];
        f1 = &f[beadSet[n * i]];
        f2 = &f[beadSet[n * i + 1]];
        f3 = &f[beadSet[n * i + 2]];

        midPointCoordinate(mp, coord1, coord2, position);
/*        Method I
        X = sqrt(scalarProduct(mp, coord2, mp, coord2));
        D = sqrt(scalarProduct(mp, coord3, mp, coord3));

        XD = X * D;
        xd = scalarProduct(mp, coord2, mp, coord3);
        invX = 1/X;
        invD = 1/D;
        A = invX*invD;
        B = invX*invX;
        C = invD*invD;*/
        //Method 2
        Xsq = (scalarProduct(mp, coord2, mp, coord2));
        Dsq = (scalarProduct(mp, coord3, mp, coord3));

        xd = scalarProduct(mp, coord2, mp, coord3);
        A = 1/sqrt(Xsq*Dsq);
        x = xd*A;
        B = x/Xsq;
        C = x/Dsq;
        //$$$

        if(abs(abs(x) - 1.0)<0.001)
            x = 0.999*x;

        if (x < -1.0) x = -1.0;
        else if (x > 1.0) x = 1.0;

        floatingpoint cosp =  x;
        posheta = 0.5*M_PI;
        floatingpoint sinp = sqrt(max<floatingpoint>((1-cosp*cosp),(floatingpoint)0.0));
        floatingpoint sinpminusq = sinp * cos(posheta) - cosp * sin(posheta);

        k = kpos[i] * sinpminusq/sinp;

	    /*if(abs(xd/XD - 1.0)<0.001){
		    xd = 0.999*XD;
	    }

        theta = safeacos(xd / XD);
        posheta = 0.5*M_PI;
        dTheta = theta-posheta;

        position = pos[i];

        k =  kpos[i] * A * sin(dTheta)/sin(theta);*/

	    //If the branching point is NOT bound to plusend, the local force variables
	    // ftemp1[3], ftemp2[3], ftemp3[3] represent forces on parent_minus, parent_plus
	    // and offspring_minus ends respectively.
	    // If the branching point IS bound to plusend, the force variables represent
	    // forces on parent_minus, parent_extendedplus and offspring_minus ends
	    // respectively. Under this condition, a transformation is necessary to realize
	    // the actual forces on the beads of interest.

	    r1x = coord2[0] - mp[0];
	    r1y = coord2[1] - mp[1];
        r1z = coord2[2] - mp[2];
        r2x = coord3[0] - mp[0];
        r2y = coord3[1] - mp[1];
        r2z = coord3[2] - mp[2];
        //Forces acting along the bond vectors, r1 and r2.
        Fr1x =  k * ( r2x*A - r1x*B );
        Fr1y =  k * ( r2y*A - r1y*B );
        Fr1z =  k * ( r2z*A - r1z*B );
        Fr2x =  k * ( r1x*A - r2x*C );
        Fr2y =  k * ( r1y*A - r2y*C );
        Fr2z =  k * ( r1z*A - r2z*C );
        //Force acting along the points mp, c2, and c3
        fmpx = -Fr1x-Fr2x;
        fmpy = -Fr1y-Fr2y;
        fmpz = -Fr1z-Fr2z;
        fc2x = Fr1x;
        fc2y = Fr1y;
        fc2z = Fr1z;
        fc3x = Fr2x;
        fc3y = Fr2y;
        fc3z = Fr2z;

        //Let us transform the forces from U(mp,c2,c3) to Utilde(c1,c2,c3)
        //bead1
        f1tempx = (1-position)*fmpx;
        f1tempy = (1-position)*fmpy;
        f1tempz = (1-position)*fmpz;

        //bead2
        f2tempx = fc2x + position*fmpx;
        f2tempy = fc2y + position*fmpy;
        f2tempz = fc2z + position*fmpz;

        //bead3
        f3tempx = fc3x;
        f3tempy = fc3y;
        f3tempz = fc3z;

        f3[0] += f3tempx;
        f3[1] += f3tempy;
        f3[2] += f3tempz;

        stretchforce[3*i] = f3tempx;
        stretchforce[3*i + 1] = f3tempy;
        stretchforce[3*i + 2] = f3tempz;

        //If you had calculated forces on the extended plus end, additional
        // transformations are needed.
        //Going from Utilda(c1, c2prime, c3) to U(c1, c2, c3)
        //c1-parent minusend | c2 - parent plusend/bindingsite | c2prime-parent
        // extendedplusend   | c3 - offspring minusend
        //[dU(c1,c2,c3)]      [dUtilda(c1,c2prime,c3)]   [dUtilda(c1,c2prime,c3)]   dc2prime
        //[------------]     =[----------------------] + [----------------------] x --------
        //[    dc1     ]c2,c3 [        dc1           ]   [        dc2prime      ]      dc1
        //__________________________________________________________________________________
        //[dU(c1,c2,c3)]       [dUtilda(c1,c2prime,c3)]   dc2prime
        //[------------]     = [----------------------] x --------
        //[    dc2     ]c1,c3  [        dc2prime      ]      dc2
        //__________________________________________________________________________________
        //We define c2 = c1 + s(c2prime - c1)
        // dc2prime    s - 1  | dc2prime    1
        // -------- = ------- | -------- = ---
        //   dc1         s    |   dc2       s
        if(areEqual(pos[i],(floatingpoint)1.0)){
            floatingpoint factor = (position-1)/position;
            f1[0] += f1tempx + f2tempx*factor;
            f1[1] += f1tempy + f2tempy*factor;
            f1[2] += f1tempz + f2tempz*factor;

            f2[0] += f2tempx*(1/position);
            f2[1] += f2tempy*(1/position);
            f2[2] += f2tempz*(1/position);
        }
        else{
            f1[0] += f1tempx;
            f1[1] += f1tempy;
            f1[2] += f1tempz;

            f2[0] += f2tempx;
            f2[1] += f2tempy;
            f2[2] += f2tempz;
        }

	    #ifdef CHECKFORCES_INF_NAN
	    if(checkNaN_INF<floatingpoint>(f1, 0, 2)||checkNaN_INF<floatingpoint>(f2,0,2)||checkNaN_INF<floatingpoint>(f3,0,2)){
		    cout<<"Branching Position Force becomes infinite. Printing data "<<endl;

            auto b = BranchingPoint::getBranchingPoints()[i];
            auto cyl1 = b->getFirstCylinder();
            auto cyl2 = b->getSecondCylinder();
            cout<<"Cylinder IDs "<<cyl1->getId()<<" "<<cyl2->getId()<<" with cIndex "
                <<cyl1->getStableIndex()<<" "<<cyl2->getStableIndex()<<" and bIndex "
                <<cyl1->getFirstBead()->getStableIndex()<<" "
                <<cyl1->getSecondBead()->getStableIndex()<<" "
                <<cyl2->getFirstBead()->getStableIndex()<<" "
                <<cyl2->getSecondBead()->getStableIndex()<<endl;
            cyl1->adjustedrelativeposition(pos[i], true);
            cout<<"Printing intermediate variables"<<endl;
            cout<<"xd="<<xd<<"A="<<A
                <<", B="<<B<<", C="<<C<<", cosp="<<cosp<<", sinp="<<sinp
                <<", sinpminusq="<<sinpminusq<<", position="<<position<<endl;

		    cout<<"Printing coords"<<endl;
		    cout<<coord1[0]<<" "<<coord1[1]<<" "<<coord1[2]<<endl;
		    cout<<coord2[0]<<" "<<coord2[1]<<" "<<coord2[2]<<endl;
		    cout<<coord3[0]<<" "<<coord3[1]<<" "<<coord3[2]<<endl;

		    cout<<"Printing force"<<endl;
		    cout<<f1[0]<<" "<<f1[1]<<" "<<f1[2]<<endl;
		    cout<<f2[0]<<" "<<f2[1]<<" "<<f2[2]<<endl;
		    cout<<f3[0]<<" "<<f3[1]<<" "<<f3[2]<<endl;

		    cout<<"Printing binary Coords"<<endl;
		    printvariablebinary(coord1,0,2);
		    printvariablebinary(coord2,0,2);
		    printvariablebinary(coord3,0,2);

		    cout<<"Printing binary Force"<<endl;
		    printvariablebinary(f1,0,2);
		    printvariablebinary(f2,0,2);
		    printvariablebinary(f3,0,2);

		    exit(EXIT_FAILURE);
	    }
	    #endif
    }
    delete[] mp;
    delete[] coord2prime;
}

} // namespace medyan

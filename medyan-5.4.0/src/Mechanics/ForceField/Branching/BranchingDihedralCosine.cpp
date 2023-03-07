
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

#include "BranchingDihedralCosine.h"
#include "BranchingDihedralCosineCUDA.h"
#include "BranchingDihedral.h"

#include "BranchingPoint.h"
#include "Bead.h"

#include "MathFunctions.h"
#include "Cylinder.h"
#ifdef CUDAACCL
#include <cuda.h>
#include <cuda_runtime.h>
#include "nvToolsExt.h"
#endif

namespace medyan {
using namespace mathfunc;
#ifdef CUDAACCL
void BranchingDihedralCosine::deallocate(){
    if(!(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamDestroy(stream));
    CUDAcommon::handleerror(cudaFree(gU_i));
    CUDAcommon::handleerror(cudaFree(gU_sum));
    CUDAcommon::handleerror(cudaFree(gFF));
    CUDAcommon::handleerror(cudaFree(ginteraction));
}
void BranchingDihedralCosine::optimalblocksnthreads( int nint){
    //CUDA stream create
    if(stream == NULL || !(CUDAcommon::getCUDAvars().conservestreams))
        CUDAcommon::handleerror(cudaStreamCreate(&stream));
    blocksnthreadse.clear();
    blocksnthreadsez.clear();
    blocksnthreadsf.clear();
    int blockSize;   // The launch configurator returned block size
    int minGridSize; // The minimum grid size needed to achieve the
    // maximum occupancy for a full device launch
//    int gridSize;    // The actual grid size needed, based on input size
//    unaryfn::argument_type blksize;
//    unaryfn::result_type result;
//    unaryfn ufn;
    if(nint>0) {
        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       BranchingDihedralCosineenergy, blockToSmem, 0);
//    std::cout<<(nint +blockSize -1) / blockSize<<" "<<blockSize<<endl;
//
//    cudaOccupancyMaxPotentialBlockSize( &minGridSize, &blockSize,
//                                        CUDAExclVolRepulsionenergy, 0, 0);
        blocksnthreadse.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadse.push_back(blockSize);
//    std::cout<<(nint +blockSize -1) / blockSize<<" "<<blockSize<<endl;
        blockSize = 0;

        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       BranchingDihedralCosineenergyz, blockToSmemez, 0);
        blocksnthreadsez.push_back((nint + blockSize - 1) / blockSize);
        blocksnthreadsez.push_back(blockSize);
        blockSize = 0;

        cudaOccupancyMaxPotentialBlockSizeVariableSMem(&minGridSize, &blockSize,
                                                       BranchingDihedralCosineforces, blockToSmem, 0);
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
        char b[] = "Branching Dihedral Cosine";
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
floatingpoint* BranchingDihedralCosine::energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                         floatingpoint *kdih, floatingpoint *pos, int *params) {
//    if(blocksnthreadse[1]>0) {
//
//        BranchingDihedralCosineenergy<<<blocksnthreadse[0], blocksnthreadse[1], (12 * blocksnthreadse[1]) * sizeof
//                (floatingpoint), stream>>>
//                          (coord, f, beadSet, kdih, pos, params, gU_i, CUDAcommon::getCUDAvars().gculpritID,
//                                  CUDAcommon::getCUDAvars().gculpritFF,
//                                  CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//        CUDAcommon::handleerror( cudaGetLastError(),"BranchingDihedralCosineenergy", "BranchingDihedralCosine.cu");
//        floatingpoint* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0,stream>>>(gU_i, params, gU_sum, gpu_Utot);
//        CUDAcommon::handleerror( cudaGetLastError() ,"BranchingDihedralCosineenergy", "BranchingDihedralCosine.cu");
//        return gU_sum;}
//    else
//        return NULL;
}


floatingpoint* BranchingDihedralCosine::energy(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                         floatingpoint *kdih, floatingpoint *pos, floatingpoint *z,
                                         int *params) {
        if(blocksnthreadse[1]>0) {

        BranchingDihedralCosineenergy<<<blocksnthreadse[0], blocksnthreadse[1], (12 * blocksnthreadse[1]) * sizeof
                (floatingpoint), stream>>>
                          (coord, f, beadSet, kdih, pos, params, gU_i, z, CUDAcommon::getCUDAvars().gculpritID,
                                  CUDAcommon::getCUDAvars().gculpritFF,
                                  CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction);
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//        CUDAcommon::handleerror( cudaGetLastError(),"BranchingDihedralCosineenergy", "BranchingDihedralCosine.cu");
//        floatingpoint* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0,stream>>>(gU_i, params, gU_sum, gpu_Utot);
        CUDAcommon::handleerror( cudaGetLastError() ,"BranchingDihedralCosineenergy", "BranchingDihedralCosine.cu");
//        return gU_sum;
        }

    if(blocksnthreadsez[1]>0) {
        BranchingDihedralCosineenergyz << < blocksnthreadsez[0], blocksnthreadsez[1], (24 * blocksnthreadsez[1]) *
                                            sizeof(floatingpoint), stream>> > (coord, f, beadSet, kdih, pos,
                                            params, gU_i, z, CUDAcommon::getCUDAvars().gculpritID,
                CUDAcommon::getCUDAvars().gculpritFF,
                CUDAcommon::getCUDAvars().gculpritinteraction, gFF, ginteraction );
//        auto cvars = CUDAcommon::getCUDAvars();
//        cvars.streamvec.push_back(&stream);
//        CUDAcommon::cudavars = cvars;
//        CUDAcommon::handleerror(cudaGetLastError(),"BranchingDihedralCosineenergyz", "BranchingDihedralCosine.cu");
//        floatingpoint* gpu_Utot = CUDAcommon::getCUDAvars().gpu_energy;
//        addvector<<<1,1,0,stream>>>(gU_i, params, gU_sum, gpu_Utot);
        CUDAcommon::handleerror(cudaGetLastError(),"BranchingDihedralCosineenergyz", "BranchingDihedralCosine.cu");

        return gU_sum;
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
void BranchingDihedralCosine::forces(floatingpoint *coord, floatingpoint *f, int *beadSet,
                                      floatingpoint *kdih,  floatingpoint *pos, int *params) {
    if (blocksnthreadsf[1] > 0) {
        BranchingDihedralCosineforces << < blocksnthreadsf[0], blocksnthreadsf[1], (12 * blocksnthreadsf[1]) *
                                            sizeof(floatingpoint), stream >> > (coord, f, beadSet, kdih, pos, params);
        auto cvars = CUDAcommon::getCUDAvars();
        cvars.streamvec.push_back(&stream);
        CUDAcommon::cudavars = cvars;
        CUDAcommon::handleerror(cudaGetLastError(),"BranchingDihedralCosineforces", "BranchingDihedralCosine.cu");
    }
}

void BranchingDihedralCosine::checkforculprit() {
    CUDAcommon::printculprit("BranchingDihedral","BranchingDihedralCosine");
    BranchingPoint *br;
    br = (BranchingPoint::getBranchingPoints()[CUDAcommon::getCUDAvars().culpritID[0]]);
    cout<<"Printing culprit branching point information."<<endl;
    br->printSelf();
    exit(EXIT_FAILURE);
}
#endif
//This version uses a the following vectors to determine dihedral angles.
//b1 = c2 - mp;
//b2 = c3 - mp;
//b3 = c4 - c3;
//n1 = b1 x b2;
//n2 = b3 x b2;
floatingpoint BranchingDihedralCosine::energy(
    const floatingpoint *coord, size_t nint,
    unsigned int *beadSet, floatingpoint *kdih, floatingpoint *pos
) const {

    int n = BranchingDihedral<BranchingDihedralCosine>::n;

    const floatingpoint *coord1, *coord2, *coord3, *coord4;
    floatingpoint position;
    floatingpoint n1n2, U_i;
    floatingpoint *mp = new floatingpoint[3];
    floatingpoint *n1 = new floatingpoint[3];
    floatingpoint *n2 = new floatingpoint[3];
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
        coord4 = &coord[beadSet[n * i + 3]];

        midPointCoordinate(mp, coord1, coord2, position);

        vectorProduct(n1, mp, coord2, mp, coord3);
        vectorProduct(n2, coord3, coord4, mp, coord3);

        normalizeVector(n1);
        normalizeVector(n2);
        n1n2 = dotProduct(n1, n2);

        U_i = kdih[i] * ( 1 - n1n2 );

        U += U_i;
    }
//	cout<<"Cosine      "<<U<<endl;
    delete [] mp;
    delete [] n1;
    delete [] n2;
    delete [] coord2prime;

    return U;
}



void BranchingDihedralCosine::forces(
    const floatingpoint *coord, floatingpoint *f, size_t nint,
    unsigned int *beadSet, floatingpoint *kdih, floatingpoint *pos, floatingpoint *stretchforce
) const {

    int n = BranchingDihedral<BranchingDihedralCosine>::n;

    double *coord1, *coord2, *coord3, *coord4;
    floatingpoint *f1, *f2, *f3, *f4;
    double force1[3], force2[3], force3[3], force4[3], position;

    coord1 = new double[3];
    coord2 = new double[3];
    coord3 = new double[3];
    coord4 = new double[3];

    double *mp = new double[3];
    double *n1 = new double[3];
    double *n2 = new double[3];
    double *zero = new double[3]; zero[0] = 0; zero[1] = 0; zero[2] = 0;

    //@{
    double Lforce1[3], Lforce2[3], Lforce3[3], Lforce4[3], Lforcemp[3], Lforce2prime[3];
    double vb1x, vb1y, vb1z, vb2x, vb2y, vb2z, vb3x, vb3y, vb3z;
    double vb2xm, vb2ym, vb2zm;
    double ax, ay, az, bx, by, bz;
    double rasq, rbsq, rgsq, rg, rginv, ra2inv, rb2inv, rabinv;
    double c,s;
    double p, df1, ddf1;
    double fg, fga, gaa, gbb, hg, hgb;
    double dtfx, dtfy, dtfz, dtgx, dtgy, dtgz, dthx, dthy, dthz, df, sx2, sy2,
            sz2;

    //  double *Nforce;
    //  Nforce = new double[12];
    //@}

    for(int i = 0; i < nint; i += 1) {

        for(int j = 0; j < 3; j++){
            coord1[j] = coord[beadSet[n * i]+ j];
            coord2[j] = coord[beadSet[n * i + 1]+ j];
            coord3[j] = coord[beadSet[n * i + 2]+ j];
            coord4[j] = coord[beadSet[n * i + 3]+ j];
        }

        //If the branching point is at the plus end of the cylinder, we end up with
        // singularities in energy and force expressions. To avoid it, we will create a
        // virtual plus end and use that to define vectors.
        position = pos[i];
        if(areEqual(position, 1.0)){
            position = 0.5;
            coord2[0] = (1/position)*(coord2[0] - (1-position)*coord1[0]);
            coord2[1] = (1/position)*(coord2[1] - (1-position)*coord1[1]);
            coord2[2] = (1/position)*(coord2[2] - (1-position)*coord1[2]);
        }

        midPointCoordinate<double>(mp, coord1, coord2, position);

        f1 = &f[beadSet[n * i]];
        f2 = &f[beadSet[n * i + 1]];
        f3 = &f[beadSet[n * i + 2]];
        f4 = &f[beadSet[n * i + 3]];

        //@ LAMMPS version test
        // 1st bond
        vb1x = coord2[0] - mp[0];
        vb1y = coord2[1] - mp[1];
        vb1z = coord2[2] - mp[2];

        // 2nd bond
        vb2x = coord3[0] - mp[0];
        vb2y = coord3[1] - mp[1];
        vb2z = coord3[2] - mp[2];

        vb2xm = -vb2x;
        vb2ym = -vb2y;
        vb2zm = -vb2z;

        // 3rd bond
        vb3x = coord4[0] - coord3[0];
        vb3y = coord4[1] - coord3[1];
        vb3z = coord4[2] - coord3[2];

        // c,s calculation

        ax = vb1y*vb2zm - vb1z*vb2ym;
        ay = vb1z*vb2xm - vb1x*vb2zm;
        az = vb1x*vb2ym - vb1y*vb2xm;
        bx = vb3y*vb2zm - vb3z*vb2ym;
        by = vb3z*vb2xm - vb3x*vb2zm;
        bz = vb3x*vb2ym - vb3y*vb2xm;

        //|b1x-b2|
        rasq = ax*ax + ay*ay + az*az;
        //|b3x-b2|
        rbsq = bx*bx + by*by + bz*bz;
        //|b2|^2
        rgsq = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm;
        //|b2|
        rg = sqrt(rgsq);

        rginv = ra2inv = rb2inv = 0.0;
        if (rg > 0) rginv = 1.0/rg;//1/|-b2|
        if (rasq > 0) ra2inv = 1.0/rasq;//1/|b1x-b2|
        if (rbsq > 0) rb2inv = 1.0/rbsq;//1/|b3x-b2|
        rabinv = sqrt(ra2inv*rb2inv);//1/|b1x-b2||b3x-b2|

        c = (ax*bx + ay*by + az*bz)*rabinv;//(b1x-b2).(b3x-b2)/|b1x-b2||b3x-b2|
        s = rg*rabinv*(ax*vb3x + ay*vb3y + az*vb3z);//|b2|((b1x-b2).b3)/|b1x-b2||b3x-b2|=cos<n1,b3>/sin<b1,b2>

        if (c > 1.0) c = 1.0;
        if (c < -1.0) c = -1.0;

        ddf1 = c;
        df1 = s;

        p = 1.0 - c;

        fg = vb1x*vb2xm + vb1y*vb2ym + vb1z*vb2zm;//b1.-b2
        hg = vb3x*vb2xm + vb3y*vb2ym + vb3z*vb2zm;//b3.-b2
        fga = fg*ra2inv*rginv;//(b1.-b2)/(|b1x-b2||-b2|)
        hgb = hg*rb2inv*rginv;//(b3.-b2)/(|b3x-b2||-b2|)
        gaa = -ra2inv*rg;//-|-b2|/|b1x-b2|
        gbb = rb2inv*rg;//|-b2|/|b3x-b2|

        dtfx = gaa*ax;//-|-b2|(b1x-b2_x)/|b1x-b2|=-|-b2|(n1cap_x)
        dtfy = gaa*ay;//-|-b2|(b1x-b2_y)/|b1x-b2|=-|-b2|(n1cap_y)
        dtfz = gaa*az;//-|-b2|(b1x-b2_z)/|b1x-b2|=-|-b2|(n1cap_z)
        dtgx = fga*ax - hgb*bx;//-|-b2||n1cap_x|-(b3.-b2)(b3x-b2_x)/(|b3x-b2||-b2|)=-|-b2||n1cap_x|-(b3.-b2)n2cap_x/|-b2|
        dtgy = fga*ay - hgb*by;//-|-b2||n1cap_y|-(b3.-b2)(b3x-b2_y)/(|b3x-b2||-b2|)=-|-b2||n1cap_y|-(b3.-b2)n2cap_y/|-b2|
        dtgz = fga*az - hgb*bz;//-|-b2||n1cap_z|-(b3.-b2)(b3x-b2_z)/(|b3x-b2||-b2|)=-|-b2||n1cap_z|-(b3.-b2)n2cap_z/|-b2|
        dthx = gbb*bx;//|-b2|(b3x-b2_x)/|b3x-b2|=|-b2|(n2cap_x)
        dthy = gbb*by;//|-b2|(b3x-b2_y)/|b3x-b2|=|-b2|(n2cap_y)
        dthz = gbb*bz;//|-b2|(b3x-b2_z)/|b3x-b2|=|-b2|(n2cap_z)

        df = -kdih[i] * df1;//-Kd.s

        sx2 = df*dtgx;//-Kd cos<n1,b3>/sin<b1,b2>(-|-b2||n1cap_x|-(b3.-b2)n2cap_x/|-b2|)
        sy2 = df*dtgy;
        sz2 = df*dtgz;

        Lforce4[0] = df*dthx;
        Lforce4[1] = df*dthy;
        Lforce4[2] = df*dthz;

        Lforce3[0] = -sx2 - Lforce4[0];
        Lforce3[1] = -sy2 - Lforce4[1];
        Lforce3[2] = -sz2 - Lforce4[2];

        //Special case where the branching point is plus end.
        if(areEqual(pos[i], (floatingpoint)1.0)){
            //In this scenario, the energy is defined as U1(c2, c2prime, c3, c4). We are
            // trying to get U(c1, c2, c3, c4) from it.
            //c1-parent minusend | c2 - parent plusend/bindingsite | c2prime-parent
            // extendedplusend   | c3 - offspring minusend         | c4 - offspring plusend
            //[dU(c1,c2,c3,c4)]             [dUtilda(c2,c2prime,c3,c4)]   dc2prime
            //[---------------]        =    [-------------------------] x --------
            //[    dc1        ]c2,c3,c4     [        dc2prime         ]     dc1
            //______________________________________________________________________________
            //[dU(c1,c2,c3,c4)]          [   [dUtilda(c2,c2prime,c3,c4)]   dc2prime
            //[---------------]        = [   [-------------------------] x --------
            //[    dc2        ]c1,c3,c4  [   [        dc2prime         ]     dc2
            //                           [
            //                           [        [dUtilda(c2,c2prime,c3,c4)]
            //                           [  +     [-------------------------]
            //                           [        [           dc2           ]
            //We define c2 = c1 + s(c2prime - c1)
            // dc2prime    s - 1  | dc2prime    1
            // -------- = ------- | -------- = ---
            //   dc1         s    |   dc2       s

            Lforce2prime[0] = df*dtfx;
            Lforce2prime[1] = df*dtfy;
            Lforce2prime[2] = df*dtfz;

            Lforce2[0] = sx2 - Lforce2prime[0];
            Lforce2[1] = sy2 - Lforce2prime[1];
            Lforce2[2] = sz2 - Lforce2prime[2];

            //Transformation
            double factor = (position-1)/position;
            Lforce1[0] = Lforce2prime[0]*factor;
            Lforce1[1] = Lforce2prime[1]*factor;
            Lforce1[2] = Lforce2prime[2]*factor;

            factor = (1/position);
            Lforce2[0] += Lforce2prime[0]*factor;
            Lforce2[1] += Lforce2prime[1]*factor;
            Lforce2[2] += Lforce2prime[2]*factor;

        }
        //Default case.
        else {
            //In this scenario, the energy is defined as U2(mp, c2, c3, c4). We are trying
            // to get U(c1, c2, c3, c4) from it.
            Lforce2[0] = df*dtfx;
            Lforce2[1] = df*dtfy;
            Lforce2[2] = df*dtfz;

            Lforcemp[0] = sx2 - Lforce2[0];
            Lforcemp[1] = sy2 - Lforce2[1];
            Lforcemp[2] = sz2 - Lforce2[2];

            //Transformation
            Lforce1[0] = (1 - position) * Lforcemp[0];
            Lforce1[1] = (1 - position) * Lforcemp[1];
            Lforce1[2] = (1 - position) * Lforcemp[2];

            Lforce2[0] += (position) * Lforcemp[0];
            Lforce2[1] += (position) * Lforcemp[1];
            Lforce2[2] += (position) * Lforcemp[2];
        }

        //Add to force vector
        f1[0] += Lforce1[0];
        f1[1] += Lforce1[1];
        f1[2] += Lforce1[2];

        f2[0] += Lforce2[0];
        f2[1] += Lforce2[1];
        f2[2] += Lforce2[2];

        f3[0] += Lforce3[0];
        f3[1] += Lforce3[1];
        f3[2] += Lforce3[2];

        f4[0] += Lforce4[0];
        f4[1] += Lforce4[1];
        f4[2] += Lforce4[2];

        if(stretchforce) {
            stretchforce[3*i] = Lforce3[0];
            stretchforce[3*i + 1] = Lforce3[1];
            stretchforce[3*i + 2] = Lforce3[2];
        }

        //@
        #ifdef CHECKFORCES_INF_NAN
        if(checkNaN_INF<floatingpoint>(f1, 0, 2)||checkNaN_INF<floatingpoint>(f2,0,2)
        ||checkNaN_INF<floatingpoint>(f3,0,2) ||checkNaN_INF<floatingpoint>(f4,0,2)){
            cout<<"Branching Dihedral Force becomes infinite. Printing data "<<endl;

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

            cout<<"Parent filament binding fraction "<<position<<endl;
            cout<<"Parent filament binding position "<<mp[0]<<" "<<mp[1]<<" "<<mp[2]<<endl;
            cout<<"ax-bz"<<ax<<" "<<ay<<" "<<az<<" "<<bx<<" "<<by<<" "<<bz<<endl;
            cout<<"rasq "<<rasq<<" rbsq "<<rbsq<<" rgsq "<<rgsq<<" rg "<<rg<<endl;
            cout<<"rginv "<<rginv<<" ra2inv "<<ra2inv<<" rb2inv "<<rb2inv<<endl;
            cout<<"c "<<c<<" s "<<s<<endl;
            cout<<"fg "<<fg<<" hg "<<hg<<" fga "<<fga<<" hgb "<<hgb<<" gaa "<<gaa<<" gbb "
                                                                                    ""<<gbb<<endl;
            cout<<"dtfx-z "<<dtfx<<" "<<dtfy<<" "<<dtfz<<endl;
            cout<<"dtgx-z "<<dtgx<<" "<<dtgy<<" "<<dtgz<<endl;
            cout<<"dthx-z "<<dthx<<" "<<dthy<<" "<<dthz<<endl;
            cout<<"df "<<df<<endl;
            cout<<"s(x-z)2 "<<sx2<<" "<<sy2<<" "<<sz2<<endl;

            cout<<"Printing coords"<<endl;
            cout<<coord1[0]<<" "<<coord1[1]<<" "<<coord1[2]<<endl;
            cout<<coord2[0]<<" "<<coord2[1]<<" "<<coord2[2]<<endl;
            cout<<coord3[0]<<" "<<coord3[1]<<" "<<coord3[2]<<endl;
            cout<<coord4[0]<<" "<<coord4[1]<<" "<<coord4[2]<<endl;
            cout<<"Printing force"<<endl;
            cout<<f1[0]<<" "<<f1[1]<<" "<<f1[2]<<endl;
            cout<<f2[0]<<" "<<f2[1]<<" "<<f2[2]<<endl;
            cout<<f3[0]<<" "<<f3[1]<<" "<<f3[2]<<endl;
            cout<<f4[0]<<" "<<f4[1]<<" "<<f4[2]<<endl;
            cout<<"Printing binary Coords"<<endl;
            printvariablebinary(coord1,0,2);
            printvariablebinary(coord2,0,2);
            printvariablebinary(coord3,0,2);
            printvariablebinary(coord4,0,2);
            cout<<"Printing binary Force"<<endl;
            printvariablebinary(f1,0,2);
            printvariablebinary(f2,0,2);
            printvariablebinary(f3,0,2);
            printvariablebinary(f4,0,2);
            exit(EXIT_FAILURE);
        }
        #endif

        /*forcesNumericalinteraction(coord, f, nint, beadSet, kdih, pos,i, Nforce);
        //@ Calculate EmagerrorLvsN
            double Mmags[4];
            Mmags[0] = sqrt(force1[0]*force1[0] + force1[1]*force1[1]+force1[2]*force1[2]);
            Mmags[1] = sqrt(force2[0]*force2[0] + force2[1]*force2[1]+force2[2]*force2[2]);
            Mmags[2] = sqrt(force3[0]*force3[0] + force3[1]*force3[1]+force3[2]*force3[2]);
            Mmags[3] = sqrt(force4[0]*force4[0] + force4[1]*force4[1]+force4[2]*force4[2]);

            double Lmags[4];
            Lmags[0] = sqrt(Lforce1[0]*Lforce1[0] + Lforce1[1]*Lforce1[1]+Lforce1[2]*Lforce1[2]);
            Lmags[1] = sqrt(Lforce2[0]*Lforce2[0] + Lforce2[1]*Lforce2[1]+Lforce2[2]*Lforce2[2]);
            Lmags[2] = sqrt(Lforce3[0]*Lforce3[0] + Lforce3[1]*Lforce3[1]+Lforce3[2]*Lforce3[2]);
            Lmags[3] = sqrt(Lforce4[0]*Lforce4[0] + Lforce4[1]*Lforce4[1]+Lforce4[2]*Lforce4[2]);


            double Nmags[4];
            Nmags[0] = sqrt(Nforce[0]*Nforce[0] + Nforce[1]*Nforce[1]+Nforce[2]*Nforce[2]);
            Nmags[1] = sqrt(Nforce[3]*Nforce[3] + Nforce[4]*Nforce[4]+Nforce[5]*Nforce[5]);
            Nmags[2] = sqrt(Nforce[6]*Nforce[6] + Nforce[7]*Nforce[7]+Nforce[8]*Nforce[8]);
            Nmags[3] = sqrt(Nforce[9]*Nforce[9] + Nforce[10]*Nforce[10]+Nforce[11]*Nforce[11]);

            //Normalize vectors
            for(int axis = 0; axis <3; axis++){
            force1[axis] = force1[axis]/Mmags[0];
            force2[axis] = force2[axis]/Mmags[1];
            force3[axis] = force3[axis]/Mmags[2];
            force4[axis] = force4[axis]/Mmags[3];
            Lforce1[axis] = Lforce1[axis]/Lmags[0];
            Lforce2[axis] = Lforce2[axis]/Lmags[1];
            Lforce3[axis] = Lforce3[axis]/Lmags[2];
            Lforce4[axis] = Lforce4[axis]/Lmags[3];
            for(int fvec =0; fvec < 4; fvec++){
                Nforce[3*fvec + axis] = Nforce[3*fvec + axis]/Nmags[fvec];
            }
            }

            //Calculate dotProduct
            FdirerrorMvsN[0] += abs(1 - (force1[0]*Nforce[0] + force1[1]*Nforce[1] +
                force1[2]*Nforce[2]));
            FdirerrorMvsN[1] += abs(1 - (force2[0]*Nforce[3] + force2[1]*Nforce[4] +
                force2[2]*Nforce[5]));
            FdirerrorMvsN[2] += abs(1 - (force3[0]*Nforce[6] + force3[1]*Nforce[7] +
                force3[2]*Nforce[8]));
            FdirerrorMvsN[3] += abs(1 - (force4[0]*Nforce[9] + force4[1]*Nforce[10] +
                force4[2]*Nforce[11]));

            double t1, t2, t3, t4;
            t1 = (Lforce1[0]*Nforce[0] + Lforce1[1]*Nforce[1] + Lforce1[2]*Nforce[2]);
            t2 = Lforce2[0]*Nforce[3] + Lforce2[1]*Nforce[4] + Lforce2[2]*Nforce[5];
            t3 = Lforce3[0]*Nforce[6] + Lforce3[1]*Nforce[7] + Lforce3[2]*Nforce[8];
            t4 = Lforce4[0]*Nforce[9] + Lforce4[1]*Nforce[10] + Lforce4[2]*Nforce[11];
    //          cout<<t1<<" "<<t2<<" "<<t3<<" "<<t4<<endl;

            FdirerrorLvsN[0] += abs(1 - t1);
            FdirerrorLvsN[1] += abs(1 - t2);
            FdirerrorLvsN[2] += abs(1 - t3);
            FdirerrorLvsN[3] += abs(1 - t4);

            for(int comp =0; comp < 4; comp++){
            FmagerrorMvsN[comp] += abs(Nmags[comp] - Mmags[comp])/Nmags[comp];
            FmagerrorLvsN[comp] += abs(Nmags[comp] - Lmags[comp])/Nmags[comp];
            }
        counterF++;
        //@}

    }
    if(counterF>0){
    cout<<"----------------"<<endl;
    if(true) {
        cout << "Forces MAG relative error counter " << counterF << endl;
        cout << "MEDYAN vs LAMMPS" << endl;
        for (int comp = 0; comp < 4; comp++) {
            cout << "f" << comp + 1 << " " << FmagerrorMvsN[comp] / counterF << " "
                    << FmagerrorLvsN[comp] / counterF << endl;
        }
    }
    if(false) {
        cout << "Forces DIR relative error counter " << counterF << endl;
        cout << "MEDYAN vs LAMMPS" << endl;
        for (int comp = 0; comp < 4; comp++) {
            cout << "f" << comp << " " << FdirerrorMvsN[comp] / counterF << " "
                    << FdirerrorLvsN[comp] / counterF << endl;
        }
    }*/
    }

    delete [] coord1;
    delete [] coord2;
    delete [] coord3;
    delete [] coord4;

    delete [] mp;
    delete [] n1;
    delete [] n2;
    delete [] zero;

    //    delete [] Nforce;

}

void BranchingDihedralCosine::forcesNumericalinteraction(
		floatingpoint *coord, floatingpoint *f, size_t nint,
		unsigned int *beadSet, floatingpoint *kdih, floatingpoint *pos, int interactionID, double *Nforce)
		{
    int i = interactionID;

    int n = BranchingDihedral<BranchingDihedralCosine>::n;

    double delta = 0.001;

    double Ucurr = energyininteractionperturbed(coord, nint, beadSet,
		                                                   kdih, pos, i, -1, -1, delta);
//		cout<<i<<" Ucurr "<<Ucurr<<endl;

		//perturb x y z of each of the 4 coordinates
		//Calculate energy
    for(int perturbcoord = 0; perturbcoord<4; perturbcoord++){
        for(int perturbaxis = 0;perturbaxis<3;perturbaxis++){
            double Upert = energyininteractionperturbed(coord, nint, beadSet,
                    kdih, pos, i, perturbcoord, perturbaxis, delta);

//				cout<<i<<" Upert "<<Upert<<endl;
            //Calculate slope
            double slope = -(Upert - Ucurr)/delta;
            // Assign as force
            Nforce[3*perturbcoord + perturbaxis] = slope;
        }
    }
    /*cout<<"forcesN ";
    for(int i = 0; i < 12;i++)
        cout<<Nforce[i]<<" ";
    cout<<endl;*/
}

template <class dataType>
dataType BranchingDihedralCosine::energyininteractionperturbed(
		floatingpoint *coord, size_t nint,
		unsigned int *beadSet, floatingpoint *kdih, floatingpoint *pos, int
		interactionID, int perturbcoord, int perturbaxis, double delta){

    int i = interactionID;

    int n = BranchingDihedral<BranchingDihedralCosine>::n;

    dataType *coord1, *coord2, *coord3, *coord4, n1n2, U_i;
    coord1 = new dataType[3];
    coord2 = new dataType[3];
    coord3 = new dataType[3];
    coord4 = new dataType[3];
    dataType *mp = new dataType[3];
    dataType *n1 = new dataType[3];
    dataType *n2 = new dataType[3];

    //@{
    dataType vb1x, vb1y, vb1z, vb2x, vb2y, vb2z, vb3x, vb3y, vb3z;
    dataType vb2xm, vb2ym, vb2zm;
    dataType ax, ay, az, bx, by, bz;
    dataType rasq, rbsq, rgsq, rg, rginv, ra2inv, rb2inv, rabinv;
    dataType c,s;
    dataType p, df1, ddf1;
    dataType fg, fga, gaa, gbb, hg, hgb;
    //@}

    for(int j = 0; j < 3; j++){
        coord1[j] = coord[beadSet[n * i]+ j];
        coord2[j] = coord[beadSet[n * i + 1]+ j];
        coord3[j] = coord[beadSet[n * i + 2]+ j];
        coord4[j] = coord[beadSet[n * i + 3]+ j];
    }

    if(perturbcoord != -1 && perturbaxis != -1){
        if(perturbcoord == 0){
            if(perturbaxis == 0)
                coord1[0] += delta;
            else if(perturbaxis == 1)
                coord1[1] += delta;
            else
                coord1[2] += delta;
        }
        else if(perturbcoord == 1){
            if(perturbaxis == 0)
                coord2[0] += delta;
            else if(perturbaxis == 1)
                coord2[1] += delta;
            else
                coord2[2] += delta;
        }
        else if(perturbcoord == 2){
            if(perturbaxis == 0)
                coord3[0] += delta;
            else if(perturbaxis == 1)
                coord3[1] += delta;
            else
                coord3[2] += delta;
        }
        else if(perturbcoord == 3){
            if(perturbaxis == 0)
                coord4[0] += delta;
            else if(perturbaxis == 1)
                coord4[1] += delta;
            else
                coord4[2] += delta;
        }
    }

    dataType alpha = pos[i];
    dataType beta = 1 - alpha;
    mp[0] = (coord1[0] * beta + alpha * coord2[0]);
    mp[1] = (coord1[1] * beta + alpha * coord2[1]);
    mp[2] = (coord1[2] * beta + alpha * coord2[2]);

    //@ LAMMPS version test
    // 1st bond
    vb1x = coord2[0] - mp[0];
    vb1y = coord2[1] - mp[1];
    vb1z = coord2[2] - mp[2];

    // 2nd bond
    vb2x = coord3[0] - mp[0];
    vb2y = coord3[1] - mp[1];
    vb2z = coord3[2] - mp[2];

    vb2xm = -vb2x;
    vb2ym = -vb2y;
    vb2zm = -vb2z;

    // 3rd bond
    vb3x = coord4[0] - coord3[0];
    vb3y = coord4[1] - coord3[1];
    vb3z = coord4[2] - coord3[2];

    // c,s calculation

    ax = vb1y*vb2zm - vb1z*vb2ym;
    ay = vb1z*vb2xm - vb1x*vb2zm;
    az = vb1x*vb2ym - vb1y*vb2xm;
    bx = vb3y*vb2zm - vb3z*vb2ym;
    by = vb3z*vb2xm - vb3x*vb2zm;
    bz = vb3x*vb2ym - vb3y*vb2xm;

    rasq = ax*ax + ay*ay + az*az;
    rbsq = bx*bx + by*by + bz*bz;
    rgsq = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm;
    rg = sqrt(rgsq);

    rginv = ra2inv = rb2inv = 0.0;
    if (rg > 0) rginv = 1.0/rg;
    if (rasq > 0) ra2inv = 1.0/rasq;
    if (rbsq > 0) rb2inv = 1.0/rbsq;
    rabinv = sqrt(ra2inv*rb2inv);

    c = (ax*bx + ay*by + az*bz)*rabinv;
    s = rg*rabinv*(ax*vb3x + ay*vb3y + az*vb3z);

    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    ddf1 = c;
    df1 = s;

    p = 1.0 - c;
    U_i = kdih[i] * p;

    delete [] coord1;
    delete [] coord2;
    delete [] coord3;
    delete [] coord4;

    delete [] mp;
    delete [] n1;
    delete [] n2;

    return U_i;
}

void BranchingDihedralCosine::testdihedral(){

	vector<floatingpoint> coord ={100, 100, 100, 100, 100, 300, 100, 200, 200, 80, 200,
							   300};
	vector<floatingpoint> force={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
	vector<unsigned int> beadSet ={0, 3, 6, 9};
	vector<floatingpoint> kdih ={50.0};
	vector<floatingpoint> pos = {0.5};
	vector<floatingpoint> stretchForce ={0.0, 0.0, 0.0};
	int nint = 1;
	cout<<"---Test begins----"<<endl;
	floatingpoint U = energy(coord.data(), nint, beadSet.data(), kdih.data(), pos.data());
	forces(coord.data(), force.data(), nint, beadSet.data(), kdih.data(), pos.data(),
			stretchForce.data());
	cout<<"---Test ends----"<<endl;
}

} // namespace medyan

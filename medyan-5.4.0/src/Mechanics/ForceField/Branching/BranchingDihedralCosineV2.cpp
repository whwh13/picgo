
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

#include "BranchingDihedralCosineV2.h"
#include "BranchingDihedral.h"

#include "BranchingPoint.h"
#include "Bead.h"

#include "MathFunctions.h"
#include "Cylinder.h"

namespace medyan {
using namespace mathfunc;

//This version uses a the following vectors to determine dihedral angles.
// b1 = c4 - c3;
// b2 = c2 - mp;
// b3 = c3 - mp;
// n1 = b1 x b2;
// n2 = b3 x b2;

//Forces calculated are on points c5,c2,mp,c3. c5 = c4-c3+c2;

//STEP1: Transformation of forces on c5, c2, mp, and c3 TO forces on c2, mp, c3 and c4.
// Fc2 = -dE/dc2 + -dE/dc5* dc5/dc2 = -dE/dc2 + -dE/dc5 = Fc2 + Fc5
// F_mp = Fmp
// Fc3  = -dE/dc3 + -dE/dc5 * dc5/dc3 = Fc3 -Fc5
// Fc4  = -dE/dc5 * dc5/dc4 = Fc5

// STEP2: Transofmration of forces from c2, mp, c3, c4 TO c1, c2, c3, c4
// Fc1 = (1-p) * Fmp
// Fc2 = Fc2 + p * Fmp
// Fc3 = Fc3
// Fc4 = Fc4

floatingpoint BranchingDihedralCosineV2::energy(
		floatingpoint *coord, size_t nint,
		unsigned int *beadSet, floatingpoint *kdih, floatingpoint *pos){

    int n = BranchingDihedral<BranchingDihedralCosineV2>::n;

    floatingpoint *coord1, *coord2, *coord3, *coord4, n1n2, U_i, position;
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
        if(position == 1){
            position = 0.5;//assign a dummy position and extend plus end to new coordinates.
            for(int dim = 0; dim < 3; dim++) {
                coord2prime[dim] = (1 / position) * (coord2[dim] - (1 - position) *
                                                                   coord1[dim]);//extended plus end coordinate;
            }
            coord2 = &coord2prime[0];
        }
        coord3 = &coord[beadSet[n * i + 2]];
        coord4 = &coord[beadSet[n * i + 3]];

        midPointCoordinate(mp, coord1, coord2, pos[i]);

        // b1 = c4 - c3;
        // b2 = c2 - mp;
        // b3 = c3 - mp;
        // n1 = b1 x b2;
        // n2 = b3 x b2;
        vectorProduct(n1, coord3, coord4, mp, coord2);
        vectorProduct(n2, mp, coord3, mp, coord2);

        normalizeVector(n1);
        normalizeVector(n2);
        n1n2 = dotProduct(n1, n2);

        U_i = kdih[i] * ( 1 - n1n2 );

        U += U_i;
    }

    //	cout<<"CosineV2     "<<U<<endl;
    delete [] mp;
    delete [] n1;
    delete [] n2;
    delete [] coord2prime;
    return U;
}

void BranchingDihedralCosineV2::forces(
		floatingpoint *coord, floatingpoint *f, size_t nint,
		unsigned int *beadSet, floatingpoint *kdih, floatingpoint *pos, floatingpoint *stretchforce){

    int n = BranchingDihedral<BranchingDihedralCosineV2>::n;

    double *coord1, *coord2, *coord3, *coord4;
    floatingpoint *f1, *f2, *f3, *f4;

    coord1 = new double[3];
    coord2 = new double[3];
    coord3 = new double[3];
    coord4 = new double[3];

    double *mp = new double[3];
    double *n1 = new double[3];
    double *n2 = new double[3];
    double *zero = new double[3]; zero[0] = 0; zero[1] = 0; zero[2] = 0;

    //@{
    double Fc1[3], Fc2[3], Fc3[3], Fc4[3], Fc5[3], Fmp[3], Fc2prime[3], position;
    double vb1x, vb1y, vb1z, vb2x, vb2y, vb2z, vb3x, vb3y, vb3z;
    double vb2xm, vb2ym, vb2zm;
    double ax, ay, az, bx, by, bz;
    double rasq, rbsq, rgsq, rg, rginv, ra2inv, rb2inv, rabinv;
    double c,s;
    double p, df1, ddf1;
    double fg, fga, gaa, gbb, hg, hgb;
    double dtfx, dtfy, dtfz, dtgx, dtgy, dtgz, dthx, dthy, dthz, df, sx2, sy2,
            sz2;

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

        // b1 = c4 - c3;
        // b2 = c2 - mp;
        // b3 = c3 - mp;
        // n1 = b1 x b2;
        // n2 = b3 x b2;

        //@ LAMMPS version test
        // 1st bond
        vb1x = coord4[0] - coord3[0];
        vb1y = coord4[1] - coord3[1];
        vb1z = coord4[2] - coord3[2];

        // 2nd bond
        vb2x = coord2[0] - mp[0];
        vb2y = coord2[1] - mp[1];
        vb2z = coord2[2] - mp[2];

        vb2xm = vb2x;
        vb2ym = vb2y;
        vb2zm = vb2z;

        // 3rd bond
        vb3x = coord3[0] - mp[0];
        vb3y = coord3[1] - mp[1];
        vb3z = coord3[2] - mp[2];

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

        if(areEqual(pos[i], floatingpoint(1.0))) {
            //Forces calculated are on points c5,c2prime,c2,c3. c5 = c4-c3+c2;

            //I II III IV
            //c5 c2prime c2 c3

            Fc2prime[0] = df*dtfx;
            Fc2prime[1] = df*dtfy;
            Fc2prime[2] = df*dtfz;

            Fc5[0] = sx2 - Fc2prime[0];
            Fc5[1] = sy2 - Fc2prime[1];
            Fc5[2] = sz2 - Fc2prime[2];

            Fc3[0] = df*dthx;
            Fc3[1] = df*dthy;
            Fc3[2] = df*dthz;

            Fc2[0] = -sx2 - Fc3[0];
            Fc2[1] = -sy2 - Fc3[1];
            Fc2[2] = -sz2 - Fc3[2];

            //STEP1: Transformation of forces on c5, c2prime, c2, and c3 TO forces on
            // c2prime, c2, c3 and c4.
            // Fc2prime = Fc2prime + Fc5
            // Fc2 = Fc2
            // Fc3  = Fc3 - Fc5
            // Fc4  = Fc5

            Fc2prime[0] += Fc5[0];
            Fc2prime[1] += Fc5[1];
            Fc2prime[2] += Fc5[2];

            Fc3[0] += -Fc5[0];
            Fc3[1] += -Fc5[1];
            Fc3[2] += -Fc5[2];

            Fc4[0] = Fc5[0];
            Fc4[1] = Fc5[1];
            Fc4[2] = Fc5[2];

            // STEP2: Transofmration of forces from c2prime, c2, c3, c4 TO c1, c2, c3, c4
            // Fc1 = Fc2prime*(p-1)/p
            // Fc2 = Fc2prime*(1/p) + Fc2
            // Fc3 = Fc3
            // Fc4 = Fc4
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
            double factor = (position - 1)/position;
            Fc1[0] = factor*Fc2prime[0];
            Fc1[1] = factor*Fc2prime[1];
            Fc1[2] = factor*Fc2prime[2];

            factor = (1/position);
            Fc2[0] += factor*Fc2prime[0];
            Fc2[1] += factor*Fc2prime[1];
            Fc2[2] += factor*Fc2prime[2];
        }
        //Default case
        else{
            //Forces calculated are on points c5,c2,mp,c3. c5 = c4-c3+c2;

            //I II III IV
            //c5 c2 mp c3

            Fc2[0] = df*dtfx;
            Fc2[1] = df*dtfy;
            Fc2[2] = df*dtfz;

            Fc5[0] = sx2 - Fc2[0];
            Fc5[1] = sy2 - Fc2[1];
            Fc5[2] = sz2 - Fc2[2];

            Fc3[0] = df*dthx;
            Fc3[1] = df*dthy;
            Fc3[2] = df*dthz;

            Fmp[0] = -sx2 - Fc3[0];
            Fmp[1] = -sy2 - Fc3[1];
            Fmp[2] = -sz2 - Fc3[2];

            //STEP1: Transformation of forces on c5, c2, mp, and c3 TO forces on c2, mp, c3 and c4.
            // Fc2 = Fc2 + Fc5
            // F_mp = Fmp
            // Fc3  = Fc3 - Fc5
            // Fc4  = Fc5

            Fc2[0] += Fc5[0];
            Fc2[1] += Fc5[1];
            Fc2[2] += Fc5[2];

            Fc3[0] += -Fc5[0];
            Fc3[1] += -Fc5[1];
            Fc3[2] += -Fc5[2];

            Fc4[0] = Fc5[0];
            Fc4[1] = Fc5[1];
            Fc4[2] = Fc5[2];

            // STEP2: Transofmration of forces from c2, mp, c3, c4 TO c1, c2, c3, c4
            // Fc1 = (1-p) * Fmp
            // Fc2 = Fc2 + p * Fmp
            // Fc3 = Fc3
            // Fc4 = Fc4

            Fc1[0] = (1 - position) * Fmp[0];
            Fc1[1] = (1 - position) * Fmp[1];
            Fc1[2] = (1 - position) * Fmp[2];

            Fc2[0] += (position) * Fmp[0];
            Fc2[1] += (position) * Fmp[1];
            Fc2[2] += (position) * Fmp[2];
            Fc2[2] += (position) * Fmp[2];
        }

        //Add to force vector
        f1[0] += Fc1[0];
        f1[1] += Fc1[1];
        f1[2] += Fc1[2];

        f2[0] += Fc2[0];
        f2[1] += Fc2[1];
        f2[2] += Fc2[2];

        f3[0] += Fc3[0];
        f3[1] += Fc3[1];
        f3[2] += Fc3[2];

        f4[0] += Fc4[0];
        f4[1] += Fc4[1];
        f4[2] += Fc4[2];

        stretchforce[3*i] = Fc3[0];
        stretchforce[3*i + 1] = Fc3[1];
        stretchforce[3*i + 2] = Fc3[2];
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
    }

    delete [] coord1;
    delete [] coord2;
    delete [] coord3;
    delete [] coord4;

    delete [] mp;
    delete [] n1;
    delete [] n2;
    delete [] zero;

}

} // namespace medyan

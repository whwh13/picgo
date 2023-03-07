
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

#include "FilamentStretchingHarmonicandBendingHarmonic.h"
#include "FilamentStretchingandBending.h"

#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"
#include "MathFunctions.h"


namespace medyan {
using namespace mathfunc;

void FilamentStretchingHarmonicandBendingHarmonic::energy(floatingpoint *coord,
                                                        std::size_t nint, int *beadSet, floatingpoint *kstr, floatingpoint *kbend,
                                                        floatingpoint *eql, floatingpoint *eqt, floatingpoint* totalenergy, const int
                                                        startID, const int endID, int threadID){

    int n = FilamentStretchingandBending<FilamentStretchingHarmonicandBendingHarmonic>::n;

    floatingpoint *coord1, *coord2, *coord3, L1, L2, L1L2, l1l2;

    floatingpoint Ubend = 0.0, Ustr = 0.0, U_ibend, U_istr, dist;

    for(int i = startID; i < endID; i += 1) {

        coord1 = &coord[beadSet[n * i]];
        coord2 = &coord[beadSet[n * i + 1]];
        coord3 = &coord[beadSet[n * i + 2]];

        L1 = sqrt(scalarProduct(coord1, coord2,
                                coord1, coord2));
        L2 = sqrt(scalarProduct(coord2, coord3,
                                coord2, coord3));

        //calculate stretching energy between coord2 and coord3
        dist = L2 - eql[i];
        U_istr = 0.5 * kstr[i] * dist * dist;

        L1L2 = L1*L2;
        l1l2 = scalarProduct(coord1, coord2,
                                coord2, coord3);

        floatingpoint x = l1l2/L1L2;

        if (x < -1.0) x = -1.0;
        else if (x > 1.0) x = 1.0;

        U_ibend = kbend[i] * (1 - x);

        Ubend += U_ibend;
        Ustr += U_istr;
    }

    totalenergy[3*threadID] = Ustr;
    totalenergy[3*threadID + 2] = Ubend;
}

void FilamentStretchingHarmonicandBendingHarmonic::energy(floatingpoint *coord, int *beadSetsansbending,
                                                        floatingpoint *kstrsansbending, floatingpoint *eqlsansbending,
                                                        floatingpoint* totalenergy, const int startID, const int endID, int
                                                        threadID){

    int n = 2;

    floatingpoint *coord1, *coord2, dist;

    floatingpoint U_i, U = 0.0;

    for(int i = startID; i < endID; i += 1) {

        coord1 = &coord[beadSetsansbending[n * i]];
        coord2 = &coord[beadSetsansbending[n * i + 1]];
        dist = twoPointDistance(coord1, coord2) - eqlsansbending[i];

        U_i = 0.5 * kstrsansbending[i] * dist * dist;

        U += U_i;
    }

    totalenergy[3*threadID + 1] = U;
    return;
}

void FilamentStretchingHarmonicandBendingHarmonic::forces(floatingpoint *coord, floatingpoint *f, size_t nint, int *beadSet,
                                                          floatingpoint *kbend,
                                                          floatingpoint *kstr,  floatingpoint *eql, floatingpoint *eqt){

    int n = FilamentStretchingandBending<FilamentStretchingHarmonicandBendingHarmonic>::n;

    floatingpoint *coord1, *coord2, *coord3, L1, L2, l1l2, invL1, invL2, A,B,C, k;
    floatingpoint *force1, *force2, *force3;

    for(int i = 0; i < nint; i += 1) {

        coord1 = &coord[beadSet[n * i]];
        coord2 = &coord[beadSet[n * i + 1]];
        coord3 = &coord[beadSet[n * i + 2]];

        force1 = &f[beadSet[n * i]];
        force2 = &f[beadSet[n * i + 1]];
        force3 = &f[beadSet[n * i + 2]];

        L1 = sqrt(scalarProduct(coord1, coord2, coord1, coord2));
        L2 = sqrt(scalarProduct(coord2, coord3, coord2, coord3));

        l1l2 = scalarProduct(coord1, coord2, coord2, coord3);

        invL1 = 1/L1;
        invL2 = 1/L2;
        A = invL1*invL2;
        B = l1l2*invL1*A*A*L2;
        C = l1l2*invL2*A*A*L1;

        if (areEqual(eqt[i], 0.0)) k = kbend[i];

        else{
            if(abs(abs(l1l2*A) - 1.0)<0.001)
                l1l2 = 0.999*l1l2;

            floatingpoint x = l1l2 *A;
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
        //force on i-1, f = k*(-A*l2 + B*l1):
        force1[0] +=  k * ((-coord3[0] + coord2[0])*A +
                            (coord2[0] - coord1[0])*B );
        force1[1] +=  k * ((-coord3[1] + coord2[1])*A +
                            (coord2[1] - coord1[1])*B );
        force1[2] +=  k * ((-coord3[2] + coord2[2])*A +
                            (coord2[2] - coord1[2])*B );


        //force on i, f = k*(A*(l1-l2) - B*l1 + C*l2):
        force2[0] +=  k *( (coord3[0] - 2*coord2[0] + coord1[0])*A -
                            (coord2[0] - coord1[0])*B +
                            (coord3[0] - coord2[0])*C );

        force2[1] +=  k *( (coord3[1] - 2*coord2[1] + coord1[1])*A -
                            (coord2[1] - coord1[1])*B +
                            (coord3[1] - coord2[1])*C );

        force2[2] +=  k *( (coord3[2] - 2*coord2[2] + coord1[2])*A -
                            (coord2[2] - coord1[2])*B +
                            (coord3[2] - coord2[2])*C );

        //force on i-1, f = k*(A*l1 - B*l2):
        force3[0] +=  k *( (coord2[0] - coord1[0])*A -
                            (coord3[0] - coord2[0])*C );

        force3[1] +=  k *( (coord2[1] - coord1[1])*A -
                            (coord3[1] - coord2[1])*C );

        force3[2] +=  k *( (coord2[2] - coord1[2])*A -
                            (coord3[2] - coord2[2])*C );

        #ifdef CHECKFORCES_INF_NAN
        if(checkNaN_INF<floatingpoint>(force1, 0, 2)||checkNaN_INF<floatingpoint>(force2,0,2)
            ||checkNaN_INF<floatingpoint>(force3,0,2)){
            cout<<"Filament Bending Force becomes infinite. Printing data "<<endl;

            short found = 0;
            Cylinder *cyl1, *cyl2;
            for(auto cyl:Cylinder::getCylinders()){
                auto dbIndex1 = cyl->getFirstBead()->getIndex();
                auto dbIndex2 = cyl->getSecondBead()->getIndex();
                // WARNING this is unsafe because bead starting with index 0 is assumed.
                if(dbIndex1 * 3 == beadSet[n * i] && dbIndex2 * 3 == beadSet[n * i + 1]) {
                    cyl1 = cyl;
                    found++;
                    if(found>=2)
                        break;
                }
                else if(dbIndex1 * 3 == beadSet[n * i + 1] && dbIndex2 * 3 == beadSet[n * i + 2]){
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


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

#ifndef MEDYAN_RestartParams_h
#define MEDYAN_RestartParams_h
#include <vector>
#include "utility.h"
#include "Cylinder.h"

namespace medyan {
using namespace std;

struct restartSystemData{
    floatingpoint time = 0.0;
    int Nfil = -1;
    int Ncyl = -1;
    int Nbead = -1;
    int Nlink = -1;
    int Nmotor = -1;
    int Nbranch = -1;
};
/*Structures to store data parsed from Restartfile*/
struct restartBeadData{
    vector<unsigned int> bsidvec; //stableid
    vector<int> filidvec;
    vector<int> filpos;
    vector<float> coordvec;
    vector<float> forceAuxvec;
};

/*Stores CylData for each cylinder */
struct restartCylData{
    unsigned int cylsid; //stableid
    unsigned int filid;
    unsigned int filtype;
    unsigned int filpos;
    unsigned int beadsidpairvec[2];
    bool endstatusvec[2];
    unsigned int endtypevec[2];
    unsigned int endmonomerpos[2];

    unsigned int totalmonomers;
    floatingpoint eqlen;
    Cylinder* cylinderpointer;
};

/* Stores FilData for each Filament */
struct restartFilData{
    unsigned int filid;
    unsigned int filType;
    vector<unsigned int> cylsidvec;
    Filament* filamentpointer;
};

/*stores MotorData*/
struct restartMotorData{
    unsigned int motorid;
    short motorType;
    unsigned int cylid1;
    unsigned int cylid2;
    short pos1;
    short pos2;
    float eqlen;
    string diffusingspeciesname;
    bool restartcompletion = false;
    int numHeads = 0;
    floatingpoint numBoundHeads = 0;
};

/*stores LinkerData*/
struct restartLinkerData{
    unsigned int linkerid;
    short linkerType;
    unsigned int cylid1;
    unsigned int cylid2;
    short pos1;
    short pos2;
    float eqlen;
    string diffusingspeciesname;
    bool restartcompletion = false;
};
/*stores BrancherData*/
struct restartBrancherData{
    unsigned int branchid;
    short branchType;
    unsigned int cylid1;
    unsigned int cylid2;
    short pos1;
    float eqlen;
    string diffusingspeciesnamebranch;
    string diffusingspeciesnameactin;
    bool restartcompletion = false;
};
/*stores data on Diffusing species in any given compartment*/
struct restartCompartmentDiffusingData{
    int id;
    vector<string> speciesnamevec;
    vector<int> copynumvec;

};

struct restartBulkData{
    vector<string>speciesnamevec;
    vector<int>copynumvec;
};

struct restartcopynumberrallyData{
    vector<string> speciesnamevec;
    vector<int> copynumvec;
};

} // namespace medyan

#endif

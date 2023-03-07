//
// Created by aravind on 9/18/17.
//

#ifndef MEDYAN_CUDAcommon_h
#define MEDYAN_CUDAcommon_h

#include <vector>
#include <list>

#include "FilamentStretchingHarmonic.h"
#include "FilamentBendingHarmonic.h"
#include "FilamentBendingCosine.h"
#include "LinkerStretchingHarmonic.h"
#include "MotorGhostStretchingHarmonic.h"
#include "CylinderExclVolRepulsion.h"
#include "BranchingStretchingHarmonic.h"
#include "BranchingBendingCosine.h"
#include "BranchingDihedralCosine.h"
#include "BranchingPositionCosine.h"
#include "BoundaryCylinderRepulsionExp.h"
#include "CCylinder.h"
#include "common.h"
#include "string.h"
#include "MathFunctions.h"
#include "Structure/Bead.h"
#ifdef SIMDBINDINGSEARCH
#include "Util/DistModule/dist_driver.h"

namespace medyan {
template <uint D>
dist::dOut<D> SIMDoutvar(const uint dim, uint N1, std::initializer_list<float> params) {

    if (dim == 1) {
        dist::dOut<1> out_serialdim(N1, params);
        return out_serialdim;
    }
    else if (dim == 2) {
        dist::dOut<2> out_serialdim(N1, params);
        return out_serialdim;
    }
    else if (dim == 3){
        dist::dOut<3> out_serialdim(N1, params);
        return out_serialdim;
    }
}
#endif

using namespace mathfunc;
struct bin{
    int binID;
    floatingpoint bincoord[3];
    int neighbors[27] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                          -1,-1,-1,-1,-1,-1};
    int binstencilID[27] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                          -1,-1,-1,-1,-1,-1};
};

struct Callbacktime {
	floatingpoint tUpdateBrancherBindingCallback=0.0;
	floatingpoint tUpdateLinkerBindingCallback=0.0;
	floatingpoint tUpdateMotorBindingCallback=0.0;
	floatingpoint tUpdateMotorIDCallback=0.0;
	floatingpoint tFilamentExtensionPlusEndCallback=0.0;
	floatingpoint tFilamentExtensionMinusEndCallback=0.0;
	floatingpoint tFilamentRetractionPlusEndCallback=0.0;
	floatingpoint tFilamentRetractionMinusEndCallback=0.0;
	floatingpoint tFilamentPolymerizationPlusEndCallback=0.0;
	floatingpoint tFilamentPolymerizationMinusEndCallback=0.0;
	floatingpoint tFilamentDepolymerizationPlusEndCallback=0.0;
	floatingpoint tFilamentDepolymerizationMinusEndCallback=0.0;
	floatingpoint tBranchingPointUnbindingCallback=0.0;
	floatingpoint tBranchingCallback=0.0;
	floatingpoint tLinkerUnbindingCallback=0.0;
	floatingpoint tLinkerBindingCallback=0.0;
	floatingpoint tMotorUnbindingCallback=0.0;
	floatingpoint tMotorBindingCallback=0.0;
	floatingpoint tMotorWalkingCallback=0.0;
	floatingpoint tMotorMovingCylinderCallback=0.0;
	floatingpoint tFilamentCreationCallback=0.0;
	floatingpoint tFilamentSeveringCallback=0.0;
	floatingpoint tFilamentDestructionCallback=0.0;
};

struct Callbackcount {
	uint cUpdateBrancherBindingCallback=0;
	uint cUpdateLinkerBindingCallback=0;
	uint cUpdateMotorBindingCallback=0;
	uint cUpdateMotorIDCallback=0;
	uint cFilamentExtensionPlusEndCallback=0;
	uint cFilamentExtensionMinusEndCallback=0;
	uint cFilamentRetractionPlusEndCallback=0;
	uint cFilamentRetractionMinusEndCallback=0;
	uint cFilamentPolymerizationPlusEndCallback=0;
	uint cFilamentPolymerizationMinusEndCallback=0;
	uint cFilamentDepolymerizationPlusEndCallback=0;
	uint cFilamentDepolymerizationMinusEndCallback=0;
	uint cBranchingPointUnbindingCallback=0;
	uint cBranchingCallback=0;
	uint cLinkerUnbindingCallback=0;
	uint cLinkerBindingCallback=0;
	uint cMotorUnbindingCallback=0;
	uint cMotorBindingCallback=0;
	uint cMotorWalkingCallback=0;
	uint cMotorMovingCylinderCallback=0;
	uint cFilamentCreationCallback=0;
	uint cFilamentSeveringCallback=0;
	uint cFilamentDestructionCallback=0;
};


struct PolyPlusEndTemplatetime{
	floatingpoint rxntempate1 = 0.0;
	floatingpoint rxntempate2 = 0.0;
	floatingpoint rxntempate3 = 0.0;
	floatingpoint rxntempate4 = 0.0;

};

struct timeminimization{
	floatingpoint vectorize = 0.0;
	floatingpoint findlambda = 0.0;
	floatingpoint computeforces = 0.0;
	floatingpoint computeenergy = 0.0;
	floatingpoint computeenergyzero = 0.0;
	floatingpoint computeenergynonzero = 0.0;
	vector<floatingpoint> individualenergies;
	vector<floatingpoint> individualenergieszero;
	vector<floatingpoint> individualenergiesnonzero;
	vector<floatingpoint> individualforces;
	floatingpoint endminimization = 0.0;
	floatingpoint tother = 0.0;
	floatingpoint copyforces =0.0;
	floatingpoint stretchingenergy= 0.0;
	floatingpoint stretchingforces= 0.0;
	floatingpoint bendingenergy = 0.0;
	floatingpoint bendingforces = 0.0;
	int computeenergycalls = 0;
	int computeenerycallszero = 0;
	int computeenerycallsnonzero = 0;
	int computeforcescalls = 0;
	int numinteractions[10] = {0,0,0,0,0,0,0,0,0,0};
	floatingpoint timecylinderupdate = 0.0;
	floatingpoint timelinkerupdate = 0.0;
	floatingpoint timemotorupdate = 0.0;
	int callscylinderupdate = 0;
	int callslinkerupdate = 0;
	int callsmotorupdate = 0;
	//
	int motorbindingcalls = 0;
	int motorunbindingcalls = 0;
	int linkerbindingcalls = 0;
	int linkerunbindingcalls = 0;
	int motorwalkingcalls =0;

};

struct motorwalking{
	int stretchtocontract = 0;
	int stretchtostretch = 0;
	int contracttocontract = 0;
	int contracttostretch = 0;
	int equibtocontract = 0;
	int equibtostretch = 0;
	int equibtoequib = 0;
};

struct chemdetails{
    int reactioncount[17]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double dependencytime[17]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int dependentrxncount[17]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double totaltime[17]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double emitsignal[17]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    int diffusion_passivate_count = 0;
    int diffusion_activate_count = 0;
    int ccylclonerxncounter[3] = {0,0,0};
    int ccylclonecounter[3]={0,0,0};
    double ccylclonetimer[3]={0,0,0};
    double internalrxnclone = 0.0;
    double internalrxnadd = 0.0;
    double clonefindspecies= 0.0;
    double getaffectedrxns = 0.0;
};

#if defined(CUDAACCL)
struct CUDAvars {
    floatingpoint * gpu_force = NULL;
    floatingpoint * gpu_forceAux = NULL;
    floatingpoint * gpu_forceAuxP = NULL;
    floatingpoint * gpu_coord =  NULL;
    floatingpoint * gpu_lambda = NULL;
    float vectorize = 0.0;
    floatingpoint * gpu_energy = NULL;
    bool * gpu_btstate = NULL;
    cylinder* gpu_cylindervec = NULL;
#ifdef CUDAACCL
    vector<cudaStream_t*> streamvec;
    vector<cudaEvent_t> eventvec;
#endif
    int* culpritID = NULL;
    int* gculpritID = NULL;
    char* gculpritFF = NULL;
    char* gculpritinteraction = NULL;
    char* culpritFF = NULL;
    char* culpritinteraction = NULL;
    size_t memincuda = 0;
    bool conservestreams = true;
    int offset_E = 0;
    floatingpoint *gpu_energyvec = NULL;
    vector<bool*> backtrackbools;
//    cudaEvent_t *event;

//    float Ccforce = 0.0;
//    float Scforce = 0.0;
//    unsigned int  gpu_sharedMem = 0;
//    unsigned int  gpu_globalMem = 0;
//    int * motorparams;
};

struct CylCylNLvars {
    floatingpoint* gpu_coord;
    floatingpoint* gpu_coord_com;
    int * gpu_beadSet;
    int *gpu_cylID;
    int *gpu_filID;
    int *gpu_filType;
    int *gpu_cmpID;
    int *gpu_fvecpos;
//    int *gpu_cylstate;
    int *gpu_cmon_state_brancher;
    int *gpu_cmon_state_linker;
    int *gpu_cmon_state_motor;
//    int *gpu_cylvecpospercmp;

    bin* bins;

};

struct SERLtime {
    floatingpoint TvectorizeFF = 0.0;
    floatingpoint TcomputeE = 0.0;
    floatingpoint TcomputeEiter = 0.0;
    int Ecount = 0;
    floatingpoint TcomputeF= 0.0;
    floatingpoint Tlambda = 0.0;
    floatingpoint TshiftGrad= 0.0;
    floatingpoint TmaxF= 0.0;
    floatingpoint Tcalculatedot = 0.0;
    vector<floatingpoint>TvecvectorizeFF;
    vector<floatingpoint>TveccomputeE;
    vector<floatingpoint>TveccomputeF;
    vector<floatingpoint>Tlambdavec;
    vector<floatingpoint>Tlambdap;
    vector<floatingpoint>Tlambdapcount;
};

struct CUDAtime {
    floatingpoint TvectorizeFF = 0.0;
    floatingpoint TcomputeE = 0.0;
    floatingpoint TcomputeEiter = 0.0;
    int Ecount = 0;
    floatingpoint TcomputeF= 0.0;
    floatingpoint Tlambda = 0.0;
    floatingpoint TshiftGrad= 0.0;
    floatingpoint TmaxF= 0.0;
    floatingpoint Tcalculatedot = 0.0;
    floatingpoint Tstartmin = 0.0;
    vector<floatingpoint>TvecvectorizeFF;
    vector<floatingpoint>TveccomputeE;
    vector<floatingpoint>TveccomputeF;
    vector<floatingpoint>Tlambdavec;
    vector<floatingpoint>Tlambdap;
    vector<floatingpoint>Tlambdapcount;
};
#endif
class CUDAcommon{
public:
    static Callbacktime ctime;
	static Callbackcount ccount;
	static PolyPlusEndTemplatetime ppendtime;
	static timeminimization tmin;
	static motorwalking mwalk;
	static chemdetails cdetails;

#ifdef CUDAACCL
    static CylCylNLvars cylcylnlvars;
    static SERLtime serltime;
    static CUDAtime cudatime;
    static const CUDAvars& getCUDAvars(){return cudavars;}
    static const CylCylNLvars& getCylCylNLvars(){return cylcylnlvars;}
    static void handleerror(cudaError_t a){
        if(a !=cudaSuccess){
            cout<<cudaGetErrorString(a)<<endl;
            exit(EXIT_FAILURE);
        }
    }

    static void handleerror(cudaError_t a, std::string tag1, std::string tag2){
//        cout<<CUDAcommon::getCUDAvars().culpritFF[0]<<CUDAcommon::getCUDAvars().culpritFF[1]<<endl;
        if(a == cudaErrorAssert){
            cudaDeviceSynchronize();
            auto culpritFF = CUDAcommon::getCUDAvars().culpritFF;
            auto culpritinteraction = CUDAcommon::getCUDAvars().culpritinteraction;

            if(strcmp(culpritinteraction, "Filament Stretching Harmonic")==0)
                FilamentStretchingHarmonic::checkforculprit();
            else if(strcmp(culpritinteraction, "Filament Bending Harmonic")==0)
                FilamentBendingHarmonic::checkforculprit();
            else if(strcmp(culpritinteraction, "Filament Bending Cosine")==0)
                FilamentBendingCosine::checkforculprit();
            else if(strcmp(culpritinteraction, "Linker Stretching Harmonic")==0)
                LinkerStretchingHarmonic::checkforculprit();
            else if(strcmp(culpritinteraction, "Motor Stretching Harmonic")==0)
                MotorGhostStretchingHarmonic::checkforculprit();
            else if(strcmp(culpritinteraction, "Cylinder Excluded Volume")==0)
                CylinderExclVolRepulsion::checkforculprit();
            else if(strcmp(culpritinteraction, "Branching Stretching Harmonic")==0)
                BranchingStretchingHarmonic::checkforculprit();
            else if(strcmp(culpritinteraction, "Branching Bending Cosine")==0)
                BranchingBendingCosine::checkforculprit();
            else if(strcmp(culpritinteraction, "Branching Dihedral Cosine")==0)
                BranchingDihedralCosine::checkforculprit();
            else if(strcmp(culpritinteraction, "Branching Position Cosine")==0)
                BranchingPositionCosine::checkforculprit();
            else if(strcmp(culpritinteraction, "Boundary Cylinder Repulsion Exp")==0)
                BoundaryCylinderRepulsionExp::checkforculprit();
            else{
                cout<<"unknown assert error. Check code. Exiting.."<<endl;
                exit(EXIT_FAILURE);
            }
        }
        else if(a !=cudaSuccess && a != cudaErrorAssert){
            cout<<cudaGetErrorString(a)<<endl;
            cout<<"ERROR! "<<tag1<<"error in "<<tag2<<". Check vectors. Exiting.."<<endl;
            exit(EXIT_FAILURE);
        }
    }

    static void printculprit(std::string tag1, std::string tag2){
        cout << "Energy of system became infinite. Try adjusting minimization parameters." << endl;
        cout << "The culprit was ... " << tag1 << endl;
        cout << "Culprit interaction = " << tag2 << endl;
    }
#endif
};

} // namespace medyan

#endif
//CUDA_VEC_CUDACOMMON_H

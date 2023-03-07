
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------
#ifndef MEDYAN_HybridBindingSearchManager_h
#define MEDYAN_HybridBindingSearchManager_h
#ifdef SIMDBINDINGSEARCH
#include "Util/DistModule/dist_driver.h"
#endif
#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
#include <unordered_map>
#include <unordered_set>
#include <random>
#include <algorithm>

#include "common.h"

#include "NeighborListImpl.h"
#include "HybridNeighborListImpl.h"
#include "ReactionBase.h"

#include "SysParams.h"
#include "Rand.h"
#include "BindingManager.h"

namespace medyan {
//FORWARD DECLARATIONS
class SubSystem;
class ReactionBase;
class CCylinder;
class Compartment;
class Cylinder;
class FilamentBindingManager;


class HybridBindingSearchManager {

    friend class ChemManager;

private:
    static const bool CROSSCHECK_BS_SWITCH = false;
    static const uint switchfactor  = 10;

    chrono::high_resolution_clock::time_point minsSIMD, mineSIMD, minsHYBD, mineHYBD,
            minsfind, minefind, minsmap, minemap;
    Compartment* _compartment;

    vector<vector<floatingpoint>> bindingsites1;
    vector<vector<floatingpoint>> bindingsites2;
    vector<vector<float>> _rMaxsqvec; //squared maxdistance cutoff
    vector<vector<float>> _rMinsqvec; //squared mindistance cutoff
    vector<vector<short>> _filamentIDvec;//filament ID pairs considered
    static vector<short> HNLIDvec; //Hybrid NL ID to track the total number of
    // neighborilists in the system

    vector<vector<short>> bstateposvec;
    vector<vector<float>> minparamcyl2;
    vector<vector<float>> maxparamcyl2;
    vector<vector<int>> Nbindingpairs;
    int totaluniquefIDpairs = 0;

    vector<vector<FilamentBindingManager*>> fManagervec;

    //possible bindings at current state. updated according to neighbor list
    //possible bindings at current state. updated according to neighbor list stencil

//    vector<vector<unordered_multimap<tuple<CCylinder*, short>, tuple<CCylinder*, short>>>>
//            _possibleBindingsstencilvec;

    /*vector<vector<unordered_multimap<uint32_t, uint32_t>>>
            _mpossibleBindingsstencilvecuint;*/

    vector<vector<unordered_map<uint32_t, vector<uint32_t>>>>
            _possibleBindingsstencilvecuint;

//    vector<vector<unordered_map<tuple<CCylinder*, short>, vector<tuple<CCylinder*,
//    short>>>>>_reversepossibleBindingsstencilvec;

    vector<vector<unordered_map<uint32_t, vector<uint32_t>>>>
    _reversepossibleBindingsstencilvecuint;

    vector<uint32_t> linker1, linker2;
    vector<uint32_t> motor1, motor2;
    uint Npairs = 0;

    //static neighbor list
    static HybridCylinderCylinderNL* _HneighborList;

    //SIMD variables
    unsigned mask = 0;
#ifdef SIMDBINDINGSEARCH
    static constexpr dist::tag_simd<dist::simd_avx_par,  float>  t_avx_par {};
    static constexpr dist::tag_simd<dist::simd_avx,  float>   t_avx {};
    static constexpr dist::tag_simd<dist::simd_no,   float>   t_serial {};
    bool initialized = false;

    //listing 12 variables, to support upto 8 distance pairs calculations.
    inline static dist::dOut<1U,true> bspairsself[8];
    inline static dist::dOut<1U,false> bspairs[8];


	template<uint D, bool SELF>
	dist::dOut<D,SELF>& getdOut(short dOutID);

//    static dist::dOut<1U,false> bspairslinker;
    static dist::dOut<1U,true> bspairslinkerself;
//    static dist::dOut<1U,false> bspairsmotor;
    static dist::dOut<1U,true> bspairsmotorself;
    static dist::dOut<1U,false> bspairslinker2;
    static dist::dOut<1U,false> bspairsmotor2;
#endif
    vector<uint32_t> pairslinker;
    vector<uint32_t> pairsmotor;
    vector<bool> filID_fpos_pairL;
    vector<bool> filID_fpos_pairM;
    vector<vector<bool>> pairvaluespecieslinker;
    vector<vector<bool>> pairvaluespeciesmotor;

    /*Partitioned volume refers to partitioning compartment volume in to smaller sub
volumes namely self(1), halves(6), quarters(12) and 1/8ths(8). The position in the
 vector corresponds to stencil ID. the data corresponds to the position in the
 partioned coordinate vector of vectors. 27 is a dummy entry and should never be
 called. For example, stencil 4 means that the compartment gives parition 5
     (partitioned_volume_ID[4] = 5) and so on.*/
    /*Note. There is hard coded  referencing of complimentary sub volumes. For example,
 * the top plane of a compartment has paritioned_volume_ID 1 and the bottom plane
 * (complimentary) is 2. So, if permuting through C1 and C2, and if C2 is on top of C1,
 * C2's stencil ID in C1 is 14 (Hard Coded in GController.cpp). C1 will compare
 * paritioned_volume_ID 1 with C2's paritioned_volume_ID 2*/
    uint partitioned_volume_ID[27] = {27, 27, 27, 27, 5, 7, 19, 9, 21, 27, 27, 27, 27, 0, 1,
                                      11, 3, 13, 27, 27, 27, 27, 27, 15, 23, 17, 25};

#ifdef SIMDBINDINGSEARCH

    template <uint D, bool SELF, bool LinkerorMotor>
    void calculatebspairsLMselfV3(dist::dOut<D, SELF>& bspairs, short idvec[2]);

    template <uint D, bool SELF, bool LinkerorMotor>
    void calculatebspairsLMenclosedV3(dist::dOut<D, SELF>& bspairs, dist::dOut<D, SELF>&
            bspairs2, short idvec[2]);

    template <uint D, bool SELF, bool LinkerorMotor>
    void checkcontacts(dist::dOut<D, SELF>& bspairs, dist::dOut<D, SELF>& bspairs2);

    static const short nthreads = 1;

    vector<vector<uint32_t>> valuematrixvec[2*nthreads];

    template <uint D, bool SELF, bool LinkerorMotor>
    void gatherCylindercIndexV3(dist::dOut<D,SELF>& bspairsoutS, int first, int
    last, short idvec[2], Compartment* nCmp = NULL);


    static short Totallinkermotor;

#endif
   void countNpairsfound(short idvec[2]);

public:

    //constructors
     HybridBindingSearchManager(Compartment* compartment);
    ~HybridBindingSearchManager() {}

    //@{
    //@{
    /// Getters for distances
//    float getRMin() {return _rMin;}
//    float getRMax() {return _rMax;}
    //@}

    bool isConsistent();

    void addPossibleBindingsstencil(short idvec[2], CCylinder* cc, short bindingSite);
    void removePossibleBindingsstencil(short idvec[2], CCylinder* cc, short bindingSite);

    void appendPossibleBindingsstencil(short idvec[2],
                                  CCylinder* ccyl1,
                                  CCylinder* ccyl2,
                                  short site1,
                                  short site2);

    ///update all possible binding reactions that could occur using stencil NL
    void updateAllPossibleBindingsstencilSIMDV3();
    void updateAllPossibleBindingsstencilHYBD();
	void updateAllBindingReactions();

    vector<tuple<CCylinder*, short>> chooseBindingSitesstencil(short idvec[2]);

    void setbindingsearchparameter(FilamentBindingManager* fmanager, short bstatepos,
                                   short ftype1 = 0, short ftype2 = 0, float rMax = 0.0,
                                   float rMin = 0.0);

    void addtoHNeighborList();
    void copyInfotoBindingManagers();

    vector<Cylinder*> getHNeighbors(Cylinder* c, short HNLID){
        return _HneighborList->getNeighborsstencil(HNLID, c);
    }

    void checkoccupancy(short idvec[2]);

    void checkoccupancySIMD(short idvec[2]);

    void printbindingsizes(){
        int idx, idx2;
        for(idx = 0; idx<totaluniquefIDpairs; idx++){
            int countbounds = _rMaxsqvec[idx].size();
            for (idx2 = 0; idx2 < countbounds; idx2++) {
            	cout<<Nbindingpairs[idx][idx2]<<" ";
//                std::cout<<"Hybrid "<<_possibleBindingsstencilvecuint[idx][idx2].size()<<
//                " SIMD "<<Nbindingpairs[idx][idx2]<<endl;
            }
            cout<<endl;
        }
    }

    void printbindingsitesstencil(short idvec[2]);

    void resetpossibleBindings(){

        int idx, idx2;
        for(idx = 0; idx<totaluniquefIDpairs; idx++){
            int countbounds = _rMaxsqvec[idx].size();
            for (idx2 = 0; idx2 < countbounds; idx2++) {
                _possibleBindingsstencilvecuint[idx][idx2].clear();
//                _mpossibleBindingsstencilvecuint[idx][idx2].clear();
                _reversepossibleBindingsstencilvecuint[idx][idx2].clear();
            }
        }

    }

    int numBindingSitesstencil(short idvec[2]) {

        return Nbindingpairs[idvec[0]][idvec[1]];
    }

    /*static void setdOut(){
        Totallinkermotor = 2;
*//*        bspairs2self.init_dout(10000, {900.0f, 1600.0f, 30625.0f, 50625.0f});
        bspairs2.init_dout(10000, {900.0f, 1600.0f, 30625.0f, 50625.0f});*//*

*//*        bspairs2_D1.init_dout(10000,{900.0f,1600.0f});
        bspairs2_D2.init_dout(10000,{30625.0f, 50625.0f});*//*

        // V2
        bspairslinkerself.init_dout(10000,{900.0f,1600.0f});
        bspairslinker.init_dout(10000,{900.0f,1600.0f});
        bspairsmotorself.init_dout(10000,{30625.0f, 50625.0f});
        bspairsmotor.init_dout(10000,{30625.0f, 50625.0f});

        bspairslinker2.init_dout(10000,{900.0f,1600.0f});
        bspairsmotor2.init_dout(10000,{30625.0f, 50625.0f});
    }*/

    void clearPossibleBindingsstencil(short idvec[2]);
    void initializeSIMDvars();

    static floatingpoint largestlinkerdistance;
    static floatingpoint largestmotordistance;

    static floatingpoint SIMDtime;
    static floatingpoint HYBDtime;
    static floatingpoint findtime;
    static floatingpoint findtimeV2;
    static floatingpoint appendtime;
    static floatingpoint SIMDparse1;
    static floatingpoint SIMDparse2;
    static floatingpoint SIMDparse3;
    static floatingpoint SIMDcountbs;
    static floatingpoint HYBDappendtime;
    static floatingpoint SIMDV3appendtime;
    static floatingpoint findtimeV3;

};

#ifdef SIMDBINDINGSEARCH
template<>
dist::dOut<1,true>& HybridBindingSearchManager::getdOut(short dOutID);
template<>
dist::dOut<1,false>& HybridBindingSearchManager::getdOut(short dOutID);
#endif
#endif

} // namespace medyan

#endif


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
#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
#include "Compartment.h"
#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"
#include "HybridBindingSearchManager.h"
#include "HybridNeighborListImpl.h"
#include "MotorGhost.h"
#include "MathFunctions.h"
#include "Controller/GController.h"
#include "SysParams.h"
#include "CUDAcommon.h"
#include "Rand.h"

#include "Controller/CController.h"

namespace medyan {

vector<short> HybridBindingSearchManager::HNLIDvec;
using namespace mathfunc;

void HybridBindingSearchManager::setbindingsearchparameter(
	FilamentBindingManager* fmanager,
	short bstatepos, short ftype1, short ftype2, float rMax, float rMin) {

    unordered_map<uint32_t, vector<uint32_t>> tempuint;
    unordered_map<uint32_t, vector<uint32_t>> rtempuint;
    int tempNbind;

    bool isfound = false;
    vector<short> ftypepairs;
    if(ftype1 < ftype2)
        ftypepairs = {ftype1,ftype2};
    else
        ftypepairs = {ftype2, ftype1};
    for(short idx = 0; idx < totaluniquefIDpairs; idx++){
        vector<short> fIDpair = _filamentIDvec[idx];
        if(fIDpair[0] != ftypepairs[0] || fIDpair[1] != ftypepairs[1]) continue;
        else {
            isfound = true;
            _rMaxsqvec[idx].push_back(rMax *rMax);
            _rMinsqvec[idx].push_back(rMin *rMin);
            fManagervec[idx].push_back(fmanager);

            _possibleBindingsstencilvecuint[idx].push_back(tempuint);
            _reversepossibleBindingsstencilvecuint[idx].push_back(rtempuint);
            bstateposvec[idx].push_back(bstatepos);
            Nbindingpairs[idx].push_back(tempNbind);
            break;
        }
    }
    if(isfound == false){
        vector<unordered_map<uint32_t, vector<uint32_t>>> temp2uint;
        vector<unordered_map<uint32_t, vector<uint32_t>>> rtemp2uint;

        vector<int> tempNbind2={0};
        temp2uint.push_back(tempuint);
        rtemp2uint.push_back(rtempuint);

        Nbindingpairs.push_back(tempNbind2);
        vector<float> localrmaxsq ={rMax * rMax};
        vector<float> localrminsq = {rMin * rMin};
        vector<FilamentBindingManager*> localfmanager;
        vector<short> localbstateposvec = {bstatepos};
        _rMaxsqvec.push_back(localrmaxsq);
        _rMinsqvec.push_back(localrminsq);
        localfmanager.push_back(fmanager);
        fManagervec.push_back(localfmanager);
        _filamentIDvec.push_back(ftypepairs);
        _possibleBindingsstencilvecuint.push_back(temp2uint);
        _reversepossibleBindingsstencilvecuint.push_back(rtemp2uint);
        bstateposvec.push_back(localbstateposvec);
        vector<floatingpoint> bs1, bs2;
        vector<float> minvec = {(float)*(SysParams::Chemistry().bindingSites[ftypepairs[0]]
                .begin())/ SysParams::Geometry().cylinderNumMon[ftypepairs[0]],
                                (float)*(SysParams::Chemistry().bindingSites[ftypepairs[1]]
                                        .begin())/ SysParams::Geometry().cylinderNumMon[ftypepairs[1]]};
        vector<float> maxvec = {(float)*(SysParams::Chemistry().bindingSites[ftypepairs[0]]
                                                    .end() -1)/ SysParams::Geometry().cylinderNumMon[ftypepairs[0]],
                                (float)*(SysParams::Chemistry().bindingSites[ftypepairs[1]]
                                                    .end() -1)/ SysParams::Geometry().cylinderNumMon[ftypepairs[1]]};
        for(auto it1 = SysParams::Chemistry().bindingSites[ftypepairs[0]].begin();
            it1 != SysParams::Chemistry().bindingSites[ftypepairs[0]].end(); it1++) {
            bs1.push_back((float)*it1 / SysParams::Geometry().cylinderNumMon[ftypepairs[0]]);
        }
        for(auto it1 = SysParams::Chemistry().bindingSites[ftypepairs[1]].begin();
            it1 != SysParams::Chemistry().bindingSites[ftypepairs[1]].end(); it1++) {
            bs2.push_back((float)*it1 / SysParams::Geometry().cylinderNumMon[ftypepairs[1]]);
        }
        minparamcyl2.push_back(minvec);
        maxparamcyl2.push_back(maxvec);
        bindingsites1.push_back(bs1);
        bindingsites2.push_back(bs1);
        totaluniquefIDpairs++;
    }
}

HybridBindingSearchManager::HybridBindingSearchManager(Compartment* compartment){
	mask = (1 << SysParams::Chemistry().shiftbybits) - 1;

    _compartment = compartment;
    totaluniquefIDpairs = 0;
}

#ifdef SIMDBINDINGSEARCH
void HybridBindingSearchManager::initializeSIMDvars(){
	short count = 0;
	for (short idx = 0; idx < totaluniquefIDpairs; idx++) {
		int countbounds = _rMaxsqvec[idx].size();
		for (short idx2 = 0; idx2 < countbounds; idx2++) {
			if(bstateposvec[idx][idx2] == 1) //Linker
				largestlinkerdistance = max<floatingpoint>(largestlinkerdistance,
						sqrt(_rMaxsqvec[idx][idx2]));
			else//Motor
				largestmotordistance = max<floatingpoint>(largestmotordistance,
				                                           sqrt(_rMaxsqvec[idx][idx2]));

			auto coord = _compartment->coordinates();

			getdOut<1U, true>(count).init_dout(10000, {_rMinsqvec[idx][idx2],
									 _rMaxsqvec[idx][idx2]});
			getdOut<1U, false>(count).init_dout(10000, {_rMinsqvec[idx][idx2],
			                                           _rMaxsqvec[idx][idx2]});
			count++;
		}
	}
}

template<>
dist::dOut<1,true>& HybridBindingSearchManager::getdOut(short dOutID){
	if(dOutID < 8)
		return bspairsself[dOutID];
	else{
		cout<<"Illegal ID is requested from getdOut function. Code accomodates "
				  "upto 6 (linker + motor types) Exiting."<<endl;
					exit(EXIT_FAILURE);
	}
}

template<>
dist::dOut<1,false>& HybridBindingSearchManager::getdOut(short dOutID){
	if(dOutID < 8)
		return bspairs[dOutID];
	else{
		cout<<"Illegal ID is requested from getdOut function. Code accomodates "
		      "upto 6 (linker + motor types) Exiting."<<endl;
		exit(EXIT_FAILURE);
	}
}

template <uint D, bool SELF, bool LinkerorMotor>
void HybridBindingSearchManager::calculatebspairsLMselfV3(dist::dOut<D, SELF>&
bspairsoutSself, short idvec[2]) {

	auto filTypepairs =_filamentIDvec[idvec[0]];
	auto boundstate = SysParams::Mechanics().speciesboundvec;

	minsfind = chrono::high_resolution_clock::now();
	int C1size = _compartment->getSIMDcoordsV3<LinkerorMotor>(0, filTypepairs[0]).size();

	if (C1size > 0) {
		bspairsoutSself.reset_counters();
		if(CROSSCHECK_BS_SWITCH){
            if (C1size >= switchfactor * dist::get_simd_size(t_avx))
                HybridNeighborList::_crosscheckdumpFileNL << "SELF t_avx" << endl;
            else
                HybridNeighborList::_crosscheckdumpFileNL << "SELF t_serial" << endl;
        }

		if(filTypepairs[0] == filTypepairs[1]) {
			if (C1size >= switchfactor * dist::get_simd_size(t_avx))
				dist::find_distances(bspairsoutSself,
				                     _compartment->getSIMDcoordsV3<LinkerorMotor>(0, filTypepairs[0]),
				                     t_avx);
			else
				dist::find_distances(bspairsoutSself,
				                     _compartment->getSIMDcoordsV3<LinkerorMotor>(0, filTypepairs[0]), t_serial);
		}
		else{
			int C2size = _compartment->getSIMDcoordsV3<LinkerorMotor>(0, filTypepairs[1]).size();
			if (C1size >= switchfactor * dist::get_simd_size(t_avx) &&
			    C2size >= switchfactor * dist::get_simd_size(t_avx)) {

				dist::find_distances(bspairsoutSself,
				                     _compartment->getSIMDcoordsV3<LinkerorMotor>(0, filTypepairs[0]),
				                     _compartment->getSIMDcoordsV3<LinkerorMotor>(0, filTypepairs[1]),
				                     t_avx);

			} else {

				dist::find_distances(bspairsoutSself,
				                     _compartment->getSIMDcoordsV3<LinkerorMotor>(0, filTypepairs[0]),
				                     _compartment->getSIMDcoordsV3<LinkerorMotor>(0, filTypepairs[1]), t_serial);

			}

		}

		minefind = chrono::high_resolution_clock::now();
		chrono::duration<floatingpoint> elapsed_runfind(minefind - minsfind);
		findtimeV3 += elapsed_runfind.count();


		//MERGE INTO single vector
		//@{
		if(CROSSCHECK_BS_SWITCH)
		    HybridNeighborList::_crosscheckdumpFileNL << "SELF gather" << endl;
		if (true) {
			uint N = bspairsoutSself.counter[D - 1];
//		cout<<" Contacts found by SIMD Self "<<N<<endl;
			minsfind = chrono::high_resolution_clock::now();
			if (nthreads == 1)
				gatherCylindercIndexV3<D, SELF, LinkerorMotor>
						(bspairsoutSself, 0, N, idvec, _compartment);

			minefind = chrono::high_resolution_clock::now();
			chrono::duration<floatingpoint> elapsed_append(minefind - minsfind);
			SIMDV3appendtime += elapsed_append.count();
		}
	}
	//@}
}

template<uint D, bool SELF, bool LinkerorMotor>
void HybridBindingSearchManager::calculatebspairsLMenclosedV3(dist::dOut<D,SELF>&
bspairsoutS, dist::dOut<D,SELF>& bspairsoutS2, short idvec[2]){


	auto filTypepairs =_filamentIDvec[idvec[0]];
	auto boundstate = SysParams::Mechanics().speciesboundvec;
	auto upnstencilvec = _compartment->getuniquepermuteneighborsstencil();
	short i =0;

	for(auto ncmp: _compartment->getuniquepermuteNeighbours()) {

		short pos = upnstencilvec[i];
		/*Note. There is hard coded  referencing of complimentary sub volumes. For example,
		 * the top plane of a compartment has paritioned_volume_ID 1 and the bottom plane
		 * (complimentary) is 2. So, if permuting through C1 and C2, and if C2 is on top of C1,
		 * C2's stencil ID in C1 is 14 (Hard Coded in GController.cpp). C1 will compare
		 * paritioned_volume_ID 1 with C2's paritioned_volume_ID 2*/
		int C1size = _compartment->getSIMDcoordsV3<LinkerorMotor>
				(partitioned_volume_ID[pos], filTypepairs[0]).size();
		int C2size = ncmp->getSIMDcoordsV3<LinkerorMotor>
				(partitioned_volume_ID[pos] + 1, filTypepairs[1]).size();

		if (CROSSCHECK_BS_SWITCH && C1size > 0 && C2size > 0) {
			if (C1size >= switchfactor * dist::get_simd_size(t_avx) &&
			    C2size >= switchfactor * dist::get_simd_size(t_avx))
				HybridNeighborList::_crosscheckdumpFileNL << "ENCLOSED t_avx" << endl;
			else
				HybridNeighborList::_crosscheckdumpFileNL << "ENCLOSED t_serial" << endl;
		}

		if (C1size > 0 && C2size > 0) {

			minsfind = chrono::high_resolution_clock::now();

			bspairsoutS.reset_counters();
			if (C1size >= switchfactor * dist::get_simd_size(t_avx) &&
			    C2size >= switchfactor * dist::get_simd_size(t_avx)) {

				dist::find_distances(bspairsoutS,
				                     _compartment->getSIMDcoordsV3<LinkerorMotor>
						                     (partitioned_volume_ID[pos], filTypepairs[0]),
				                     ncmp->getSIMDcoordsV3<LinkerorMotor>
						                     (partitioned_volume_ID[pos] + 1, filTypepairs[1]), t_avx);

			} else {

				dist::find_distances(bspairsoutS,
				                     _compartment->getSIMDcoordsV3<LinkerorMotor>
						                     (partitioned_volume_ID[pos], filTypepairs[0]),
				                     ncmp->getSIMDcoordsV3<LinkerorMotor>
						                     (partitioned_volume_ID[pos] + 1, filTypepairs[1]),
				                     t_serial);

			}

			minefind = chrono::high_resolution_clock::now();
			chrono::duration<floatingpoint> elapsed_runfind(minefind - minsfind);
			findtimeV3 += elapsed_runfind.count();

			i++;

			//MERGE INTO single vector
			//@{
            if(CROSSCHECK_BS_SWITCH)
                HybridNeighborList::_crosscheckdumpFileNL << "ENCLOSED gather" << endl;
			if (true) {
				minsfind = chrono::high_resolution_clock::now();
				uint N = bspairsoutS.counter[D-1];
				if (nthreads == 1)
					gatherCylindercIndexV3<D, SELF, LinkerorMotor>
							(bspairsoutS, 0, N, idvec, ncmp);
				else {
					LOG(ERROR)<<"Multiple thread support not available in MEDYAN. Exiting"
				 "."<<endl;

					std::vector<std::thread> threads_avx;
					uint nt = nthreads;
					threads_avx.reserve(nt);
					uint prev = 0;
					uint frac = N / nt;
					uint next = frac + N % nt;
					for (uint i = 0; i < nt; ++i) {
						threads_avx.push_back(std::thread
								                      (&HybridBindingSearchManager::gatherCylindercIndexV3<D, SELF, LinkerorMotor>,
								                       this, std::ref(bspairsoutS), prev, next, idvec, ncmp));
						prev = next;
						next = min(N, prev + frac);
					}

					//Join
					for (auto &t : threads_avx)
						t.join();
					threads_avx.clear();
				}//else
				minefind = chrono::high_resolution_clock::now();
				chrono::duration<floatingpoint> elapsed_append(minefind - minsfind);
				SIMDV3appendtime += elapsed_append.count();
			}

		}
		//@}
	}
}

template <uint D, bool SELF, bool LinkerorMotor>
void HybridBindingSearchManager::gatherCylindercIndexV3(dist::dOut<D,SELF>&
bspairsoutS, int first, int last, short idvec[2], Compartment* nCmp){

	const auto& cylinderInfoData = Cylinder::getDbData();
	unsigned int count64 = 0;//counter to bits in random integer
	bitset<64> randInt = 0;
	short idx = idvec[0];
	short idx2 = idvec[1];
	auto &nCmppbs = nCmp->getHybridBindingSearchManager()
			->_possibleBindingsstencilvecuint;
	auto &nCmprpbs = nCmp->getHybridBindingSearchManager()
			->_reversepossibleBindingsstencilvecuint;

	for(uint pid = first; pid < last; pid++) {
		uint32_t t1 = bspairsoutS.dout[2 * (D - 1)][pid];
		uint32_t t2 = bspairsoutS.dout[2 * (D - 1) + 1][pid];

		if(SELF == true){
			_possibleBindingsstencilvecuint[idx][idx2][t1].push_back(t2);

			_reversepossibleBindingsstencilvecuint[idx][idx2][t2].push_back(t1);
		}
		else {
			//Generate random number of 64 bits
			if(count64 == 0) {
				randInt = bitset<64>(Rand::randUInt64bit());
//				cout<<randInt<<endl;
			}
			if(randInt[count64]){
				_possibleBindingsstencilvecuint[idx][idx2][t1].push_back(t2);

				_reversepossibleBindingsstencilvecuint[idx][idx2][t2].push_back(t1);
			}
			else{
				nCmppbs[idx][idx2][t2].push_back(t1);
				nCmprpbs[idx][idx2][t1].push_back(t2);
			}
			//Keep count
			count64++;
//			cout<<count64<<" "<<first<<" "<<last<<endl;
			//reset if you reach 64 bits
			if(count64>63)
				count64=0;

			/*if(t1>t2) {
				_possibleBindingsstencilvecuint[idvec[0]][idvec[1]][t1].push_back(t2);

				_reversepossibleBindingsstencilvecuint[idvec[0]][idvec[1]][t2].push_back(t1);
			}
			else{
				nCmp->getHybridBindingSearchManager()
				->_possibleBindingsstencilvecuint[idvec[0]][idvec[1]][t2].push_back(t1);
				nCmp->getHybridBindingSearchManager()
				->_reversepossibleBindingsstencilvecuint[idvec[0]][idvec[1]][t1].push_back(t2);
			}*/
		}
	}
}
#endif

void HybridBindingSearchManager::addPossibleBindingsstencil(short idvec[2],
                                 CCylinder* cc, short bindingSite) {
	if (SysParams::INITIALIZEDSTATUS ) {
		short idx = idvec[0];
		short idx2 = idvec[1];

		auto fIDpair = _filamentIDvec[idx].data();
		short bstatepos = bstateposvec[idx][idx2];
//		cout<<"bstatepos "<<bstatepos<<endl;
		short _filamentType = cc->getType();
		short HNLID = HNLIDvec[idx];
		short complimentaryfID;

		if (_filamentType != fIDpair[0] && _filamentType != fIDpair[1]) return;
		else if (_filamentType == fIDpair[0]) complimentaryfID = fIDpair[1];
		else complimentaryfID = fIDpair[0];

		bool bstatecheck = false;
		float _rMaxsq = _rMaxsqvec[idx][idx2];
		float _rMinsq = _rMinsqvec[idx][idx2];

		int bindingsitestep = 1;

		if (bstatepos == 1) {
			bindingsitestep = SysParams::Chemistry().linkerbindingskip-1;
			bstatecheck = areEqual(cc->getCMonomer(bindingSite)->speciesBound(
					SysParams::Chemistry().linkerBoundIndex[_filamentType])->getN(), 1.0);
		}

		if (bstatepos == 2) {
			bstatecheck = areEqual(cc->getCMonomer(bindingSite)->speciesBound(
					SysParams::Chemistry().motorBoundIndex[_filamentType])->getN(), 1.0);
		}

		/*if Minus End, check if binding site exists*/
		auto cylinder = cc->getCylinder();
		if(cylinder->isMinusEnd()){
			auto sf = cc->getCMonomer(bindingSite)->activeSpeciesFilament();
			if(sf == -1)
				return;
		}

		//add valid binding sites
		if (bstatecheck) {

			uint32_t shiftedIndex1 = cc->getCylinder()->getStableIndex();
			shiftedIndex1 = shiftedIndex1 << SysParams::Chemistry().shiftbybits;

			uint32_t pos = find(SysParams::Chemistry().bindingSites[_filamentType].begin(),
			                    SysParams::Chemistry().bindingSites[_filamentType].end(),
			                    bindingSite)
			               - SysParams::Chemistry().bindingSites[_filamentType].begin();

			// Do not add if the linker to be added does not satisfy this condition.
			if(pos % bindingsitestep != 0) return;

			uint32_t t1 = shiftedIndex1|pos;
			//loop through neighbors
			//now re add valid based on CCNL
			vector<Cylinder *> Neighbors;

			Neighbors = _HneighborList->getNeighborsstencil(HNLID, cc->getCylinder());

			for (auto cn : Neighbors) {
				Cylinder *c = cc->getCylinder();
				short _nfilamentType = cn->getType();

				if (_nfilamentType != complimentaryfID) return;
				if (cn->getParent() == c->getParent()) continue;

				auto ccn = cn->getCCylinder();
				int k = 0;

				for (int itI = 0; itI < SysParams::Chemistry().bindingSites[_nfilamentType].size(); itI += bindingsitestep) {

					auto it = SysParams::Chemistry().bindingSites[_nfilamentType].begin() + itI;

					bool filstatecheckn = true;
					if(cn->isMinusEnd()) {
						auto sfn = ccn->getCMonomer(*it)->activeSpeciesFilament();
						if (sfn == -1)
							filstatecheckn = false;
					}

					bool bstatecheckn = false;

					if (bstatepos == 1)
						bstatecheckn = areEqual(ccn->getCMonomer(*it)->speciesBound(
								SysParams::Chemistry().linkerBoundIndex[_nfilamentType])->getN(),
						                        1.0);
					if (bstatepos == 2)
						bstatecheckn = areEqual(ccn->getCMonomer(*it)->speciesBound(
								SysParams::Chemistry().motorBoundIndex[_nfilamentType])->getN(),
						                        1.0);

					if (bstatecheckn && filstatecheckn) {
//						cout<<"Adding neighbor "<<cn->getID()<<" "<<*it<<" "<<k<<endl;
						//check distances..
						auto mp1 = (float) bindingSite /
						           SysParams::Geometry().cylinderNumMon[_filamentType];
						auto mp2 = (float) *it /
						           SysParams::Geometry().cylinderNumMon[_nfilamentType];

						auto x1 = c->getFirstBead()->vcoordinate();
						auto x2 = c->getSecondBead()->vcoordinate();
						auto x3 = cn->getFirstBead()->vcoordinate();
						auto x4 = cn->getSecondBead()->vcoordinate();

						auto m1 = midPointCoordinate(x1, x2, mp1);
						auto m2 = midPointCoordinate(x3, x4, mp2);

						floatingpoint distsq = twoPointDistancesquared(m1, m2);

						if (distsq > _rMaxsq || distsq < _rMinsq) {k = k + bindingsitestep;continue;}

						uint32_t shiftedIndex2 = cn->getStableIndex() << SysParams::Chemistry().shiftbybits;

						uint32_t t2 = shiftedIndex2|k;

						//add in correct order
//						_mpossibleBindingsstencilvecuint[idx][idx2].emplace(t1,t2);
						_possibleBindingsstencilvecuint[idx][idx2][t1].push_back(t2);
						_reversepossibleBindingsstencilvecuint[idx][idx2][t2].push_back(t1);
					}
					k = k + bindingsitestep;
				}
			}

			//Calculate N
			countNpairsfound(idvec);
			fManagervec[idx][idx2]->updateBindingReaction(Nbindingpairs[idx][idx2]);

			minefind = chrono::high_resolution_clock::now();
			chrono::duration<floatingpoint> elapsed_countsites(minefind - minsfind);
			SIMDcountbs += elapsed_countsites.count();
		}
//		checkoccupancySIMD(idvec);
	}
}

void HybridBindingSearchManager::removePossibleBindingsstencil(short idvec[2], CCylinder*
                                    cc, short bindingSite) {
    if(CROSSCHECK_BS_SWITCH) {
        CController::_crosscheckdumpFilechem << "Removing site " << cc->getCylinder()->getId()
                                             << " " << cc->getCylinder()->getStableIndex() << " " << bindingSite
                                             << " with idvec " << idvec[0] << " " << idvec[1] << endl;
    }

    short idx = idvec[0];
    short idx2 = idvec[1];
    auto fIDpair = _filamentIDvec[idx].data();
    short _filamentType = cc->getType();

    if(_filamentType != fIDpair[0] && _filamentType != fIDpair[1] ) return;

    //remove all tuples which have this ccylinder as key
    uint32_t t = cc->getCylinder()->getStableIndex()<<SysParams::Chemistry().shiftbybits;

    uint32_t pos = find (SysParams::Chemistry().bindingSites[_filamentType].begin(),
                      SysParams::Chemistry().bindingSites[_filamentType].end(),
                      bindingSite) -
                      SysParams::Chemistry().bindingSites[_filamentType].begin();
    //Key
    t = t|pos;

    if(CROSSCHECK_BS_SWITCH) {
        CController::_crosscheckdumpFilechem <<"Removing by key"<<endl;
    }
    _possibleBindingsstencilvecuint[idx][idx2].erase(t);

    if(CROSSCHECK_BS_SWITCH) {
        CController::_crosscheckdumpFilechem <<"Removing by value"<<endl;
    }
    //remove all tuples which have this as value
    //Iterate through the reverse map
    auto keys = _reversepossibleBindingsstencilvecuint[idx][idx2][t];//keys that contain
    // t as value in possiblebindings

    for(auto k:keys){
        //get the iterator range that corresponds to this key.
        auto range = _possibleBindingsstencilvecuint[idx][idx2].equal_range(k);
//	    auto range = _mpossibleBindingsstencilvecuint[idx][idx2].equal_range(k);
        //iterate through the range
        for(auto it = range.first; it != range.second;){
            //Go through the value vector and delete entries that match
            it->second.erase(remove(it->second.begin(), it->second.end(), t), it->second.end());
            it++;
        }
    }
    //remove from the reverse map.
	_reversepossibleBindingsstencilvecuint[idx][idx2][t].clear();

    if(CROSSCHECK_BS_SWITCH) {
        CController::_crosscheckdumpFilechem <<"Update rxn"<<endl;
    }
    countNpairsfound(idvec);
    fManagervec[idx][idx2]->updateBindingReaction(Nbindingpairs[idx][idx2]);

    //remove all neighbors which have this binding site pair
    //Go through enclosing compartments to remove all entries with the current cylinder
    // and binding site as values.
    auto nencl = _compartment->getenclosingNeighbours().size();
    for(auto nc: _compartment->getenclosingNeighbours()){
        if(nc != _compartment) {
            auto m = nc->getHybridBindingSearchManager();
	        if(CROSSCHECK_BS_SWITCH)
	            CController::_crosscheckdumpFilechem <<"Remove by value from neighbor "
												""<<nc->getId()<<" total "<<nencl<<endl;

            //Iterate through the reverse map
            auto keys = m->_reversepossibleBindingsstencilvecuint[idx][idx2][t];//keys that
            // contain t as value in possiblebindings

            for(auto k:keys){
                //get the iterator range that corresponds to this key.
                auto range = m->_possibleBindingsstencilvecuint[idx][idx2].equal_range(k);
//	            auto range = m->_mpossibleBindingsstencilvecuint[idx][idx2].equal_range(k);
                //iterate through the range
                for(auto it = range.first; it != range.second;){
                    //Go through the value vector and delete entries that match
                    it->second.erase(remove(it->second.begin(), it->second.end(), t),
                            it->second.end());
	                it++;
                }
            }
            //remove from the reverse map.
            m->_reversepossibleBindingsstencilvecuint[idx][idx2][t].clear();
            if(CROSSCHECK_BS_SWITCH) {
                CController::_crosscheckdumpFilechem <<"Update rxn"<<endl;
            }
            m->countNpairsfound(idvec);
            m->fManagervec[idx][idx2]->updateBindingReaction(m->Nbindingpairs[idx][idx2]);

        }
    }
//	checkoccupancySIMD(idvec);
}

void HybridBindingSearchManager::appendPossibleBindingsstencil(short idvec[2],
																CCylinder* ccyl1,
																CCylinder* ccyl2,
																short site1,
																short site2){
    short idx = idvec[0];
    short idx2 = idvec[1];
    short _filamentType = ccyl1->getType();
    short _nfilamentType = ccyl2->getType();
    uint32_t shiftedIndex1 = ccyl1->getCylinder()->getStableIndex();
    shiftedIndex1 = shiftedIndex1 << SysParams::Chemistry().shiftbybits;
    uint32_t pos1 = find(SysParams::Chemistry().bindingSites[_filamentType].begin(),
                            SysParams::Chemistry().bindingSites[_filamentType].end(), site1)
                    - SysParams::Chemistry().bindingSites[_filamentType].begin();
    uint32_t t1 = shiftedIndex1|pos1;

    uint32_t shiftedIndex2 = ccyl2->getCylinder()->getStableIndex();
    shiftedIndex2 = shiftedIndex2 << SysParams::Chemistry().shiftbybits;
    uint32_t pos2 = find(SysParams::Chemistry().bindingSites[_nfilamentType].begin(),
                            SysParams::Chemistry().bindingSites[_nfilamentType].end(), site2)
                    - SysParams::Chemistry().bindingSites[_nfilamentType].begin();
    uint32_t t2 = shiftedIndex2|pos2;

    _possibleBindingsstencilvecuint[idx][idx2][t1].push_back(t2);
    _reversepossibleBindingsstencilvecuint[idx][idx2][t2].push_back(t1);

    countNpairsfound(idvec);
    fManagervec[idx][idx2]->updateBindingReaction(Nbindingpairs[idx][idx2]);

}

void HybridBindingSearchManager::checkoccupancySIMD(short idvec[2]){

    short idx = idvec[0];
    short idx2 = idvec[1];
    short _filamentType = 0;
    short _complimentaryfilamentType = 0;
    auto speciesboundvec = SysParams::Mechanics().speciesboundvec;

    vector<int> CIDvec(Cylinder::rawNumStableElements());
    const auto& cylinderInfoData = Cylinder::getDbData();
    short bstatepos = bstateposvec[idx][idx2];

    for(auto cyl: Cylinder::getCylinders())
        CIDvec[cyl->getStableIndex()] = cyl->getId();

    auto pbs = _possibleBindingsstencilvecuint[idx][idx2];

    for(auto pair = pbs.begin(); pair != pbs.end(); pair++){

        //Key
        uint32_t leg1 = pair->first;
        vector<uint32_t> leg2 = pair->second;

        uint32_t cIndex1 = leg1 >> SysParams::Chemistry().shiftbybits;
        uint32_t bsite1 = mask & leg1;

        CCylinder* ccyl1 = cylinderInfoData[cIndex1].chemCylinder;

        short it1 = SysParams::Chemistry().bindingSites[_filamentType][bsite1];

        bool boundstate1 = false;

        if(bstatepos == 1) {
            boundstate1 = areEqual(ccyl1->getCMonomer(it1)->speciesBound(
                    SysParams::Chemistry().linkerBoundIndex[_filamentType])->getN(), 1.0);

        }
        if(bstatepos == 2) {
            boundstate1 = areEqual(ccyl1->getCMonomer(it1)->speciesBound(
                    SysParams::Chemistry().motorBoundIndex[_filamentType])->getN(), 1.0);
        }


        //Values
        for(auto V:leg2){
            uint32_t cIndex2 = V >> SysParams::Chemistry().shiftbybits;
            uint32_t bsite2 = mask & V;
            CCylinder* ccyl2 = cylinderInfoData[cIndex2].chemCylinder;

            short it2 = SysParams::Chemistry().bindingSites[_complimentaryfilamentType][bsite2];

            bool boundstate2 = false;
//            bool vecboundstate2 = speciesboundvec[bstatepos][maxnbs * cIndex2 + bsite2];
            if(bstatepos == 1) {
                boundstate2 = areEqual(ccyl2->getCMonomer(it2)->speciesBound(
                        SysParams::Chemistry().linkerBoundIndex[_complimentaryfilamentType])->getN(),
                                       1.0);
            }
            if(bstatepos == 2) {
                boundstate2 = areEqual(ccyl2->getCMonomer(it2)->speciesBound(
                        SysParams::Chemistry().motorBoundIndex[_complimentaryfilamentType])->getN(),
                                       1.0);
            }
			auto ncmpcoord = ccyl2->getCompartment()->coordinates();
            if(boundstate1 == false || boundstate2 == false) {

                std::cout << "OOPS occupied species exist " << boundstate1 << " "
                          << boundstate2 << endl;
                cout<<"bstate pos (1-Linker, 2-Motor)= "<<bstatepos<<endl;

                std::cout<<"Cmp "<<_compartment->coordinates()[0]<<" "<<
                                   _compartment->coordinates()[1]<<" "<<
                                   _compartment->coordinates()[2]<<
                        " nCmp "<<ncmpcoord[0]<<" "<< ncmpcoord[1]<<" "<< ncmpcoord[2]<<
                        " L/M Cylindex bs "<< cIndex1<<" "<<it1<<" "<<
                                  cIndex2<<" "<<it2<<" cIDs "
                                  <<ccyl1->getCylinder()->getId()<<" "
                                  <<ccyl2->getCylinder()->getId()<<" "<<endl;

                cout<<"uint32_t "<<leg1<<" "<<V<<endl;
                exit(EXIT_FAILURE);

            }

        }
    }
}

void HybridBindingSearchManager::updateAllPossibleBindingsstencilHYBD() {

	//Delete all entries in the binding pair maps
    for (int idx = 0; idx < totaluniquefIDpairs; idx++){
        int countbounds = _rMaxsqvec[idx].size();
        for (int idx2 = 0; idx2 < countbounds; idx2++) {
            _possibleBindingsstencilvecuint[idx][idx2].clear();
            _reversepossibleBindingsstencilvecuint[idx][idx2].clear();
        }
    }

    floatingpoint min1,min2,max1,max2;
    bool status1 = true;
    bool status2 = true;
    floatingpoint minveca[2];
    floatingpoint maxveca[2];
    //vector of squared cylinder length
    floatingpoint* cylsqmagnitudevector = SysParams::Mechanics().cylsqmagnitudevector;
    int Ncylincmp = _compartment->getCylinders().size();
    int* cindexvec = new int[Ncylincmp]; //stores cindex of cylinders in this compartment
    vector<vector<int>> ncindices; //cindices of cylinders in neighbor list.
    vector<int> ncindex; //helper vector
	//vector storing information on state (bound/free) of each binding site
    auto boundstate = SysParams::Mechanics().speciesboundvec;
    int maxnbs = SysParams::Chemistry().maxbindingsitespercylinder;

    const auto& cylinderInfoData = Cylinder::getDbData();

    int idx; int idx2;
	//Go through all filament types in our simulation
    for(idx = 0; idx<totaluniquefIDpairs;idx++) {

        long id = 0;
        ncindices.clear();
        auto fpairs = _filamentIDvec[idx].data();

        int nbs1 = SysParams::Chemistry().bindingSites[fpairs[0]].size();
        int nbs2 = SysParams::Chemistry().bindingSites[fpairs[1]].size();
        int Nbounds = _rMaxsqvec[idx].size();

        //Go through cylinders in the compartment, get the half of neighbors.
        for (auto c : _compartment->getCylinders()) {
            cindexvec[id] = c->getStableIndex();
            id++;
            //Get neighors corresponding to the cylinder.
            auto Neighbors = _HneighborList->getNeighborsstencil(HNLIDvec[idx], c);
            ncindex.reserve(Neighbors.size());
            for (auto cn : Neighbors) {
                if (c->getId() > cn->getId())
                    ncindex.push_back(cn->getStableIndex());
            }
            ncindices.push_back(ncindex);
            ncindex.clear();
        }

        //Go through cylinder
        for (int i = 0; i < Ncylincmp; i++) {

            int cindex = cindexvec[i];
            short complimentaryfID;
            const auto& c = cylinderInfoData[cindex];
            const auto& x1 = Bead::getStableElement(c.beadIndices[0])->coord;
            const auto& x2 = Bead::getStableElement(c.beadIndices[1])->coord;
            floatingpoint X1X2[3] = {x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2]};
            int *cnindices = ncindices[i].data();

            if (c.type != fpairs[0] && c.type != fpairs[1]) continue;
            else if(c.type == fpairs[0]) complimentaryfID = fpairs[1];
            else complimentaryfID = fpairs[0];

            //Go through the neighbors of the cylinder
            for (int arraycount = 0; arraycount < ncindices[i].size(); arraycount++) {

                int cnindex = cnindices[arraycount];
                const auto& cn = cylinderInfoData[cnindex];

//            if(c.ID < cn.ID) {counter++; continue;} commented as the above vector does
//              not contain ncs that will fail this cndn.
                if (c.filamentId == cn.filamentId) continue;
                if(c.type != complimentaryfID) continue;

                const auto& x3 = Bead::getStableElement(cn.beadIndices[0])->coord;
                const auto& x4 = Bead::getStableElement(cn.beadIndices[1])->coord;
                floatingpoint X1X3[3] = {x3[0] - x1[0], x3[1] - x1[1], x3[2] - x1[2]};
                floatingpoint X3X4[3] = {x4[0] - x3[0], x4[1] - x3[1], x4[2] - x3[2]};
                floatingpoint X1X3squared = sqmagnitude(X1X3);
                floatingpoint X1X2squared = cylsqmagnitudevector[cindex];
                floatingpoint X1X3dotX1X2 = scalarprojection(X1X3, X1X2);
                floatingpoint X3X4squared = cylsqmagnitudevector[cnindex];
                floatingpoint X1X3dotX3X4 = scalarprojection(X1X3, X3X4);
                floatingpoint X3X4dotX1X2 = scalarprojection(X3X4, X1X2);

                //Number of binding sites on the cylinder/filament of Type A.
                for (int pos1 = 0; pos1 < nbs1; pos1++) {

                    //Number of binding distance pairs we are looking at
                    for (idx2 = 0; idx2 < Nbounds; idx2++) {

                        short bstatepos = bstateposvec[idx][idx2];
						// Check for linkers if the binding site is acceptable.
	                    if (bstatepos == 1) {
		                    int bindingsitestep = SysParams::Chemistry().linkerbindingskip - 1;
		                    if(pos1 % bindingsitestep != 0) continue;
	                    }

                        //now re add valid binding sites
                        if (areEqual(boundstate[bstatepos][maxnbs * cindex + pos1], 1.0)) {

                            auto mp1 = bindingsites1[idx][pos1];
                            floatingpoint A = X3X4squared;
                            floatingpoint B = 2.0 * X1X3dotX3X4 - 2.0 * mp1 * X3X4dotX1X2;
                            floatingpoint C = X1X3squared + mp1 * mp1 * X1X2squared -
                                       2.0 * mp1 * X1X3dotX1X2;
                            floatingpoint Bsq = B*B;
                            floatingpoint C1 = C - _rMinsqvec[idx][idx2];
                            floatingpoint C2 = C - _rMaxsqvec[idx][idx2];
                            floatingpoint b2m4ac1 = Bsq - 4 * A * C1;
                            floatingpoint b2m4ac2 = Bsq - 4 * A * C2;


                            status1 = b2m4ac1 < 0;
                            status2 = b2m4ac2 < 0;

                            if (status1 && status2) continue;

                            if (!status1) {
                                min1 = (-B + sqrt(b2m4ac1)) / (2 * A);
                                min2 = (-B - sqrt(b2m4ac1)) / (2 * A);
                                if (min1 < min2) {
                                    minveca[0] = (min1);
                                    minveca[1] = (min2);
                                } else {
                                    minveca[0] = (min2);
                                    minveca[1] = (min1);
                                }
                                //Compare the MIN solutions are within the first and last
                                // binding sites in filament/cylinder of type B.
                                if (minveca[0] < minparamcyl2[idx][1] &&
                                    minveca[1] > maxparamcyl2[idx][1]) continue;
                            }

                            if (!status2) {
                                max1 = (-B + sqrt(b2m4ac2)) / (2 * A);
                                max2 = (-B - sqrt(b2m4ac2)) / (2 * A);
                                if (max1 < max2) {
                                    maxveca[0] = (max1);
                                    maxveca[1] = (max2);
                                } else {
                                    maxveca[0] = (max2);
                                    maxveca[1] = (max1);
                                }
                                //Compare the mAX solutions are within the first and last
                                // binding sites in filament/cylinder of type B.
                                if (maxveca[0] > maxparamcyl2[idx][1] ||
                                    maxveca[1] < minparamcyl2[idx][1]) continue;
                            }

                            for (int pos2 = 0; pos2 < nbs2; pos2++) {

	                            // Check for linkers if the binding site is acceptable.
	                            if (bstatepos == 1) {
		                            int bindingsitestep = SysParams::Chemistry().linkerbindingskip - 1;
		                            if(pos2 % bindingsitestep != 0) continue;
	                            }

                                if (areEqual(boundstate[bstatepos][maxnbs * cnindex + pos2], 1.0)) {

                                    //check distances..
                                    auto mp2 = bindingsites2[idx][pos2];
                                    if (!status2)
                                        if (mp2 < maxveca[0] || mp2 > maxveca[1]) continue;
                                    if (!status1)
                                        if (mp2 > minveca[0] && mp2 < minveca[1]) continue;


                                    auto it1 = SysParams::Chemistry().bindingSites[fpairs[0]][pos1];
                                    auto it2 = SysParams::Chemistry().bindingSites[fpairs[1]][pos2];

	                                uint32_t shiftedIndex1 = cindex;
	                                shiftedIndex1 = shiftedIndex1 << SysParams::Chemistry().shiftbybits;
	                                uint32_t t1 = shiftedIndex1|pos1;

	                                uint32_t shiftedIndex2 = cnindex << SysParams::Chemistry().shiftbybits;
	                                uint32_t t2 = shiftedIndex2|pos2;

	                                //add in correct order
	                                _possibleBindingsstencilvecuint[idx][idx2][t1].push_back(t2);
	                                _reversepossibleBindingsstencilvecuint[idx][idx2][t2].push_back(t1);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    //Place in the appropriate BindingManager
    /*for(short idx = 0; idx<totaluniquefIDpairs; idx++){
        int countbounds = _rMaxsqvec[idx].size();
        for (short idx2 = 0; idx2 < countbounds; idx2++) {
        	short idvec[2];
        	idvec = {idx, idx2};
	        countNpairsfound(idvec);
	        fManagervec[idx][idx2]->updateBindingReaction(Nbindingpairs[idx][idx2]);
        }
    }*/
    delete[] cindexvec;
}

void HybridBindingSearchManager::countNpairsfound(short idvec[2]){
    short idx = idvec[0];
    short idx2 = idvec[1];
    int N = 0;
    Nbindingpairs[idx][idx2] = 0;

    for (auto iter = _possibleBindingsstencilvecuint[idx][idx2].begin(); iter !=
		    _possibleBindingsstencilvecuint[idx][idx2].end(); iter++) {
        N += iter->second.size();
    }
//    cout<<Nbindingpairs[idx][idx2]<<" "<<N<<endl;
    Nbindingpairs[idx][idx2] = N;
}

void HybridBindingSearchManager::updateAllPossibleBindingsstencilSIMDV3() {

#ifdef SIMDBINDINGSEARCH
    short count = 0;
	for (short idx = 0; idx < totaluniquefIDpairs; idx++) {
		int countbounds = _rMaxsqvec[idx].size();
		for (short idx2 = 0; idx2 < countbounds; idx2++) {
			short idvec[2] = {idx, idx2};
			bool LinkerorMotor = false; //Motor
			if (bstateposvec[idx][idx2] == 1)
				LinkerorMotor = true;//Linker
			if(LinkerorMotor) {
				//Linker
				calculatebspairsLMselfV3<1,true, true>(getdOut<1U, true>(count), idvec);
				calculatebspairsLMenclosedV3<1,false, true>(getdOut<1U, false>(count),
				        bspairslinker2, idvec);
			}
			else{
				//Motor
				calculatebspairsLMselfV3<1,true, false>(getdOut<1U, true>(count), idvec);
				calculatebspairsLMenclosedV3<1,false, false>(getdOut<1U, false>(count),
				        bspairsmotor2, idvec);
			}
			count++;
//			checkoccupancySIMD(idvec);
		}
	}

#endif
}

void HybridBindingSearchManager::updateAllBindingReactions() {
	for (short idx = 0; idx < totaluniquefIDpairs; idx++) {
		int countbounds = _rMaxsqvec[idx].size();
		for (short idx2 = 0; idx2 < countbounds; idx2++) {
			short idvec[2] = {idx, idx2};
			countNpairsfound(idvec);
			fManagervec[idx][idx2]->updateBindingReaction(Nbindingpairs[idx][idx2]);
		}
	}
}

void HybridBindingSearchManager::addtoHNeighborList(){
    HNLIDvec.resize(totaluniquefIDpairs);
    for(short idx = 0; idx < totaluniquefIDpairs; idx++) {
        vector<float> temprMaxsq = _rMaxsqvec[idx];
        vector<float> temprMinsq = _rMinsqvec[idx];
        vector<short> ftypepair = _filamentIDvec[idx];
        float localmaxcylsize = max(SysParams::Geometry().cylinderSize[ftypepair[0]],
                                    SysParams::Geometry().cylinderSize[ftypepair[1]]);
        float maxrMaxsq = *max_element(temprMaxsq.begin(), temprMaxsq.end());
        float minrMinsq = *min_element(temprMinsq.begin(), temprMinsq.end());
	    HNLIDvec[idx] = _HneighborList->setneighborsearchparameters(ftypepair[0], ftypepair[1],
		                                                   /*uniquestatus*/false,
		                                                     /*fullstatus*/ true,
		                                              localmaxcylsize + sqrt(maxrMaxsq),
	                                                       max( sqrt(minrMinsq) - localmaxcylsize,float(0.0)));

        /*for(short idx2 = 0; idx2<temprMaxsq.size();idx2++) {
            short idvec[2] = {idx, idx2};
            fManagervec[idx][idx2]->setHNLID(HNLIDvec[idx], idvec);
        }*/
    }
}

void HybridBindingSearchManager::copyInfotoBindingManagers() {
	for (short idx = 0; idx < totaluniquefIDpairs; idx++) {
		vector<float> temprMaxsq = _rMaxsqvec[idx];
		for (short idx2 = 0; idx2 < temprMaxsq.size(); idx2++) {
			short idvec[2] = {idx, idx2};
			fManagervec[idx][idx2]->setHNLID(HNLIDvec[idx], idvec);
		}
	}
}

vector<tuple<CCylinder*, short>>
HybridBindingSearchManager::chooseBindingSitesstencil(short idvec[2]){

	if(CROSSCHECK_BS_SWITCH)
	    CController::_crosscheckdumpFilechem <<"Choosing site"<<endl;

    short idx = idvec[0];
    short idx2 = idvec[1];
    auto fpairs = _filamentIDvec[idx].data();
    int pbsSize = Nbindingpairs[idx][idx2];

    const auto& cylinderInfoData = Cylinder::getDbData();

    assert((pbsSize!= 0)
           && "Major bug: Linker/Motor binding manager should not have zero binding \
                   sites when called to choose a binding site.");
    if(true) {
	    uint randomIndex = Rand::randInteger(1, pbsSize);

	    auto it = _possibleBindingsstencilvecuint[idx][idx2].begin();
	    uint cumulativesum = 0;
	    uint prevstepsum = 0;

	    while (it != _possibleBindingsstencilvecuint[idx][idx2].end() &&
	           cumulativesum < randomIndex) {
		    prevstepsum = cumulativesum;
		    cumulativesum += it->second.size();
		    if (cumulativesum < randomIndex)
			    it++;
	    }
	    Nbindingpairs[idx][idx2]--;

	    uint position = randomIndex - prevstepsum - 1;

	    uint32_t site1 = it->first;
	    uint32_t site2 = it->second[position];

	    uint32_t cIndex1 = site1 >> SysParams::Chemistry().shiftbybits;
	    uint32_t cIndex2 = site2 >> SysParams::Chemistry().shiftbybits;

	    short bsitepos1 = mask & site1;
	    short bsitepos2 = mask & site2;

	    if(CROSSCHECK_BS_SWITCH)
	        CController::_crosscheckdumpFilechem <<"Chosen cindices, pos "<<cIndex1<<" "
	                <<cIndex2<<" "<<bsitepos1<<" "<<bsitepos2<<endl;

	    CCylinder *ccyl1;
	    CCylinder *ccyl2;

	    ccyl1 = cylinderInfoData[cIndex1].chemCylinder;
	    ccyl2 = cylinderInfoData[cIndex2].chemCylinder;

	    short bindingSite1 = SysParams::Chemistry().bindingSites[ccyl1->getType
			    ()][bsitepos1];
	    short bindingSite2 = SysParams::Chemistry().bindingSites[ccyl2->getType
			    ()][bsitepos2];

	    tuple<CCylinder *, short> t1 = make_tuple(ccyl1, bindingSite1);
	    tuple<CCylinder *, short> t2 = make_tuple(ccyl2, bindingSite2);

	    if(CROSSCHECK_BS_SWITCH)
	        CController::_crosscheckdumpFilechem <<"Chosen!"<<endl;
	    return vector<tuple<CCylinder *, short>>{t1, t2};
    }
}

void HybridBindingSearchManager::clearPossibleBindingsstencil(short idvec[2]){
	short idx = idvec[0];
	short idx2 = idvec[1];
	_possibleBindingsstencilvecuint[idx][idx2].clear();
	_reversepossibleBindingsstencilvecuint[idx][idx2].clear();
	countNpairsfound(idvec);
	fManagervec[idx][idx2]->updateBindingReaction(Nbindingpairs[idx][idx2]);
}

void HybridBindingSearchManager::printbindingsitesstencil(short idvec[2]) {
    cout<<"BINDINGSITES: CYL1(SIDX) CYL2(SIDX) SITE1 SITE2"<<endl;
    short idx = idvec[0];
    short idx2 = idvec[1];

    auto pbs = _possibleBindingsstencilvecuint[idx][idx2];

    for (auto pair = pbs.begin(); pair != pbs.end(); pair++) {

        //Key
        uint32_t leg1 = pair->first;
        vector<uint32_t> leg2 = pair->second;

        uint32_t cIndex1 = leg1 >> SysParams::Chemistry().shiftbybits;
        uint32_t bsite1 = mask & leg1;
        //Values
        for (auto V:leg2) {
            uint32_t cIndex2 = V >> SysParams::Chemistry().shiftbybits;
            uint32_t bsite2 = mask & V;
            cout<<cIndex1<<" "<<cIndex2<<" "<<bsite1<<" "<<bsite2<<endl;
        }
    }
}

HybridCylinderCylinderNL* HybridBindingSearchManager::_HneighborList;
bool initialized = false;
#ifdef SIMDBINDINGSEARCH
//D = 1
//dist::dOut<1U,false> HybridBindingSearchManager::bspairslinker;
dist::dOut<1U,true> HybridBindingSearchManager::bspairslinkerself;
//dist::dOut<1U,false> HybridBindingSearchManager::bspairsmotor;
dist::dOut<1U,true> HybridBindingSearchManager::bspairsmotorself;
dist::dOut<1U,false> HybridBindingSearchManager::bspairsmotor2;
dist::dOut<1U,false> HybridBindingSearchManager::bspairslinker2;

floatingpoint HybridBindingSearchManager::largestlinkerdistance = 0.0;
floatingpoint HybridBindingSearchManager::largestmotordistance = 0.0;
//D = 2
/*dist::dOut<2U,true> HybridBindingSearchManager::bspairs2self;
dist::dOut<2U,false> HybridBindingSearchManager::bspairs2;
dist::dOut<1U,false> HybridBindingSearchManager::bspairs2_D1;
dist::dOut<1U,false> HybridBindingSearchManager::bspairs2_D2;*/


short HybridBindingSearchManager::Totallinkermotor = 0;
#endif
floatingpoint HybridBindingSearchManager::SIMDtime = 0.0;
floatingpoint HybridBindingSearchManager::HYBDtime = 0.0;
floatingpoint HybridBindingSearchManager::findtime = 0.0;
floatingpoint HybridBindingSearchManager::appendtime = 0.0;
floatingpoint HybridBindingSearchManager::findtimeV2 = 0.0;
floatingpoint HybridBindingSearchManager::SIMDparse1 = 0.0;
floatingpoint HybridBindingSearchManager::SIMDparse2 = 0.0;
floatingpoint HybridBindingSearchManager::SIMDparse3 = 0.0;
floatingpoint HybridBindingSearchManager::SIMDcountbs = 0.0;
floatingpoint HybridBindingSearchManager::HYBDappendtime = 0.0;
floatingpoint HybridBindingSearchManager::SIMDV3appendtime = 0.0;
floatingpoint HybridBindingSearchManager::findtimeV3 = 0.0;
#endif

} // namespace medyan

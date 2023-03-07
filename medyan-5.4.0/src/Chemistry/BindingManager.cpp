
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
#include "HybridBindingSearchManager.h"
#include "BindingManager.h"

#include "Compartment.h"
#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"
#include "Structure/SurfaceMesh/Membrane.hpp"

#include "MotorGhost.h"

#include "SubSystem.h"
#include "Boundary.h"
#include "CompartmentGrid.h"

#include "ChemCallbacks.h"
#include "MathFunctions.h"
#include "Controller/GController.h"
#include "SysParams.h"
#include "CUDAcommon.h"
#include <algorithm>
#include "SubSystem.h"



namespace medyan {
using namespace mathfunc;

//BRANCHER

BranchingManager::BranchingManager(ReactionBase* reaction,
                                   Compartment* compartment,
                                   short boundInt,
                                   vector<short> filamentIDvec,
                                   NucleationZoneType zone, floatingpoint nucleationDistance)

: FilamentBindingManager(reaction, compartment, boundInt, filamentIDvec),
_nucleationZone(zone), _nucleationDistance(nucleationDistance) {

    //find the single binding species
    RSpecies** rs = reaction->rspecies();
    string name = rs[B_RXN_INDEX]->getSpecies().getName();

    _bindingSpecies = _compartment->findSpeciesByName(name);

}
#ifdef NLORIGINAL
void BranchingManager::addPossibleBindings(CCylinder* cc, short bindingSite) {

    if (SysParams::INITIALIZEDSTATUS ) {
        short complimentaryfID;
        short _filamentType = cc->getType();
        if (_filamentType != _filamentIDvec[0] && _filamentType != _filamentIDvec[1])
            return;
        else if (_filamentType == _filamentIDvec[0]) complimentaryfID = _filamentIDvec[1];
        else complimentaryfID = _filamentIDvec[0];

        bool inZone = true;
        //see if in nucleation zone
        if (_nucleationZone != NucleationZoneType::ALL) {

            auto mp = (float) bindingSite /
                      SysParams::Geometry().cylinderNumMon[_filamentType];

            auto x1 = cc->getCylinder()->getFirstBead()->vcoordinate();
            auto x2 = cc->getCylinder()->getSecondBead()->vcoordinate();

            auto coord = midPointCoordinate(x1, x2, mp);

        //set nucleation zone
        // For membrane acting as boundaries, only the 0th membrane will be considered.
        if(_nucleationZone == NucleationZoneType::MEMBRANE) {
            if(_subSystem->membranes.size()) {
                if(cc->getCompartment()->isActivated()) {
                    if(cc->getCompartment()->getVolumeFrac() < 1.0) // Not fully activated
                        if(!medyan::contains(*_subSystem, _subSystem->membranes.begin()->getMesh(), vector2Vec<3, floatingpoint>(coord)))
                            inZone = false;
                }
                else inZone = false;
            } // else no membrane exists, always "in zone".
        }
        else if(_subSystem->getBoundary()->distance(coord) < _nucleationDistance) {

            //set nucleation zone
            if (_subSystem->getBoundary()->distance(coord) < _nucleationDistance) {

                //if top boundary, check if we are above the center coordinate in z
                if (_nucleationZone == NucleationZoneType::TOPBOUNDARY) {

                    if (coord[2] >= GController::getCenter()[2])
                        inZone = true;
                    else
                        inZone = false;
                }
                    //add SIDEBOUNDARY that check the distance to the side of cylinder
                else if (_nucleationZone == NucleationZoneType::SIDEBOUNDARY) {
                    if (_subSystem->getBoundary()->sidedistance(coord) <
                        _nucleationDistance)
                        inZone = true;
                    else
                        inZone = false;
                } else inZone = true;
            } else
                inZone = false;
        }

        //add valid site
        if (areEqual(cc->getCMonomer(bindingSite)->speciesBound(
                SysParams::Chemistry().brancherBoundIndex[_filamentType])->getN(),
                     (floatingpoint) 1.0) &&
            inZone) {

            auto t = tuple<CCylinder *, short>(cc, bindingSite);
            _possibleBindings.insert(t);
        }

        int oldN = _bindingSpecies->getN();
        int newN = numBindingSites();

        updateBindingReaction(oldN, newN);
    }
}

void BranchingManager::addPossibleBindings(CCylinder* cc) {

    short _filamentType = cc->getType();
    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
        bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++) {
#ifdef NLORIGINAL
        addPossibleBindings(cc, *bit);
#endif
#if defined(NLSTENCILLIST) || defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
        addPossibleBindingsstencil(cc, *bit);
#endif
    }
}

void BranchingManager::removePossibleBindings(CCylinder* cc, short bindingSite) {

    short complimentaryfID;
    short _filamentType = cc->getType();
    if (_filamentType != _filamentIDvec[0] && _filamentType != _filamentIDvec[1]) return;

    //remove tuple which has this ccylinder
    _possibleBindings.erase(tuple<CCylinder*, short>(cc, bindingSite));

    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();

    updateBindingReaction(oldN, newN);
}

void BranchingManager::removePossibleBindings(CCylinder* cc) {

    short _filamentType = cc->getType();

    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
        bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++) {
#ifdef NLORIGINAL
        removePossibleBindings(cc, *bit);
#endif
    }
}

void BranchingManager::updateAllPossibleBindings() {


    //clear all
    _possibleBindings.clear();
    auto boundstate = SysParams::Mechanics().speciesboundvec;
    int offset = 0;

    for(auto &c : _compartment->getCylinders()) {

        short _filamentType = c->getType();
        if (_filamentType != _filamentIDvec[0] && _filamentType != _filamentIDvec[1])
        	return;

        auto cc = c->getCCylinder();
        int j = -1;
        //now re add valid binding sites
        for(auto it = SysParams::Chemistry().bindingSites[_filamentType].begin();
            it != SysParams::Chemistry().bindingSites[_filamentType].end(); it++) {
            j++;
            bool inZone = true;
            //see if in nucleation zone
            if(_nucleationZone != NucleationZoneType::ALL) {

                auto mp = (float)*it / SysParams::Geometry().cylinderNumMon[_filamentType];

                auto x1 = cc->getCylinder()->getFirstBead()->vcoordinate();
                auto x2 = cc->getCylinder()->getSecondBead()->vcoordinate();

                auto coord = midPointCoordinate(x1, x2, mp);

                //set nucleation zone
                // For membrane acting as boundaries, only the 0th membrane will be considered.
                if(_nucleationZone == NucleationZoneType::MEMBRANE) {
                    if(_subSystem->membranes.size()) {
                        if(cc->getCompartment()->isActivated()) {
                            if(cc->getCompartment()->getVolumeFrac() < 1.0) // Not fully activated
                                if(!medyan::contains(*_subSystem, _subSystem->membranes.begin()->getMesh(), vector2Vec<3, floatingpoint>(coord)))
                                    inZone = false;
                        }
                        else inZone = false;
                    } // else no membrane exists, always "in zone".
                }
                else if(_subSystem->getBoundary()->distance(coord) < _nucleationDistance) {
                    
                    //if top boundary, check if we are above the center coordinate in z
                    if(_nucleationZone == NucleationZoneType::TOPBOUNDARY) {

                        if(coord[2] >= GController::getCenter()[2])
                            inZone = true;
                        else
                            inZone = false;
                    }
                    //add SIDEBOUNDARY that check the distance to the side of cylinder
                    else if(_nucleationZone == NucleationZoneType::SIDEBOUNDARY){
                        if(_subSystem->getBoundary()->sidedistance(coord) < _nucleationDistance){
                            inZone = true;
                            //cout << "x= " << coord[1] << "y= " << coord[2] << endl;
                        }


                        else
                            inZone = false;
                    }
                    else inZone = true;
                }
                else
                    inZone = false;
            }
            if (areEqual(boundstate[0][offset + SysParams::Chemistry()
                                       .bindingSites[_filamentType]
                                       .size()*c->getStableIndex() + j], 1.0) && inZone) {

                auto t = tuple<CCylinder*, short>(cc, *it);
                _possibleBindings.insert(t);
            }
        }
    }
//        std::cout<<_possibleBindings.size()<<endl;
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();

    updateBindingReaction(oldN, newN);
}

void BranchingManager::appendpossibleBindings(
    CCylinder* ccyl1, CCylinder* ccyl2,
    short site1, short site2
) {
    floatingpoint oldN=numBindingSites();
    auto t1 = make_tuple(ccyl1, site1);
    auto t2 = make_tuple(ccyl2, site2);
    _possibleBindings.insert(t1);
    _branchrestarttuple.push_back(make_tuple(t1,t2));
    floatingpoint newN=numBindingSites();
    updateBindingReaction(oldN,newN);}

    void BranchingManager::printbindingsites(){
    cout<<"BINDINGSITES: CYL1(SIDX) SITE1"<<endl;
    for(auto it2 = _possibleBindings.begin(); it2!=_possibleBindings.end(); it2++) {
        auto cyl1 = get<0>(*it2)->getCylinder();
        auto bs1 = get<1>(*it2);
        cout<<cyl1->getStableIndex()<<" "<<cyl2->getStableIndex()<<" "<<pos1<<endl;
    }

}
#endif
bool BranchingManager::isConsistent() {
#ifdef NLORIGINAL
    auto bindinglist = _possibleBindings;
#elif defined(NLSTENCILLIST) || defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
    auto bindinglist = _possibleBindingsstencil;
#endif
    for (auto it = bindinglist.begin(); it != bindinglist.end(); it++) {

        CCylinder* cc = get<0>(*it);
        Cylinder* c   = cc->getCylinder();

        short bindingSite = get<1>(*it);

        bool flag = true;
        short _filamentType = cc->getType();
        if(_filamentType != _filamentIDvec[0] && _filamentType != _filamentIDvec[1]){
            flag = false;
        }
        //check site empty
        if(!areEqual(cc->getCMonomer(bindingSite)->speciesBound(
                SysParams::Chemistry().brancherBoundIndex[_filamentType])->getN(), (floatingpoint)1.0))
            flag = false;

        if(!flag) {
            cout << "Binding site in branching manager is inconsistent. " << endl;
            cout << "Binding site = " << bindingSite << endl;

            cout << "Cylinder info ..." << endl;
            c->printSelf();

            return false;
        }
    }
    return true;
}

#if defined(NLSTENCILLIST) || defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
void BranchingManager::addPossibleBindingsstencil(CCylinder* cc) {

  short _filamentType = cc->getType();

    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
        bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++) {
        addPossibleBindingsstencil(cc, *bit);
    }
}
void BranchingManager::addPossibleBindingsstencil(CCylinder* cc, short bindingSite) {

	if (SysParams::INITIALIZEDSTATUS ) {
    short complimentaryfID;
    short _filamentType = cc->getType();
    if (_filamentType != _filamentIDvec[0] && _filamentType != _filamentIDvec[1]) return;
    else if (_filamentType == _filamentIDvec[0]) complimentaryfID = _filamentIDvec[1];
    else complimentaryfID = _filamentIDvec[0];

    bool inZone = true;
    //see if in nucleation zone
    if(_nucleationZone != NucleationZoneType::ALL) {

        auto mp = (float)bindingSite / SysParams::Geometry().cylinderNumMon[_filamentType];

        auto x1 = cc->getCylinder()->getFirstBead()->vcoordinate();
        auto x2 = cc->getCylinder()->getSecondBead()->vcoordinate();

        auto coord = midPointCoordinate(x1, x2, mp);

        //set nucleation zone
        if(_subSystem->getBoundary()->distance(coord) < _nucleationDistance) {

            //if top boundary, check if we are above the center coordinate in z
            if(_nucleationZone == NucleationZoneType::TOPBOUNDARY) {

                if(coord[2] >= GController::getCenter()[2])
                    inZone = true;
                else
                    inZone = false;
            }
                //add SIDEBOUNDARY that check the distance to the side of cylinder
            else if(_nucleationZone == NucleationZoneType::SIDEBOUNDARY){
                if(_subSystem->getBoundary()->sidedistance(coord) < _nucleationDistance)
                    inZone = true;
                else
                    inZone = false;
            }
            else if(_nucleationZone == NucleationZoneType::RIGHTBOUNDARY){
                floatingpoint dis = _subSystem->getBoundary()->getboundaryelementcoord(1) -
                        coord[0];
                if(dis > 0 && dis < _nucleationDistance)
                    inZone = true;
                else
                    inZone = false;
            }

            else inZone = true;
        }
        else
            inZone = false;
    }

    //add valid site
    if (areEqual(cc->getCMonomer(bindingSite)->speciesBound(
            SysParams::Chemistry().brancherBoundIndex[_filamentType])->getN(),
                    (floatingpoint)1.0) && inZone) {

        auto t = tuple<CCylinder*, short>(cc, bindingSite);
        _possibleBindingsstencil.insert(t);
    }

    int oldN = _bindingSpecies->getN();
    int newN = numBindingSitesstencil();

    updateBindingReaction(oldN, newN);
    }
}
void BranchingManager::updateAllPossibleBindingsstencil() {
    //clear all
    _possibleBindingsstencil.clear();
    auto boundstate = SysParams::Mechanics().speciesboundvec;
    int maxnbs = SysParams::Chemistry().maxbindingsitespercylinder;

    for(auto &c : _compartment->getCylinders()) {

        short _filamentType = c->getType();
        if (_filamentType != _filamentIDvec[0] && _filamentType != _filamentIDvec[1])
        	return;

        auto cc = c->getCCylinder();
        int j = -1;
        //now re add valid binding sites
        for(auto it = SysParams::Chemistry().bindingSites[_filamentType].begin();
            it != SysParams::Chemistry().bindingSites[_filamentType].end(); it++) {
            j++;
            bool inZone = true;
            //see if in nucleation zone
            if(_nucleationZone != NucleationZoneType::ALL) {

                auto mp = (float)*it / SysParams::Geometry().cylinderNumMon[_filamentType];

                auto x1 = cc->getCylinder()->getFirstBead()->vcoordinate();
                auto x2 = cc->getCylinder()->getSecondBead()->vcoordinate();

                auto coord = midPointCoordinate(x1, x2, mp);

                //set nucleation zone
                if(_subSystem->getBoundary()->distance(coord) < _nucleationDistance) {

                    //if top boundary, check if we are above the center coordinate in z
                    if(_nucleationZone == NucleationZoneType::TOPBOUNDARY) {

                        if(coord[2] >= GController::getCenter()[2])
                            inZone = true;
                        else
                            inZone = false;
                    }
                    //add SIDEBOUNDARY that check the distance to the side of cylinder
                    else if(_nucleationZone == NucleationZoneType::SIDEBOUNDARY){
                        if(_subSystem->getBoundary()->sidedistance(coord) < _nucleationDistance){
                            inZone = true;
                            //cout << "x= " << coord[1] << "y= " << coord[2] << endl;
                        }

                        else
                            inZone = false;
                    }
                    else if(_nucleationZone == NucleationZoneType::RIGHTBOUNDARY){
                        floatingpoint dis = _subSystem->getBoundary()->getboundaryelementcoord(1) -
                                     coord[0];
                        if(dis > 0 && dis < _nucleationDistance)
                            inZone = true;
                        else
                            inZone = false;
                    }

                    else inZone = true;
                }
                else
                    inZone = false;
            }
            if (areEqual(boundstate[0][maxnbs * c->getStableIndex() + j], (floatingpoint)1.0) && inZone) {
//                output test
//                auto mp = (float)*it / SysParams::Geometry().cylinderNumMon[_filamentType];
//                auto x1 = cc->getCylinder()->getFirstBead()->vcoordinate();
//                auto x2 = cc->getCylinder()->getSecondBead()->vcoordinate();
//
//                auto coord = midPointCoordinate(x1, x2, mp);
//                std::cout<<c->_dcIndex<<" "<<*it<<" "<<_subSystem->getBoundary()->distance(coord)<<endl;
//                end
                auto t = tuple<CCylinder*, short>(cc, *it);
                _possibleBindingsstencil.insert(t);
            }
        }
    }
    //        std::cout<<_possibleBindings.size()<<endl;
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSitesstencil();

    updateBindingReaction(oldN, newN);
    /*std::cout<<"Branching consistency "<<isConsistent()<<endl;*/
}
void BranchingManager::appendpossibleBindingsstencil(short boundInt, CCylinder* ccyl1, CCylinder* ccyl2, short site1,
                                   short site2){
    if(boundInt != _boundInt) return;
	floatingpoint oldN=numBindingSitesstencil();
	auto t1 = make_tuple(ccyl1, site1);
	auto t2 = make_tuple(ccyl2, site2);
	_possibleBindingsstencil.insert(t1);
	_branchrestarttuple.push_back(make_tuple(t1,t2));
	floatingpoint newN=numBindingSitesstencil();
	updateBindingReaction(oldN,newN);
}
void BranchingManager::removePossibleBindingsstencil(CCylinder* cc) {

   short _filamentType = cc->getType();

    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
        bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++)

    removePossibleBindingsstencil(cc, *bit);
}
void BranchingManager::removePossibleBindingsstencil(CCylinder* cc, short bindingSite) {

  short _filamentType = cc->getType();
  if (_filamentType != _filamentIDvec[0] && _filamentType != _filamentIDvec[1]) return;

    //remove tuple which has this ccylinder
    _possibleBindingsstencil.erase(tuple<CCylinder*, short>(cc, bindingSite));

    int oldN = _bindingSpecies->getN();
    int newN = numBindingSitesstencil();

    updateBindingReaction(oldN, newN);
}
void BranchingManager::crosscheck(){
    //cout<<"Branching NLORIGINAL size "<<_possibleBindings.size()<<" NLSTENCIL size "
    //       <<_possibleBindingsstencil.size()<<endl;
    if(_possibleBindings.size() != _possibleBindingsstencil.size())
    cout<<"Branching.. The two methods compared do not yield the same number of "
    "binding sites"<<endl;
    short matches = 0;

    for(auto it1 = _possibleBindings.begin();it1!=_possibleBindings.end();it1++){
        short matchperelement = 0;

        auto cyl1_o = get<0>(*it1)->getCylinder();
        auto bs1_o =  get<1>(*it1);

        short sum = 0;
        for(auto it2 = _possibleBindingsstencil.begin();it2!=_possibleBindingsstencil.end();
            it2++){
            auto cyl1_s = get<0>(*it2)->getCylinder();
            auto bs1_s =  get<1>(*it2);
            sum = 0;
            if(cyl1_o->getId() == cyl1_s->getId() )
                sum++;

            if(bs1_o == bs1_s )
                sum++;

            if(sum == 2) {
//                cout << "match found" << endl;
                if(matchperelement == 0 ) {
                    matches++;
                    matchperelement++;
                } else {
                    cout << "ERROR. Multiple matches for chosen binding site pair in "
                    "STENCILLIST. Check stencillist.." << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
        }
    std::cout<<"Branching possible bindings size "<<_possibleBindings.size()<<" Total "
            "matches "<<
                                                                                     matches<<endl;
    if(_possibleBindings.size() != matches || _possibleBindings.size() !=
       _possibleBindingsstencil.size()){
        cout<<"Branching. All binding site pairs are not found in Stencil list"<<endl;
        exit(EXIT_FAILURE);
    }
}
void BranchingManager::printbindingsitesstencil(){
	cout<<"BINDINGSITES: CYL1(SIDX) SITE1"<<endl;
    for(auto it2 = _possibleBindingsstencil.begin();it2!=_possibleBindingsstencil.end();
        it2++) {
        auto cyl1 = get<0>(*it2)->getCylinder();
        auto bs1 = get<1>(*it2);
        cout<<cyl1->getStableIndex()<<" "<<bs1<<endl;
    }
    cout<<"BRANCHINGRESTARTTUPLE"<<endl;
    for(auto it2 = _branchrestarttuple.begin();it2!=_branchrestarttuple.end();
        it2++) {
        auto t1 = get<0>(*it2);
        auto t2 = get<1>(*it2);
        auto cyl1 = get<0>(t1)->getCylinder();
        auto cyl2 = get<0>(t2)->getCylinder();
        auto site1 = get<1>(t1);
        auto site2 = get<1>(t2);
        cout<<cyl1->getStableIndex()<<" "<<site1<<" "<<cyl2->getStableIndex()<<" "
                                                                               ""<<site2<<endl;
    }

}
#endif
#ifdef CUDAACCL_NL
void BranchingManager::assigncudavars() {
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_distance, 2 * sizeof(floatingpoint)),"cuda data transfer", " "
            "BindingManager.cu");
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_numpairs, sizeof(int)),"cuda data transfer", " "
                            "BindingManager.cu");
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_zone, sizeof(int)),"cuda data transfer", " "
            "BindingManager.cu");
    floatingpoint dist[2];
    dist[0] = _nucleationDistance;
    dist[1] = GController::getCenter()[2];
    int n[1];
    n[0] = 0;
    int zone[1];
    zone[0] = -1;
    if(_nucleationZone == NucleationZoneType::ALL)
        zone[0] =0;
    else if(_nucleationZone == NucleationZoneType::BOUNDARY)
        zone[0] =1;
    else if(_nucleationZone == NucleationZoneType::TOPBOUNDARY)
        zone[0] =2;
    CUDAcommon::handleerror(cudaMemcpy(gpu_distance, dist, 2 * sizeof(floatingpoint), cudaMemcpyHostToDevice));
    CUDAcommon::handleerror(cudaMemcpy(gpu_numpairs, n, sizeof(int), cudaMemcpyHostToDevice));
    CUDAcommon::handleerror(cudaMemcpy(gpu_zone, zone, sizeof(int), cudaMemcpyHostToDevice));
    //    delete dist;
}

void BranchingManager::freecudavars() {
    CUDAcommon::handleerror(cudaFree(gpu_distance),"cudaFree", "BindingManager");
    CUDAcommon::handleerror(cudaFree(gpu_numpairs),"cudaFree", "BindingManager");
    CUDAcommon::handleerror(cudaFree(gpu_zone),"cudaFree", "BindingManager");
}


int* BranchingManager::getzoneCUDA(){
    return gpu_zone;
}
int* BranchingManager::getnumpairsCUDA(){
    return gpu_numpairs;
}
#endif

//LINKER
LinkerBindingManager::LinkerBindingManager(
    ReactionBase* reaction,
    Compartment* compartment,
    short linkerType,
    vector<short> filamentIDvec,
    int linkerSpeciesIndex1,
    int linkerSpeciesIndex2,
    float rMax, float rMin)

: FilamentBindingManager(reaction, compartment, linkerType, filamentIDvec),
    _rMin(rMin), _rMax(rMax),
    linkerSpeciesIndices_ { linkerSpeciesIndex1, linkerSpeciesIndex2 }
{
    _rMinsq =_rMin * _rMin;
    _rMaxsq = _rMax * _rMax;

    //find the pair binding species
    RSpecies** rs = reaction->rspecies();
    string name = rs[ML_RXN_INDEX]->getSpecies().getName();
    _bindingSpecies = _compartment->findSpeciesByName(name);
    _rMaxsq = rMax*rMax;
    _rMinsq = rMin*rMin;
    for(auto it1 = SysParams::Chemistry().bindingSites[filamentIDvec[0]].begin();
        it1 != SysParams::Chemistry().bindingSites[filamentIDvec[0]].end(); it1++) {
        bindingsites1.push_back((float)*it1 / SysParams::Geometry()
                .cylinderNumMon[filamentIDvec[0]]);
    }
    for(auto it1 = SysParams::Chemistry().bindingSites[filamentIDvec[1]].begin();
        it1 != SysParams::Chemistry().bindingSites[filamentIDvec[1]].end(); it1++) {
        bindingsites2.push_back((float)*it1 / SysParams::Geometry()
                .cylinderNumMon[filamentIDvec[1]]);
    }
}

#ifdef NLORIGINAL
void LinkerBindingManager::addPossibleBindings(CCylinder* cc, short bindingSite) {

    if (SysParams::INITIALIZEDSTATUS ) {
        short complimentaryfID;
        short _filamentType = cc->getType();
        if (_filamentType != _filamentIDvec[0] && _filamentType != _filamentIDvec[1])
            return;
        else if (_filamentType == _filamentIDvec[0]) complimentaryfID = _filamentIDvec[1];
        else complimentaryfID = _filamentIDvec[0];

        //if we change other managers copy number
        vector<LinkerBindingManager *> affectedManagers;

        //add valid binding sites
        if (areEqual(cc->getCMonomer(bindingSite)->speciesBound(SysParams::Chemistry()
                                                                        .linkerBoundIndex[_filamentType])->getN(),
                        (floatingpoint) 1.0)) {

            //loop through neighbors
            //now re add valid based on CCNL
            vector<Cylinder *> nList = _neighborLists[_nlIndex]->getNeighbors(
                    cc->getCylinder());
            for (auto cn : nList) {
                Cylinder *c = cc->getCylinder();

                if (cn->getParent() == c->getParent()) continue;
                short _nfilamentType = cn->getType();
                if (_nfilamentType != complimentaryfID) return;
                if (c->getId() < cn->getId()) continue;

                auto ccn = cn->getCCylinder();
                dBInt = 2;
                for (auto it = SysParams::Chemistry().bindingSites[_nfilamentType].begin();
                        it !=
                        SysParams::Chemistry().bindingSites[_nfilamentType].end(); it++) {
                    // DifBind
                    if (dBInt % dBI != 0) {
                        dBInt += 1;
                        continue;
                    } else {
                        dBInt = 1;
                    }

                    if (areEqual(ccn->getCMonomer(*it)->speciesBound(
                            SysParams::Chemistry().linkerBoundIndex[_nfilamentType])->getN(),
                                    (floatingpoint) 1.0)) {

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

                        floatingpoint distSq = twoPointDistancesquared(m1, m2);

                        if (distSq > _rMaxsq || distSq < _rMinsq) continue;

                        auto t1 = tuple<CCylinder *, short>(cc, bindingSite);
                        auto t2 = tuple<CCylinder *, short>(ccn, *it);

                        //add in correct order
                        if (c->getId() > cn->getId()) {
                            _possibleBindings.emplace(t1, t2);
                            _reversePossibleBindings[t2].push_back(t1);
                        } else {
                            //add in this compartment
                            if (cn->getCompartment() == _compartment) {

                                _possibleBindings.emplace(t1, t2);
                                _reversePossibleBindings[t2].push_back(t1);
                            }
                                //add in other
                            else {
                                auto m = (LinkerBindingManager *) cn->getCompartment()->
                                        getFilamentBindingManagers()[_mIndex].get();

                                affectedManagers.push_back(m);

                                m->_possibleBindings.emplace(t2, t1);
                                m->_reversePossibleBindings[t1].push_back(t2);
                            }
                        }
                    }
                }
            }
        }

        //update affected
        for (auto m : affectedManagers) {

            int oldNOther = m->_bindingSpecies->getN();
            int newNOther = m->numBindingSites();

            m->updateBindingReaction(oldNOther, newNOther);
        }

        //update this manager
        int oldN = _bindingSpecies->getN();
        int newN = numBindingSites();

        updateBindingReaction(oldN, newN);
    }
}

void LinkerBindingManager::addPossibleBindings(CCylinder* cc) {

    short _filamentType = cc->getType();
    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
        bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++){
#ifdef NLORIGINAL
        addPossibleBindings(cc, *bit);
#endif
#if defined(NLSTENCILLIST) || defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
        addPossibleBindingsstencil(cc, *bit);
#endif
    }
}

void LinkerBindingManager::removePossibleBindings(CCylinder* cc, short bindingSite) {

    short complimentaryfID;
    short _filamentType = cc->getType();
    if (_filamentType != _filamentIDvec[0] && _filamentType != _filamentIDvec[1]) return;
    else if (_filamentType == _filamentIDvec[0]) complimentaryfID = _filamentIDvec[1];
    else complimentaryfID = _filamentIDvec[0];

    //if we change other managers copy number
    vector<LinkerBindingManager*> affectedManagers;

    //remove all tuples which have this ccylinder as key
    auto t = tuple<CCylinder*, short>(cc, bindingSite);
    _possibleBindings.erase(t);

    //remove all tuples which have this as value

    //Iterate through the reverse map
       auto keys = _reversePossibleBindings[t];//keys that contain t as
       // value in possiblebindings
       for(auto k:keys){
               //get the iterator range that corresponds to this key.
               auto range = _possibleBindings.equal_range(k);
               //iterate through the range
               for(auto it = range.first; it != range.second;){
                       if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite) {
                            _possibleBindings.erase(it++);
                           }
                       else ++it;
               }
       }

    //remove from the reverse map.
    _reversePossibleBindings[t].clear();

    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();

    updateBindingReaction(oldN, newN);

    //remove all neighbors which have this binding site pair
    for (auto cn : _neighborLists[_nlIndex]->getNeighbors(cc->getCylinder())) {

        if(cn->getType() != complimentaryfID) continue;

        if(cn->getCompartment() != _compartment) {

            auto m = (LinkerBindingManager*)cn->getCompartment()->
            getFilamentBindingManagers()[_mIndex].get();

            if(find(affectedManagers.begin(), affectedManagers.end(), m) == affectedManagers.end())
            affectedManagers.push_back(m);
        }
    }

    //remove, update affected
    for(auto m : affectedManagers) {

        //Iterate through the reverse map
        auto keys = m->_reversePossibleBindings[t];//keys that contain t as
        // value in possiblebindings
        for(auto k:keys){
            //get the iterator range that corresponds to this key.
            auto range = m->_possibleBindings.equal_range(k);
            //iterate through the range
            for(auto it = range.first; it != range.second;){
                if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite) {
                    m->_possibleBindings.erase(it++);
                }
                else ++it;
            }
        }

        //remove from the reverse map.
        m->_reversePossibleBindings[t].clear();


        int oldNOther = m->_bindingSpecies->getN();
        int newNOther = m->numBindingSites();

        m->updateBindingReaction(oldNOther, newNOther);
    }
}

void LinkerBindingManager::removePossibleBindings(CCylinder* cc) {

    short _filamentType = cc->getType();
    if (_filamentType != _filamentIDvec[0] && _filamentType != _filamentIDvec[1]) return;

    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
        bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++){
#ifdef NLORIGINAL
        removePossibleBindings(cc, *bit);
#endif
    }
}

void LinkerBindingManager::updateAllPossibleBindings() {

    _possibleBindings.clear();
    _reversePossibleBindings.clear();
    int offset = 0;

    floatingpoint min1,min2,max1,max2;
    chrono::high_resolution_clock::time_point mins, mine, mins2, mine2,mints,minte;
    floatingpoint sqdisttermswithjustalpha;
    bool status1 = true;
    bool status2 = true;
    vector<floatingpoint> maxvec;
    vector<floatingpoint> minvec;
    int rejects16 = 0;
    int rejectsnavail =0;
    mints = chrono::high_resolution_clock::now();

    floatingpoint* cylsqmagnitudevector = SysParams::Mechanics().cylsqmagnitudevector;
    auto boundstate = SysParams::Mechanics().speciesboundvec;


    //lowest cylinder length fraction allowed for a binding site
    floatingpoint minparamcyl2;
    //highest cylinder length fraction allowed for a binding site
    floatingpoint maxparamcyl2;


    _possibleBindings.clear();
    mints = chrono::high_resolution_clock::now();

    for(auto c : _compartment->getCylinders()) {

      short _filamentType = c->getType();
      short complimentaryfID;
      if (_filamentType != _filamentIDvec[0] && _filamentType != _filamentIDvec[1]) return;
      else if (_filamentType == _filamentIDvec[0]) complimentaryfID = _filamentIDvec[1];
      else complimentaryfID = _filamentIDvec[0];

        auto x1 = c->getFirstBead()->vcoordinate();
        auto x2 = c->getSecondBead()->vcoordinate();
        auto cc = c->getCCylinder();
        vector<floatingpoint> X1X2 = {x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2]};
        //Loop through neighboring cylinders to c.
        for (auto cn : _neighborLists[_nlIndex]->getNeighbors(cc->getCylinder())) {

            if(cn->getParent() == c->getParent()) continue;
            if(cn->getType() != complimentaryfID) continue;
            if(c->getId() < cn->getId()) continue;
            auto ccn = cn->getCCylinder();

            auto x3 = cn->getFirstBead()->vcoordinate();
            auto x4 = cn->getSecondBead()->vcoordinate();

            vector<floatingpoint> X1X3 = {x3[0] - x1[0], x3[1] - x1[1], x3[2] - x1[2]};
            vector<floatingpoint> X3X4 = {x4[0] - x3[0], x4[1] - x3[1], x4[2] - x3[2]};
            floatingpoint maxdistsq = maxdistbetweencylinders(x1,x2,x3,x4);

            floatingpoint mindistsq = scalarprojection(X1X3, normalizeVector(vectorProduct(x1,x2,
                                                                                    x3,x4)));
            mindistsq = mindistsq * mindistsq;
            if(mindistsq > _rMaxsq || maxdistsq < _rMinsq) continue;

            floatingpoint X1X3squared = sqmagnitude(X1X3);
            floatingpoint X1X2squared = cylsqmagnitudevector[c->getStableIndex()];
            floatingpoint X1X3dotX1X2 = scalarprojection(X1X3, X1X2);
            floatingpoint X3X4squared = cylsqmagnitudevector[cn->getStableIndex()];
            floatingpoint X1X3dotX3X4 = scalarprojection(X1X3,X3X4);
            floatingpoint X3X4dotX1X2 = scalarprojection(X3X4, X1X2);
            mins2 = chrono::high_resolution_clock::now();

            int i = -1;
            dBInt = 2;
            //Loop through binding sites on cylinder c
            for(auto it1 = SysParams::Chemistry().bindingSites[_filamentType].begin();
                it1 != SysParams::Chemistry().bindingSites[_filamentType].end(); it1++) {
                i++;
                // DifBind
                if(dBInt % dBI != 0) {
                    dBInt += 1 ;
                    continue;
                } else {
                    dBInt = 1;}
                //now re add valid binding sites (check if binding site on cylinder c is empty)
                if (areEqual(boundstate[1][offset + SysParams::Chemistry()
                                           .bindingSites[_filamentType].size()
                                           *c->getStableIndex() + i], 1.0)) {

                    floatingpoint mp1, minparamcyl2, maxparamcyl2;
                    if (_filamentType == _filamentIDvec[0]){
                      mp1 = bindingsites1.at(i);
                      minparamcyl2 = bindingsites2.front();
                      maxparamcyl2 = bindingsites2.back();
                    }
                    else {
                      mp1 = bindingsites2.at(i);
                      minparamcyl2 = bindingsites1.front();
                      maxparamcyl2 = bindingsites1.back();
                    }

                    floatingpoint A = X3X4squared;
                    floatingpoint B = 2 * X1X3dotX3X4 - 2 * mp1 * X3X4dotX1X2;
                    floatingpoint C = X1X3squared + mp1 * mp1 * X1X2squared - 2 * mp1 *
                    X1X3dotX1X2;
                    floatingpoint C1 = C - _rMinsq;
                    floatingpoint C2 = C - _rMaxsq;
                    floatingpoint b2m4ac1 = B*B - 4*A*C1;
                    floatingpoint b2m4ac2 = B*B - 4*A*C2;
                    status1 = b2m4ac1 < 0;
                    status2 = b2m4ac2 < 0;
                    if(status1 && status2) continue;
                    maxvec.clear();
                    minvec.clear();
                    if(!status1){
                        min1 = (-B + sqrt(b2m4ac1))/(2*A);
                        min2 = (-B - sqrt(b2m4ac1))/(2*A);
                        if(min1<min2) {
                            minvec.push_back(min1);
                            minvec.push_back(min2);
                        }
                        else{
                            minvec.push_back(min2);
                            minvec.push_back(min1);
                        }
                        if(minvec.at(0)< minparamcyl2 && minvec.at(1) > maxparamcyl2) {
                            continue;
                        }
                    }
                    if(!status2){
                        max1 = (-B + sqrt(b2m4ac2))/(2*A);
                        max2 = (-B - sqrt(b2m4ac2))/(2*A);
                        if(max1<max2) {
                            maxvec.push_back(max1);
                            maxvec.push_back(max2);
                        }
                        else{
                            maxvec.push_back(max2);
                            maxvec.push_back(max1);
                        }
                        if(maxvec.at(0) > maxparamcyl2 || maxvec.at(1) <minparamcyl2){
                            continue;
                        }
                    }

                    dBInt = 2;
                    int j = -1;
                    //Loop through binding sites of complimentary cylinder
                    for(auto it2 = SysParams::Chemistry().bindingSites[complimentaryfID].begin();
                        it2 != SysParams::Chemistry().bindingSites[complimentaryfID].end(); it2++) {
                        j++;

                        // DifBind
                        if(dBInt % dBI != 0) {
                            dBInt += 1 ;
                            continue;
                        } else {
                            dBInt = 1;}

                        bool check2 = true;
                        //check if the binding site on cylinder nc is empty
                        if (areEqual(boundstate[1][offset + SysParams::Chemistry()
                                                   .bindingSites[complimentaryfID]
                                                   .size()*cn->getStableIndex() + j], 1.0)) {

                            //check distances..
                            floatingpoint mp2;
                            if (complimentaryfID == _filamentIDvec[0])
                              mp2 = bindingsites1.at(j);
                            else
                              mp2 = bindingsites2.at(j);

                            if(!status2) {
                                if (mp2 < maxvec.at(0) || mp2 > maxvec.at(1)) {
                                    continue;
                                }
                            }
                            if(!status1){
                                if (mp2 > minvec.at(0) && mp2 < minvec.at(1)) {
                                    continue;
                                }
                            }

                            if(check2) {
                                mins = chrono::high_resolution_clock::now();
                                auto t1 = tuple<CCylinder *, short>(cc, *it1);
                                auto t2 = tuple<CCylinder *, short>(ccn, *it2);

                                //add in correct order
                                if(c->getId() > cn->getId()) {
                                    _possibleBindings.emplace(t1, t2);
                                    _reversePossibleBindings[t2].push_back(t1);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();
    updateBindingReaction(oldN, newN);
}

/// Choose random binding sites based on current state
vector<tuple<CCylinder*, short>> LinkerBindingManager::chooseBindingSites() {

    assert((_possibleBindings.size() != 0)
           && "Major bug: Linker binding manager should not have zero binding \
           sites when called to choose a binding site.");

    int randomIndex = Rand::randInteger(0, _possibleBindings.size() - 1);
        auto it = _possibleBindings.begin();
    advance(it, randomIndex);
    return vector<tuple<CCylinder*, short>>{it->first, it->second};
}
void LinkerBindingManager::appendpossibleBindings(short boundInt, CCylinder* ccyl1,
                                                  CCylinder* ccyl2, short site1, short site2){
	short _filamentType = ccyl1->getType();
	short _nfilamentType = ccyl2->getType();

	bool cndn1 = _filamentType == _filamentIDvec[0] && _nfilamentType == _filamentIDvec[1];
	bool cndn2 = _filamentType == _filamentIDvec[1] && _nfilamentType == _filamentIDvec[0];

	if(!cndn1 && !cndn2 && boundInt != _boundInt) return;

	auto t1 = make_tuple(ccyl1, site1);
	auto t2 = make_tuple(ccyl2, site2);
	floatingpoint oldN=numBindingSites();
	_possibleBindings.emplace(t1,t2);
	floatingpoint newN=numBindingSites();
	updateBindingReaction(oldN,newN);
}

void LinkerBindingManager::printbindingsites() {
	cout<<"BINDINGSITES: CYL1(SIDX) CYL2(SIDX) SITE1 SITE2"<<endl;
	for(auto it1 = _possibleBindings.begin();it1!=_possibleBindings.end();
	it1++) {
		auto cyl1 = get<0>(it1->first)->getCylinder();
		auto bs1 = get<1>(it1->first);
		auto cyl2 = get<0>(it1->second)->getCylinder();
		auto bs2 = get<1>(it1->second);
		cout<<cyl1->getStableIndex()<<" "
		                      <<cyl2->getStableIndex()<<" "<<bs1<<" "<<bs2<<endl;
	}
}
#endif
//Deprecated. Used to compare stencil based search with original search.
bool LinkerBindingManager::isConsistent() {
    return true;
}
#if defined(NLSTENCILLIST) || defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
void LinkerBindingManager::addPossibleBindingsstencil(CCylinder* cc) {
	short _filamentType = cc->getType();
    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
        bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++) {
        addPossibleBindingsstencil(cc, *bit);
    }
}
void LinkerBindingManager::addPossibleBindingsstencil(CCylinder* cc, short bindingSite) {
#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
    auto HManager = _compartment->getHybridBindingSearchManager();
    HManager->addPossibleBindingsstencil(_Hbsmidvec,cc,bindingSite);
#else
	if (SysParams::INITIALIZEDSTATUS ) {
		short complimentaryfID;
		short _filamentType = cc->getType();
		if (_filamentType != _filamentIDvec[0] && _filamentType != _filamentIDvec[1])
			return;
		else if (_filamentType == _filamentIDvec[0]) complimentaryfID = _filamentIDvec[1];
		else complimentaryfID = _filamentIDvec[0];

		//if we change other managers copy number
		vector<LinkerBindingManager *> affectedManagers;

		//add valid binding sites
		if (areEqual(cc->getCMonomer(bindingSite)->speciesBound(
				SysParams::Chemistry().linkerBoundIndex[_filamentType])->getN(),
		             (floatingpoint) 1.0)) {

			//loop through neighbors
			//now re add valid based on CCNL
			vector<Cylinder *> Neighbors;
//#ifdef HYBRID_NLSTENCILLIST
//        Neighbors = cc->getCompartment()->getHybridBindingSearchManager()->getHNeighbors
//                (cc->getCylinder(),HNLID);
//#else

            Neighbors = _neighborLists[_nlIndex]->getNeighborsstencil(cc->getCylinder());
//#endif
            for (auto cn : Neighbors) {

                Cylinder *c = cc->getCylinder();

                if (cn->getParent() == c->getParent()) continue;
                if (cn->getType() != complimentaryfID) continue;

                auto ccn = cn->getCCylinder();

                for (auto it = SysParams::Chemistry().bindingSites[complimentaryfID].begin();
                        it !=
                        SysParams::Chemistry().bindingSites[complimentaryfID].end(); it++) {

                    if (areEqual(ccn->getCMonomer(*it)->speciesBound(
                            SysParams::Chemistry().linkerBoundIndex[complimentaryfID])->getN(),
                                    (floatingpoint) 1.0)) {

                        //check distances..
                        auto mp1 = (float) bindingSite /
                                    SysParams::Geometry().cylinderNumMon[_filamentType];
                        auto mp2 = (float) *it /
                                    SysParams::Geometry().cylinderNumMon[complimentaryfID];

                        auto x1 = c->getFirstBead()->vcoordinate();
                        auto x2 = c->getSecondBead()->vcoordinate();
                        auto x3 = cn->getFirstBead()->vcoordinate();
                        auto x4 = cn->getSecondBead()->vcoordinate();

                        auto m1 = midPointCoordinate(x1, x2, mp1);
                        auto m2 = midPointCoordinate(x3, x4, mp2);

                        floatingpoint distsq = twoPointDistancesquared(m1, m2);

                        if (distsq > _rMaxsq || distsq < _rMinsq) continue;

                        auto t1 = tuple<CCylinder *, short>(cc, bindingSite);
                        auto t2 = tuple<CCylinder *, short>(ccn, *it);

                        //add in correct order
                        if (c->getId() > cn->getId())
                            _possibleBindingsstencil.emplace(t1, t2);
                        else {
                            //add in this compartment
                            if (cn->getCompartment() == _compartment) {

                                _possibleBindingsstencil.emplace(t2, t1);
                            }
                                //add in other
                            else {
                                auto m = (LinkerBindingManager *) cn->getCompartment()->
                                        getFilamentBindingManagers()[_mIndex].get();

                                affectedManagers.push_back(m);

                                m->_possibleBindingsstencil.emplace(t2, t1);
                            }
                        }
                    }
                }
            }
        }

        //update affected
        for (auto m : affectedManagers) {

            int oldNOther = m->_bindingSpecies->getN();
            int newNOther = m->numBindingSitesstencil();

            m->updateBindingReaction(oldNOther, newNOther);
        }

        //update this manager
        int oldN = _bindingSpecies->getN();
        int newN = numBindingSitesstencil();

        updateBindingReaction(oldN, newN);
    }
#endif
}
void LinkerBindingManager::updateAllPossibleBindingsstencil() {
#ifdef NLSTENCILLIST
    _possibleBindingsstencil.clear();
    floatingpoint min1,min2,max1,max2;
    bool status1 = true;
    bool status2 = true;
    floatingpoint minveca[2];
    floatingpoint maxveca[2];

    int maxnbs = SysParams::Chemistry().maxbindingsitespercylinder;
    floatingpoint* cylsqmagnitudevector = SysParams::Mechanics().cylsqmagnitudevector;
    auto boundstate = SysParams::Mechanics().speciesboundvec;

    const auto& cylinderInfoData = Cylinder::getDbData();

    int Ncylincmp =  _compartment->getCylinders().size();
    int* cindexvec = new int[Ncylincmp]; //stores cindex of cylinders in this compartment
    vector<vector<int>> ncindices; //cindices of cylinders in neighbor list.
    vector<int>ncindex; //helper vector
    long id = 0;
    for(auto c : _compartment->getCylinders()){
        cindexvec[id] = c->getStableIndex();
        id++;
        for (auto cn : _neighborLists[_nlIndex]->getNeighborsstencil(c)) {
            if(c->getId() > cn->getId())
                ncindex.push_back(cn->getStableIndex());
        }
        ncindices.push_back(ncindex);
        ncindex.clear();
    }

    for(int i=0;i<Ncylincmp;i++){
        int cindex = cindexvec[i];
        const auto& c = cylinderInfoData[cindex];

        int nbs1, nbs2;
        short complimentaryfID;

        short _filamentType = c.type;
        if (_filamentType != _filamentIDvec[0] && _filamentType != _filamentIDvec[1])
            return;
        else if (_filamentType == _filamentIDvec[0])
            complimentaryfID = _filamentIDvec[1];
        else
            complimentaryfID = _filamentIDvec[0];

        nbs1 = SysParams::Chemistry().bindingSites[_filamentType].size();
        nbs2 = SysParams::Chemistry().bindingSites[complimentaryfID].size();

        const auto& x1 = Bead::getStableElement(c.beadIndices[0])->coord;
        const auto& x2 = Bead::getStableElement(c.beadIndices[1])->coord;
        floatingpoint X1X2[3] ={x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2]};

        int* cnindices = ncindices[i].data();
        for(int arraycount = 0; arraycount < ncindices[i].size();arraycount++){
            int cnindex = cnindices[arraycount];
            const auto& cn = cylinderInfoData[cnindex];
            // if(c.ID < cn.ID) {counter++; continue;} commented as the above vector does
            // not contain ncs that will fail this cndn.
            if(c.filamentId == cn.filamentId){
                continue;}
            if(cn.type != complimentaryfID){
                continue;}

            const auto& x3 = Bead::getStableElement(cn.beadIndices[0])->coord;
            const auto& x4 = Bead::getStableElement(cn.beadIndices[1])->coord;
            floatingpoint X1X3[3] = {x3[0] - x1[0], x3[1] - x1[1], x3[2] - x1[2]};
            floatingpoint X3X4[3] = {x4[0] - x3[0], x4[1] - x3[1], x4[2] - x3[2]};
            floatingpoint X1X3squared = sqmagnitude(X1X3);
            floatingpoint X1X2squared = cylsqmagnitudevector[cindex];
            floatingpoint X1X3dotX1X2 = scalarprojection(X1X3, X1X2);
            floatingpoint X3X4squared = cylsqmagnitudevector[cnindex];
            floatingpoint X1X3dotX3X4 = scalarprojection(X1X3,X3X4);
            floatingpoint X3X4dotX1X2 = scalarprojection(X3X4, X1X2);
            for(int posfil =0; posfil<nbs1;posfil++) {
                //now re add valid binding sites
                if (areEqual(boundstate[1][maxnbs * cindex + posfil], (floatingpoint)1.0)) {

                    floatingpoint mp1, minparamcyl2, maxparamcyl2;
                    if (_filamentType == _filamentIDvec[0]){
                        mp1 = bindingsites1.at(posfil);
                        minparamcyl2 = bindingsites2.front();
                        maxparamcyl2 = bindingsites2.back();
                    }
                    else {
                        mp1 = bindingsites2.at(posfil);
                        minparamcyl2 = bindingsites1.front();
                        maxparamcyl2 = bindingsites1.back();
                    }

                    floatingpoint A = X3X4squared;
                    floatingpoint B = 2 * X1X3dotX3X4 - 2 * mp1 * X3X4dotX1X2;
                    floatingpoint C = X1X3squared + mp1 * mp1 * X1X2squared - 2 * mp1 * X1X3dotX1X2;
                    floatingpoint C1 = C - _rMinsq;
                    floatingpoint C2 = C - _rMaxsq;
                    floatingpoint b2m4ac1 = B*B - 4*A*C1;
                    floatingpoint b2m4ac2 = B*B - 4*A*C2;
                    status1 = b2m4ac1 < 0;
                    status2 = b2m4ac2 < 0;
                    if(status1 && status2) {
                        continue;}
                    if(!status1){
                        min1 = (-B + sqrt(b2m4ac1))/(2*A);
                        min2 = (-B - sqrt(b2m4ac1))/(2*A);
                        if(min1<min2) {
                            minveca[0] = (min1);
                            minveca[1] = (min2);
                        }
                        else{
                            minveca[0] = (min2);
                            minveca[1] = (min1);
                        }
                        if(minveca[0]< minparamcyl2 && minveca[1] > maxparamcyl2) {
                            continue;
                        }
                    }
                    if(!status2){
                        max1 = (-B + sqrt(b2m4ac2))/(2*A);
                        max2 = (-B - sqrt(b2m4ac2))/(2*A);
                        if(max1<max2) {
                            maxveca[0] = (max1);
                            maxveca[1] = (max2);
                        }
                        else{
                            maxveca[0] = (max2);
                            maxveca[1] = (max1);
                        }
                        if(maxveca[0] > maxparamcyl2 || maxveca[1] <minparamcyl2){
                            continue;
                        }
                    }
                    for(int posCfil = 0; posCfil<nbs2;posCfil++){
                        if (areEqual(boundstate[1][maxnbs * cnindex + posCfil], (floatingpoint)1.0)) {

                            //check distances..
                            floatingpoint mp2;
                            if (complimentaryfID == _filamentIDvec[0])
                                mp2 = bindingsites1.at(posCfil);
                            else
                                mp2 = bindingsites2.at(posCfil);
                            if(!status2) {
                                if (mp2 < maxveca[0] || mp2 > maxveca[1])
                                    continue;
                            }
                            if(!status1){
                                if (mp2 > minveca[0] && mp2 < minveca[1])
                                    continue;
                            }

                            auto it1 = SysParams::Chemistry()
                                    .bindingSites[_filamentType][posfil];
                            auto it2 = SysParams::Chemistry().bindingSites[complimentaryfID][posCfil];

                            auto t1 = tuple<CCylinder *, short>(cylinderInfoData[cindex].chemCylinder, it1);
                            auto t2 = tuple<CCylinder *, short>(cylinderInfoData[cnindex].chemCylinder, it2);
                            //add in correct order
                            _possibleBindingsstencil.emplace(t1, t2);
                        }
                    }
                }
            }
        }
    }
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSitesstencil();
    updateBindingReaction(oldN, newN);
    delete[] cindexvec;
#else
    cout<<"Erroneous function call. Please check compiler macros. Exiting."<<endl;
    exit(EXIT_FAILURE);
#endif
}
void LinkerBindingManager::appendpossibleBindingsstencil(short boundInt, CCylinder* ccyl1,
                                                        CCylinder* ccyl2, short site1, short site2){

	short _filamentType = ccyl1->getType();
	short _nfilamentType = ccyl2->getType();

	bool cndn1 = _filamentType == _filamentIDvec[0] && _nfilamentType == _filamentIDvec[1];
	bool cndn2 = _filamentType == _filamentIDvec[1] && _nfilamentType == _filamentIDvec[0];

	if(!cndn1 && !cndn2 && boundInt != _boundInt) return;

	#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
	auto HManager = _compartment->getHybridBindingSearchManager();
	HManager->appendPossibleBindingsstencil(_Hbsmidvec, ccyl1, ccyl2, site1, site2);
	#else
	auto t1 = make_tuple(ccyl1, site1);
	auto t2 = make_tuple(ccyl2, site2);
	floatingpoint oldN=numBindingSitesstencil();
	_possibleBindingsstencil.emplace(t1,t2);
	floatingpoint newN=numBindingSitesstencil();
	updateBindingReaction(oldN,newN);
	#endif
}
void LinkerBindingManager::removePossibleBindingsstencil(CCylinder* cc) {

	short _filamentType = cc->getType();

    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
        bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++)
        removePossibleBindingsstencil(cc, *bit);
}
void LinkerBindingManager::removePossibleBindingsstencil(CCylinder* cc, short bindingSite) {
#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
    auto HManager = _compartment->getHybridBindingSearchManager();
    HManager->removePossibleBindingsstencil(_Hbsmidvec, cc, bindingSite);
/*    for(auto C:SubSystem::getstaticgrid()->getCompartments()){
        C->getHybridBindingSearchManager()->checkoccupancySIMD(_Hbsmidvec);
    }*/

#else

	short _filamentType = cc->getType();
    short complimentaryfID;
    if (_filamentType != _filamentIDvec[0] && _filamentType != _filamentIDvec[1]) return;
    else if (_filamentType == _filamentIDvec[0]) complimentaryfID = _filamentIDvec[1];
    else complimentaryfID = _filamentIDvec[0];

    //if we change other managers copy number
    vector<LinkerBindingManager*> affectedManagers;

    //remove all tuples which have this ccylinder as key
    auto t = tuple<CCylinder*, short>(cc, bindingSite);
    _possibleBindingsstencil.erase(t);

    //remove all tuples which have this as value
    for (auto it = _possibleBindingsstencil.begin(); it != _possibleBindingsstencil.end();
         ) {

        if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite)
            _possibleBindingsstencil.erase(it++);

        else ++it;
    }

    int oldN = _bindingSpecies->getN();
    int newN = numBindingSitesstencil();

    updateBindingReaction(oldN, newN);

    //remove all neighbors which have this binding site pair
    for (auto cn : _neighborLists[_nlIndex]->getNeighborsstencil(cc->getCylinder())) {

        if(cn->getType() != complimentaryfID) continue;

        if(cn->getCompartment() != _compartment) {

            auto m = (LinkerBindingManager*)cn->getCompartment()->
                    getFilamentBindingManagers()[_mIndex].get();

            if(find(affectedManagers.begin(), affectedManagers.end(), m) == affectedManagers.end())
                affectedManagers.push_back(m);
        }
    }
//#endif

    //remove, update affected
    for(auto m : affectedManagers) {

        for (auto it = m->_possibleBindingsstencil.begin(); it !=
             m->_possibleBindingsstencil.end(); ) {

            if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite)
            m->_possibleBindingsstencil.erase(it++);

            else ++it;
        }

        int oldNOther = m->_bindingSpecies->getN();
        int newNOther = m->numBindingSitesstencil();

        m->updateBindingReaction(oldNOther, newNOther);
    }
#endif
}
void LinkerBindingManager::crosscheck(){
    cout<<"Linker NLORIGINAL size "<<_possibleBindings.size()<<" NLSTENCIL size "
        <<_possibleBindingsstencil.size()<<endl;
    if(_possibleBindings.size() != _possibleBindingsstencil.size())
        cout<<"Linker.. The two methods compared do not yield the same number of "
                "binding sites"<<endl;
    short matches = 0;

    for(auto it1 = _possibleBindings.begin();it1!=_possibleBindings.end();it1++){
        short matchperelement = 0;

        auto cyl1_o = get<0>(it1->first)->getCylinder();
        auto bs1_o =  get<1>(it1->first);
        auto cyl2_o = get<0>(it1->second)->getCylinder();
        auto bs2_o =  get<1>(it1->second);
        //
        short sum = 0;
        for(auto it2 = _possibleBindingsstencil.begin();it2!=_possibleBindingsstencil.end();
            it2++){
            auto cyl1_s = get<0>(it2->first)->getCylinder();
            auto bs1_s =  get<1>(it2->first);
            auto cyl2_s = get<0>(it2->second)->getCylinder();
            auto bs2_s =  get<1>(it2->second);
            sum = 0;
            if(cyl1_o->getId() == cyl1_s->getId() )
                sum++;
            if(cyl2_o->getId() == cyl2_s->getId() )
                sum++;
            if(bs1_o == bs1_s )
                sum++;
            if(bs2_o == bs2_s )
                sum++;
            if(sum == 4) {
//                cout << "match found" << endl;
                if(matchperelement == 0 ) {
                    matches++;
                    matchperelement++;
                } else {
                    cout << "ERROR. Multiple matches for chosen binding site pair in "
                    "STENCILLIST. Check stencillist.." << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
    std::cout<<"Linker possible bindings size "<<_possibleBindings.size()<<" Total matches"
            " "<<
             matches<<endl;
        if(_possibleBindings.size() != matches || _possibleBindings.size() !=
                                                      _possibleBindingsstencil.size()){
        cout<<"Linker. All binding site pairs are not found in Stencil list"<<endl;
        exit(EXIT_FAILURE);
    }
}
vector<tuple<CCylinder*, short>> LinkerBindingManager::chooseBindingSitesstencil() {
#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
    auto HManager = _compartment->getHybridBindingSearchManager();
    return HManager->chooseBindingSitesstencil(_Hbsmidvec);
#else
    assert((_possibleBindingsstencil.size() != 0)
           && "Major bug: Linker binding manager should not have zero binding \
                   sites when called to choose a binding site.");

    int randomIndex = Rand::randInteger(0, _possibleBindingsstencil.size() - 1);
    auto it = _possibleBindingsstencil.begin();

    advance(it, randomIndex);

    return vector<tuple<CCylinder*, short>>{it->first, it->second};
#endif
}
void LinkerBindingManager::clearpossibleBindingsstencil() {
    #ifdef NLSTENCILLIST
    floatingpoint oldN=numBindingSitesstencil();
        _possibleBindingsstencil.clear();
        updateBindingReaction(oldN,0);
    #else
    auto HManager = _compartment->getHybridBindingSearchManager();
    HManager->clearPossibleBindingsstencil(_Hbsmidvec);
    #endif
}
int LinkerBindingManager::numBindingSitesstencil() {
#ifdef NLSTENCILLSIT
    return _possibleBindingsstencil.size();
#else
    auto HManager = _compartment->getHybridBindingSearchManager();
    return HManager->numBindingSitesstencil(_Hbsmidvec);
#endif

}

void LinkerBindingManager::printbindingsitesstencil() {
#ifdef NLSTENCILLSIT
	cout<<"BINDINGSITES: CYL1(SIDX) CYL2(SIDX) SITE1 SITE2"<<endl;
	for(auto it1 = _possibleBindingsstencil.begin();it1!=_possibleBindingsstencil.end();
	it1++) {
		auto cyl1 = get<0>(it1->first)->getCylinder();
		auto bs1 = get<1>(it1->first);
		auto cyl2 = get<0>(it1->second)->getCylinder();
		auto bs2 = get<1>(it1->second);
		cout<<cyl1->getStableIndex()<<" "
		                      <<cyl2->getStableIndex()<<" "<<bs1<<" "<<bs2<<endl;
	}
#else
	auto HManager = _compartment->getHybridBindingSearchManager();
	HManager->printbindingsitesstencil(_Hbsmidvec);
#endif
}
#endif
#ifdef CUDAACCL_NL
void LinkerBindingManager::assigncudavars() {
    //    if(gpu_rminmax == NULL) {
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_rminmax, 2 * sizeof(floatingpoint)), "cuda data transfer", " "
                            "BindingManager.cu");
    floatingpoint dist[2];
    dist[0] = _rMin;
    dist[1] = _rMax;
    CUDAcommon::handleerror(cudaMemcpy(gpu_rminmax, dist, 2 *sizeof(floatingpoint), cudaMemcpyHostToDevice));
    //    }
    //    if(gpu_numpairs == NULL) {
    CUDAcommon::handleerror(cudaMalloc((void **) &gpu_numpairs, sizeof(int)), "cuda data transfer", " "
                            "BindingManager.cu");
    //    int n[1];
    //    n[0] = 0;
    //    CUDAcommon::handleerror(cudaMemcpy(gpu_numpairs, n, sizeof(int), cudaMemcpyHostToDevice));
    //    }
    //    delete dist;
}
void LinkerBindingManager::freecudavars() {
    CUDAcommon::handleerror(cudaFree(gpu_rminmax),"cudaFree", "BindingManager");
    CUDAcommon::handleerror(cudaFree(gpu_numpairs),"cudaFree", "BindingManager");
}
#endif

//MOTOR
MotorBindingManager::MotorBindingManager(
    ReactionBase* reaction,
    Compartment* compartment,
    short linkerType,
    vector<short> filamentIDvec,
    int motorSpeciesIndex1,
    int motorSpeciesIndex2,
    float rMax, float rMin)

: FilamentBindingManager(reaction, compartment, linkerType, filamentIDvec),
    _rMin(rMin), _rMax(rMax),
    motorSpeciesIndices_ { motorSpeciesIndex1, motorSpeciesIndex2 }
{

    _rMinsq =_rMin * _rMin;
    _rMaxsq = _rMax * _rMax;

    //find the pair binding species
    RSpecies** rs = reaction->rspecies();

    string name = rs[ML_RXN_INDEX]->getSpecies().getName();

    _bindingSpecies = _compartment->findSpeciesByName(name);

    //initialize ID's based on number of species in compartment
    //    int numSpecies = rs[ML_RXN_INDEX + 1]->getSpecies().getN();

    //DEPRECATED AS OF 9/22/16
    //    for(int i = 0; i < numSpecies; i++)
    //        _unboundIDs.push_back(MotorGhost::_motorGhosts.getId());


    //attach an rspecies callback to this species
    Species* sd = &(rs[ML_RXN_INDEX + 1]->getSpecies());

    UpdateMotorIDCallback mcallback(linkerType);
    sd->connect(mcallback);
    _rMaxsq = rMax*rMax;
    _rMinsq = rMin*rMin;

    for(auto it1 = SysParams::Chemistry().bindingSites[filamentIDvec[0]].begin();
        it1 != SysParams::Chemistry().bindingSites[filamentIDvec[0]].end(); it1++) {
        bindingsites1.push_back((float)*it1 / SysParams::Geometry()
                .cylinderNumMon[filamentIDvec[0]]);
    }
	for(auto it1 = SysParams::Chemistry().bindingSites[filamentIDvec[1]].begin();
	    it1 != SysParams::Chemistry().bindingSites[filamentIDvec[1]].end(); it1++) {
		bindingsites2.push_back((float)*it1 / SysParams::Geometry()
				.cylinderNumMon[filamentIDvec[1]]);
	}
}

#ifdef NLORIGINAL
void MotorBindingManager::addPossibleBindings(CCylinder* cc, short bindingSite) {

    if (SysParams::INITIALIZEDSTATUS ) {

        short complimentaryfID;
        short _filamentType = cc->getType();
        if (_filamentType != _filamentIDvec[0] && _filamentType != _filamentIDvec[1])
            return;
        else if (_filamentType == _filamentIDvec[0]) complimentaryfID = _filamentIDvec[1];
        else complimentaryfID = _filamentIDvec[0];

        //if we change other managers copy number
        vector<MotorBindingManager *> affectedManagers;

        //add valid binding sites
        if (areEqual(cc->getCMonomer(bindingSite)->speciesBound(
                SysParams::Chemistry().motorBoundIndex[_filamentType])->getN(),
                        (floatingpoint) 1.0)) {

            //loop through neighbors
            //now re add valid based on CCNL
            vector<Cylinder *> nList = _neighborLists[_nlIndex]->getNeighbors(
                    cc->getCylinder());

            //loop through neighbors
            //now re add valid based on CCNL
            for (auto cn : nList) {

                Cylinder *c = cc->getCylinder();

                if (cn->getParent() == c->getParent()) continue;

                short _nfilamentType = cn->getType();
                if (_nfilamentType != complimentaryfID) return;
                if (c->getId() < cn->getId()) continue;

                auto ccn = cn->getCCylinder();

                for (auto it = SysParams::Chemistry().bindingSites[complimentaryfID].begin();
                        it !=
                        SysParams::Chemistry().bindingSites[complimentaryfID].end(); it++) {

                    if (areEqual(ccn->getCMonomer(*it)->speciesBound(
                            SysParams::Chemistry().motorBoundIndex[complimentaryfID])->getN(),
                                    (floatingpoint) 1.0)) {

                        //check distances..
                        auto mp1 = (float) bindingSite /
                                    SysParams::Geometry().cylinderNumMon[_filamentType];
                        auto mp2 = (float) *it /
                                    SysParams::Geometry().cylinderNumMon[complimentaryfID];

                        auto x1 = c->getFirstBead()->vcoordinate();
                        auto x2 = c->getSecondBead()->vcoordinate();
                        auto x3 = cn->getFirstBead()->vcoordinate();
                        auto x4 = cn->getSecondBead()->vcoordinate();

                        auto m1 = midPointCoordinate(x1, x2, mp1);
                        auto m2 = midPointCoordinate(x3, x4, mp2);

                        floatingpoint distSq = twoPointDistancesquared(m1, m2);

                        if (distSq > _rMaxsq || distSq < _rMinsq) continue;

                        auto t1 = tuple<CCylinder *, short>(cc, bindingSite);
                        auto t2 = tuple<CCylinder *, short>(ccn, *it);

                        //add in correct order
                        if (c->getId() > cn->getId()) {
                            _possibleBindings.emplace(t1, t2);
                            _reversePossibleBindings[t2].push_back(t1);

                        } else {
                            //add in this compartment
                            if (cn->getCompartment() == _compartment) {
                                _possibleBindings.emplace(t1, t2);
                                _reversePossibleBindings[t2].push_back(t1);

                            }
                                //add in other
                            else {

                                auto m = (MotorBindingManager *) cn->getCompartment()->
                                        getFilamentBindingManagers()[_mIndex].get();

                                affectedManagers.push_back(m);

                                m->_possibleBindings.emplace(t2, t1);
                                m->_reversePossibleBindings[t1].push_back(t2);

                            }
                        }
                    }
                }
            }
        }

        //update affected
        for (auto m : affectedManagers) {

            int oldNOther = m->_bindingSpecies->getN();
            int newNOther = m->numBindingSites();

            m->updateBindingReaction(oldNOther, newNOther);
        }

        //update this manager
        int oldN = _bindingSpecies->getN();
        int newN = numBindingSites();

        updateBindingReaction(oldN, newN);
    }
}

void MotorBindingManager::addPossibleBindings(CCylinder* cc) {

    short _filamentType = cc->getType();

    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
        bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++){
#ifdef NLORIGINAL
        addPossibleBindings(cc, *bit);
#endif
#if defined(NLSTENCILLIST) || defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
        addPossibleBindingsstencil(cc, *bit);
#endif
    }
}


void MotorBindingManager::removePossibleBindings(CCylinder* cc, short bindingSite) {

    short complimentaryfID;
    short _filamentType = cc->getType();
    if (_filamentType != _filamentIDvec[0] && _filamentType != _filamentIDvec[1]) return;
    else if (_filamentType == _filamentIDvec[0]) complimentaryfID = _filamentIDvec[1];
    else complimentaryfID = _filamentIDvec[0];

    //if we change other managers copy number
    vector<MotorBindingManager*> affectedManagers;

    //remove all tuples which have this ccylinder as key
    auto t = tuple<CCylinder*, short>(cc, bindingSite);
    _possibleBindings.erase(t);

    //remove all tuples which have this as value
    //Iterate through the reverse map
    auto keys = _reversePossibleBindings[t];//keys that contain t as
    // value in possiblebindings
    for(auto k:keys){
        //get the iterator range that corresponds to this key.
        auto range = _possibleBindings.equal_range(k);
        //iterate through the range
        for(auto it = range.first; it != range.second;){
            if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite) {
                _possibleBindings.erase(it++);
            }
            else ++it;
        }
    }

    //remove from the reverse map.
    _reversePossibleBindings[t].clear();


    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();

    updateBindingReaction(oldN, newN);

    //remove all neighbors which have this binding site pair
    for (auto cn : _neighborLists[_nlIndex]->getNeighbors(cc->getCylinder())) {

        if(cn->getType() != complimentaryfID) continue;

        if(cn->getCompartment() != _compartment) {

            auto m = (MotorBindingManager*)cn->getCompartment()->
            getFilamentBindingManagers()[_mIndex].get();

            if(find(affectedManagers.begin(), affectedManagers.end(), m) == affectedManagers.end())
            affectedManagers.push_back(m);
        }
    }

    //remove, update affected
    for(auto m : affectedManagers) {

        //Iterate through the reverse map
        auto keys = m->_reversePossibleBindings[t];//keys that contain t as
        // value in possiblebindings
        for(auto k:keys){
            //get the iterator range that corresponds to this key.
            auto range = m->_possibleBindings.equal_range(k);
            //iterate through the range
            for(auto it = range.first; it != range.second;){
                if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite) {
                    m->_possibleBindings.erase(it++);
                }
                else ++it;
            }

        }

        //remove from the reverse map.
        m->_reversePossibleBindings[t].clear();


        int oldNOther = m->_bindingSpecies->getN();
        int newNOther = m->numBindingSites();

        m->updateBindingReaction(oldNOther, newNOther);
    }
}

void MotorBindingManager::removePossibleBindings(CCylinder* cc) {

    short _filamentType = cc->getType();

    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
        bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++){
#ifdef NLORIGINAL
        removePossibleBindings(cc, *bit);
#endif
    }
}

void MotorBindingManager::updateAllPossibleBindings() {

    _possibleBindings.clear();
    _reversePossibleBindings.clear();

    int offset = 0;

    floatingpoint min1,min2,max1,max2;
    chrono::high_resolution_clock::time_point mins, mine, mins2, mine2,mints,minte;

    floatingpoint sqdisttermswithjustalpha;
    bool status1 = true;
    bool status2 = true;
    vector<floatingpoint> maxvec;
    vector<floatingpoint> minvec;

    mints = chrono::high_resolution_clock::now();
    vector<floatingpoint> bindingsites;
    floatingpoint* cylsqmagnitudevector = SysParams::Mechanics().cylsqmagnitudevector;
    auto boundstate = SysParams::Mechanics().speciesboundvec;
    _possibleBindings.clear();
    mints = chrono::high_resolution_clock::now();
    for(auto c : _compartment->getCylinders()) {

        short _filamentType = c->getType();
        short complimentaryfID;
        if (_filamentType != _filamentIDvec[0] && _filamentType != _filamentIDvec[1])
        	return;
        else if (_filamentType == _filamentIDvec[0]) complimentaryfID = _filamentIDvec[1];
        else complimentaryfID = _filamentIDvec[0];

        auto x1 = c->getFirstBead()->vcoordinate();
        auto x2 = c->getSecondBead()->vcoordinate();
        auto cc = c->getCCylinder();
        vector<floatingpoint> X1X2 = {x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2]};

        for (auto cn : _neighborLists[_nlIndex]->getNeighbors(cc->getCylinder())) {

            if(cn->getParent() == c->getParent()) continue;
            if(cn->getType() != complimentaryfID) continue;
            if(c->getId() < cn->getId()) continue;
            auto ccn = cn->getCCylinder();
            auto x3 = cn->getFirstBead()->vcoordinate();
            auto x4 = cn->getSecondBead()->vcoordinate();

            vector<floatingpoint> X1X3 = {x3[0] - x1[0], x3[1] - x1[1], x3[2] - x1[2]};
            vector<floatingpoint> X3X4 = {x4[0] - x3[0], x4[1] - x3[1], x4[2] - x3[2]};
            floatingpoint maxdistsq = maxdistbetweencylinders(x1,x2,x3,x4);

            floatingpoint mindistsq = scalarprojection(X1X3, normalizeVector(vectorProduct(x1,x2,
                                                                                    x3,x4)));
            mindistsq = mindistsq * mindistsq;
            if(mindistsq > _rMaxsq || maxdistsq < _rMinsq) continue;

            floatingpoint X1X3squared = sqmagnitude(X1X3);
            floatingpoint X1X2squared = cylsqmagnitudevector[c->getStableIndex()];
            floatingpoint X1X3dotX1X2 = scalarprojection(X1X3, X1X2);
            floatingpoint X3X4squared = cylsqmagnitudevector[cn->getStableIndex()];
            floatingpoint X1X3dotX3X4 = scalarprojection(X1X3,X3X4);
            floatingpoint X3X4dotX1X2 = scalarprojection(X3X4, X1X2);
            mins2 = chrono::high_resolution_clock::now();
            int i = -1;
            //For every binding site in cylinder c
            for(auto it1 = SysParams::Chemistry().bindingSites[_filamentType].begin();
                it1 != SysParams::Chemistry().bindingSites[_filamentType].end(); it1++) {
                i++;
                //Check if binding site is empty
                if (areEqual(boundstate[2][offset + SysParams::Chemistry()
                                           .bindingSites[_filamentType].size()
                                           *c->getStableIndex() + i], 1.0)) {
                    floatingpoint mp1, minparamcyl2, maxparamcyl2;
                    if (_filamentType == _filamentIDvec[0]){
                        mp1 = bindingsites1.at(i);
                        minparamcyl2 = bindingsites2.front();
                        maxparamcyl2 = bindingsites2.back();
                    }
                    else {
                        mp1 = bindingsites2.at(i);
                        minparamcyl2 = bindingsites1.front();
                        maxparamcyl2 = bindingsites1.back();
                    }
                    floatingpoint A = X3X4squared;
                    floatingpoint B = 2 * X1X3dotX3X4 - 2 * mp1 * X3X4dotX1X2;
                    floatingpoint C = X1X3squared + mp1 * mp1 * X1X2squared - 2 * mp1 *
                    X1X3dotX1X2;
                    floatingpoint C1 = C - _rMinsq;
                    floatingpoint C2 = C - _rMaxsq;
                    floatingpoint b2m4ac1 = B*B - 4*A*C1;
                    floatingpoint b2m4ac2 = B*B - 4*A*C2;
                    status1 = b2m4ac1 < 0;
                    status2 = b2m4ac2 < 0;
                    if(status1 && status2) continue;
                    maxvec.clear();
                    minvec.clear();
                    if(!status1){
                        min1 = (-B + sqrt(b2m4ac1))/(2*A);
                        min2 = (-B - sqrt(b2m4ac1))/(2*A);
                        if(min1<min2) {
                            minvec.push_back(min1);
                            minvec.push_back(min2);
                        }
                        else{
                            minvec.push_back(min2);
                            minvec.push_back(min1);
                        }
                        if(minvec.at(0)< minparamcyl2 && minvec.at(1) > maxparamcyl2) {
                            continue;
                        }
                    }
                    if(!status2){
                        max1 = (-B + sqrt(b2m4ac2))/(2*A);
                        max2 = (-B - sqrt(b2m4ac2))/(2*A);
                        if(max1<max2) {
                            maxvec.push_back(max1);
                            maxvec.push_back(max2);
                        }
                        else{
                            maxvec.push_back(max2);
                            maxvec.push_back(max1);
                        }
                        if(maxvec.at(0) > maxparamcyl2 || maxvec.at(1) <minparamcyl2){
                            continue;
                        }
                    }
                    int j =-1;
                    for(auto it2 = SysParams::Chemistry().bindingSites[complimentaryfID].begin();
                        it2 != SysParams::Chemistry().bindingSites[complimentaryfID].end(); it2++) {
                        j++;
                        bool check2 = true;
                        if (areEqual(boundstate[2][offset + SysParams::Chemistry()
                                                   .bindingSites[complimentaryfID]
                                                   .size()*cn->getStableIndex() + j],
                                                   		(floatingpoint)1.0)) {
                            //total++;
                            //check distances..
                            floatingpoint mp2;
                            if (complimentaryfID == _filamentIDvec[0])
                                mp2 = bindingsites1.at(j);
                            else
                                mp2 = bindingsites2.at(j);

                            if(!status2) {
                                if (mp2 < maxvec.at(0) || mp2 > maxvec.at(1))
                                    continue;
                            }
                            if(!status1){
                                if (mp2 > minvec.at(0) && mp2 < minvec.at(1))
                                    continue;
                            }

                            auto t1 = tuple<CCylinder*, short>(cc, *it1);
                            auto t2 = tuple<CCylinder*, short>(ccn, *it2);

                            //add in correct order
                            if(c->getId() > cn->getId()) {
                                _possibleBindings.emplace(t1, t2);
                                _reversePossibleBindings[t2].push_back(t1);
                            }
                        }
                    }
                }
            }
        }
    }
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSites();
    updateBindingReaction(oldN, newN);
}

/// Choose random binding sites based on current state
vector<tuple<CCylinder*, short>> MotorBindingManager::chooseBindingSites() {

    assert((_possibleBindings.size() != 0)
           && "Major bug: Motor binding manager should not have zero binding \
                   sites when called to choose a binding site.");

    int randomIndex = Rand::randInteger(0, _possibleBindings.size() - 1);
    auto it = _possibleBindings.begin();

    advance(it, randomIndex);

    return vector<tuple<CCylinder*, short>>{it->first, it->second};
}

void MotorBindingManager::appendpossibleBindings(short boundInt, CCylinder* ccyl1,
                                                 CCylinder* ccyl2, short site1, short site2){
	short _filamentType = ccyl1->getType();
	short _nfilamentType = ccyl2->getType();

	bool cndn1 = _filamentType == _filamentIDvec[0] && _nfilamentType == _filamentIDvec[1];
	bool cndn2 = _filamentType == _filamentIDvec[1] && _nfilamentType == _filamentIDvec[0];

	if(!cndn1 && !cndn2 && boundInt != _boundInt) return;

	auto t1 = make_tuple(ccyl1, site1);
	auto t2 = make_tuple(ccyl2, site2);
    floatingpoint oldN=numBindingSites();
    _possibleBindings.emplace(t1,t2);
    floatingpoint newN=numBindingSites();
    updateBindingReaction(oldN,newN);
}

void MotorBindingManager::printbindingsites() {
	cout<<"BINDINGSITES: CYL1(SIDX) CYL2(SIDX) SITE1 SITE2"<<endl;
	for(auto it1 = _possibleBindings.begin();it1!=_possibleBindings.end();
	    it1++) {
		auto cyl1 = get<0>(it1->first)->getCylinder();
		auto bs1 = get<1>(it1->first);
		auto cyl2 = get<0>(it1->second)->getCylinder();
		auto bs2 = get<1>(it1->second);
		cout<<cyl1->getStableIndex()<<" "
		    <<cyl2->getStableIndex()<<" "<<bs1<<" "<<bs2<<endl;
	}
}
#endif
//Deprecated.
bool MotorBindingManager::isConsistent() {
    return true;
}
#if defined(NLSTENCILLIST) || defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
void MotorBindingManager::addPossibleBindingsstencil(CCylinder* cc) {
    short _filamentType = cc->getType();
    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
        bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++) {
        addPossibleBindingsstencil(cc, *bit);
    }
}
void MotorBindingManager::addPossibleBindingsstencil(CCylinder* cc, short bindingSite) {
#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
//    cout<<"Adding "<<cc->getCylinder()->getID()<<" "<<bindingSite<<endl;
    auto HManager = _compartment->getHybridBindingSearchManager();
    HManager->addPossibleBindingsstencil(_Hbsmidvec,cc,bindingSite);
    //    HManager->checkoccupancySIMD(_Hbsmidvec);
#else
    if (SysParams::INITIALIZEDSTATUS ) {
        short complimentaryfID;
        short _filamentType = cc->getType();
        if (_filamentType != _filamentIDvec[0] && _filamentType != _filamentIDvec[1])
            return;
        else if (_filamentType == _filamentIDvec[0]) complimentaryfID = _filamentIDvec[1];
        else complimentaryfID = _filamentIDvec[0];

        //if we change other managers copy number
        vector<MotorBindingManager *> affectedManagers;

        //add valid binding sites
        if (areEqual(cc->getCMonomer(bindingSite)->speciesBound(
                SysParams::Chemistry().motorBoundIndex[_filamentType])->getN(),
                        (floatingpoint) 1.0)) {

            //loop through neighbors
            //now re add valid based on CCNL
            vector<Cylinder *> Neighbors;

            Neighbors = _neighborLists[_nlIndex]->getNeighborsstencil(cc->getCylinder());
            for (auto cn : Neighbors) {

                Cylinder *c = cc->getCylinder();

                if (cn->getParent() == c->getParent()) continue;
                if (cn->getType() != complimentaryfID) continue;

                auto ccn = cn->getCCylinder();

                for (auto it = SysParams::Chemistry().bindingSites[complimentaryfID].begin();
                        it !=
                        SysParams::Chemistry().bindingSites[complimentaryfID].end(); it++) {

                    if (areEqual(ccn->getCMonomer(*it)->speciesBound(
                            SysParams::Chemistry().motorBoundIndex[complimentaryfID])->getN(),
                                    (floatingpoint) 1.0f)) {

                        //check distances..
                        auto mp1 = (float) bindingSite /
                                    SysParams::Geometry().cylinderNumMon[_filamentType];
                        auto mp2 = (float) *it /
                                    SysParams::Geometry().cylinderNumMon[complimentaryfID];

                        auto x1 = c->getFirstBead()->vcoordinate();
                        auto x2 = c->getSecondBead()->vcoordinate();
                        auto x3 = cn->getFirstBead()->vcoordinate();
                        auto x4 = cn->getSecondBead()->vcoordinate();

                        auto m1 = midPointCoordinate(x1, x2, mp1);
                        auto m2 = midPointCoordinate(x3, x4, mp2);

                        floatingpoint distsq = twoPointDistancesquared(m1, m2);

                        if (distsq > _rMaxsq || distsq < _rMinsq) continue;

                        auto t1 = tuple<CCylinder *, short>(cc, bindingSite);
                        auto t2 = tuple<CCylinder *, short>(ccn, *it);

                        //add in correct order
                        if (c->getId() > cn->getId()) {
                            _possibleBindingsstencil.emplace(t1, t2);
                        } else {
                            //add in this compartment
                            if (cn->getCompartment() == _compartment) {

                                _possibleBindingsstencil.emplace(t2, t1);
                            }
                                //add in other
                            else {

                                auto m = (MotorBindingManager *) cn->getCompartment()->
                                        getFilamentBindingManagers()[_mIndex].get();

                                affectedManagers.push_back(m);

                                m->_possibleBindingsstencil.emplace(t2, t1);
                            }
                        }
                    }
                }
            }
        }

        //update affected
        for (auto m : affectedManagers) {

            int oldNOther = m->_bindingSpecies->getN();
            int newNOther = m->numBindingSitesstencil();

            m->updateBindingReaction(oldNOther, newNOther);
        }

        //update this manager
        int oldN = _bindingSpecies->getN();
        int newN = numBindingSitesstencil();

        updateBindingReaction(oldN, newN);
    }
#endif
}
void MotorBindingManager::updateAllPossibleBindingsstencil() {
#ifdef NLSTENCILLIST
    _possibleBindingsstencil.clear();
    int offset = 0;
            //SysParams::Mechanics().bsoffsetvec.at(_filamentType);
    floatingpoint min1,min2,max1,max2;
    bool status1 = true;
    bool status2 = true;
    floatingpoint minveca[2];
    floatingpoint maxveca[2];

    int maxnbs = SysParams::Chemistry().maxbindingsitespercylinder;
    floatingpoint* cylsqmagnitudevector = SysParams::Mechanics().cylsqmagnitudevector;
    auto boundstate = SysParams::Mechanics().speciesboundvec;

    const auto& cylinderInfoData = Cylinder::getDbData();

    int counter1 = 0;
    int totalneighbors = 0;

    int Ncylincmp =  _compartment->getCylinders().size();
    int* cindexvec = new int[Ncylincmp];
    vector<vector<int>> ncindices;
    vector<int>ncindex;
    long id = 0;
    for(auto c : _compartment->getCylinders()){
        cindexvec[id] = c->getStableIndex();
        totalneighbors += _neighborLists[_nlIndex]->getNeighborsstencil(c).size();
        id++;
        for (auto cn : _neighborLists[_nlIndex]->getNeighborsstencil(c)) {
            if(c->getId() > cn->getId())
                ncindex.push_back(cn->getStableIndex());
            else
                counter1++;
        }
        ncindices.push_back(ncindex);
        ncindex.clear();
    }
    for(int i=0;i<Ncylincmp;i++){
        int cindex = cindexvec[i];
        const auto& c = cylinderInfoData[cindex];

        short _filamentType = c.type;
        short complimentaryfID;
        if (_filamentType != _filamentIDvec[0] && _filamentType != _filamentIDvec[1])
        	return;
        else if (_filamentType == _filamentIDvec[0]) complimentaryfID = _filamentIDvec[1];
        else complimentaryfID = _filamentIDvec[0];
        short nbs1 = SysParams::Chemistry().bindingSites[_filamentType].size();
        short nbs2 = SysParams::Chemistry().bindingSites[complimentaryfID].size();

        const auto& x1 = Bead::getStableElement(c.beadIndices[0])->coord;
        const auto& x2 = Bead::getStableElement(c.beadIndices[1])->coord;
        floatingpoint X1X2[3] ={x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2]};
        int* cnindices = ncindices[i].data();
        for(int arraycount = 0; arraycount < ncindices[i].size();arraycount++){
            int cnindex = cnindices[arraycount];
            const auto& cn = cylinderInfoData[cnindex];
//            if(c.ID < cn.ID) {counter++; continue;} commented as the above vector does
// not contain ncs that will fail this cndn.
            if(c.filamentId == cn.filamentId){
                continue;}
            if(cn.type != complimentaryfID){
                 continue;}

            const auto& x3 = Bead::getStableElement(cn.beadIndices[0])->coord;
            const auto& x4 = Bead::getStableElement(cn.beadIndices[1])->coord;
            floatingpoint X1X3[3] = {x3[0] - x1[0], x3[1] - x1[1], x3[2] - x1[2]};
            floatingpoint X3X4[3] = {x4[0] - x3[0], x4[1] - x3[1], x4[2] - x3[2]};
            floatingpoint X1X3squared = sqmagnitude(X1X3);
            floatingpoint X1X2squared = cylsqmagnitudevector[cindex];
            floatingpoint X1X3dotX1X2 = scalarprojection(X1X3, X1X2);
            floatingpoint X3X4squared = cylsqmagnitudevector[cnindex];
            floatingpoint X1X3dotX3X4 = scalarprojection(X1X3,X3X4);
            floatingpoint X3X4dotX1X2 = scalarprojection(X3X4, X1X2);
            for(int posfil =0; posfil<nbs1;posfil++) {
            //now re add valid binding sites
                if (areEqual(boundstate[2][offset + maxnbs * cindex + posfil],
                        (floatingpoint)1.0)) {

                    floatingpoint mp1, minparamcyl2, maxparamcyl2;
                    if (_filamentType == _filamentIDvec[0]){
                        mp1 = bindingsites1.at(posfil);
                        minparamcyl2 = bindingsites2.front();
                        maxparamcyl2 = bindingsites2.back();
                    }
                    else {
                        mp1 = bindingsites2.at(posfil);
                        minparamcyl2 = bindingsites1.front();
                        maxparamcyl2 = bindingsites1.back();
                    }
                floatingpoint A = X3X4squared;
                floatingpoint B = 2 * X1X3dotX3X4 - 2 * mp1 * X3X4dotX1X2;
                floatingpoint C = X1X3squared + mp1 * mp1 * X1X2squared - 2 * mp1 * X1X3dotX1X2;
                floatingpoint C1 = C - _rMinsq;
                floatingpoint C2 = C - _rMaxsq;
                floatingpoint b2m4ac1 = B*B - 4*A*C1;
                floatingpoint b2m4ac2 = B*B - 4*A*C2;
                status1 = b2m4ac1 < 0;
                status2 = b2m4ac2 < 0;
                if(status1 && status2) {
                    continue;}
                if(!status1){
                    min1 = (-B + sqrt(b2m4ac1))/(2*A);
                    min2 = (-B - sqrt(b2m4ac1))/(2*A);
                    if(min1<min2) {
                        minveca[0] = (min1);
                        minveca[1] = (min2);
                    }
                    else{
                        minveca[0] = (min2);
                        minveca[1] = (min1);
                    }
                    if(minveca[0]< minparamcyl2 && minveca[1] > maxparamcyl2) {
                        continue;
                    }
                }
                if(!status2){
                    max1 = (-B + sqrt(b2m4ac2))/(2*A);
                    max2 = (-B - sqrt(b2m4ac2))/(2*A);
                    if(max1<max2) {
                        maxveca[0] = (max1);
                        maxveca[1] = (max2);
                    }
                    else{
                        maxveca[0] = (max2);
                        maxveca[1] = (max1);
                    }
                    if(maxveca[0] > maxparamcyl2 || maxveca[1] <minparamcyl2){
                        continue;
                    }
                }
                    for(int posCfil = 0; posCfil<nbs2;posCfil++){
                    if (areEqual(boundstate[2][offset + maxnbs * cnindex + posCfil], (floatingpoint)1.0)) {

                        //check distances..
                        floatingpoint mp2;
                        if (complimentaryfID == _filamentIDvec[0])
                            mp2 = bindingsites1.at(posCfil);
                        else
                            mp2 = bindingsites2.at(posCfil);

                        if(!status2) {
                            if (mp2 < maxveca[0] || mp2 > maxveca[1]) {
                                {
                                continue;}
                            }
                        }
                        if(!status1){
                            if (mp2 > minveca[0] && mp2 < minveca[1]) {
                                {
                                    continue;}
                            }
                        }

                        auto it1 = SysParams::Chemistry().bindingSites[_filamentType][posfil];
                        auto it2 = SysParams::Chemistry()
                                .bindingSites[complimentaryfID][posCfil];
                        auto t1 = tuple<CCylinder *, short>(cylinderInfoData[cindex].chemCylinder, it1);
                        auto t2 = tuple<CCylinder *, short>(cylinderInfoData[cnindex].chemCylinder, it2);
                        //add in correct order
                            _possibleBindingsstencil.emplace(t1, t2);
                    }
                }
            }
            }
        }
    }
    int oldN = _bindingSpecies->getN();
    int newN = numBindingSitesstencil();
    updateBindingReaction(oldN, newN);
    /*std::cout<<"Motor consistency "<<isConsistent()<<endl;*/
    delete[] cindexvec;
#else
    cout<<"Erroneous function call. Please check compiler macros. Exiting."<<endl;
    exit(EXIT_FAILURE);
#endif
}
void MotorBindingManager::appendpossibleBindingsstencil(short boundInt, CCylinder* ccyl1,
		CCylinder* ccyl2, short site1, short site2){

	short _filamentType = ccyl1->getType();
	short _nfilamentType = ccyl2->getType();

	bool cndn1 = _filamentType == _filamentIDvec[0] && _nfilamentType == _filamentIDvec[1];
	bool cndn2 = _filamentType == _filamentIDvec[1] && _nfilamentType == _filamentIDvec[0];

	if(!cndn1 && !cndn2 && boundInt != _boundInt) return;

	#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
	auto HManager = _compartment->getHybridBindingSearchManager();
	HManager->appendPossibleBindingsstencil(_Hbsmidvec, ccyl1, ccyl2, site1, site2);
	#else
	auto t1 = make_tuple(ccyl1, site1);
	auto t2 = make_tuple(ccyl2, site2);
	floatingpoint oldN=numBindingSitesstencil();
	_possibleBindingsstencil.emplace(t1,t2);
	floatingpoint newN=numBindingSitesstencil();
	updateBindingReaction(oldN,newN);
	#endif
}
void MotorBindingManager::removePossibleBindingsstencil(CCylinder* cc) {

    short _filamentType = cc->getType();

    for(auto bit = SysParams::Chemistry().bindingSites[_filamentType].begin();
        bit != SysParams::Chemistry().bindingSites[_filamentType].end(); bit++)
        removePossibleBindingsstencil(cc, *bit);
}
void MotorBindingManager::removePossibleBindingsstencil(CCylinder* cc, short bindingSite) {
#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
//    cout<<"Removing "<<cc->getCylinder()->getID()<<" "<<bindingSite<<endl;
    auto HManager = _compartment->getHybridBindingSearchManager();
    HManager->removePossibleBindingsstencil(_Hbsmidvec, cc, bindingSite);
/*    for(auto C:SubSystem::getstaticgrid()->getCompartments()){
        C->getHybridBindingSearchManager()->checkoccupancySIMD(_Hbsmidvec);
    }*/
#else

    short _filamentType = cc->getType();
    short complimentaryfID;
    if (_filamentType != _filamentIDvec[0] && _filamentType != _filamentIDvec[1]) return;
    else if (_filamentType == _filamentIDvec[0]) complimentaryfID = _filamentIDvec[1];
    else complimentaryfID = _filamentIDvec[0];

    //if we change other managers copy number
    vector<MotorBindingManager*> affectedManagers;

    //remove all tuples which have this ccylinder as key
    auto t = tuple<CCylinder*, short>(cc, bindingSite);
    _possibleBindingsstencil.erase(t);

    //remove all tuples which have this as value
    for (auto it = _possibleBindingsstencil.begin(); it != _possibleBindingsstencil.end();) {

        if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite)
        _possibleBindingsstencil.erase(it++);
        else ++it;
    }

    int oldN = _bindingSpecies->getN();
    int newN = numBindingSitesstencil();

    updateBindingReaction(oldN, newN);

    //remove all neighbors which have this binding site pair
    for (auto cn : _neighborLists[_nlIndex]->getNeighborsstencil(cc->getCylinder())) {

        if(cn->getType() != complimentaryfID) continue;

        if(cn->getCompartment() != _compartment) {

            auto m = (MotorBindingManager*)cn->getCompartment()->
                    getFilamentBindingManagers()[_mIndex].get();

            if(find(affectedManagers.begin(), affectedManagers.end(), m) == affectedManagers.end())
                affectedManagers.push_back(m);
        }
    }
//#endif
    //remove, update affected
    for(auto m : affectedManagers) {

        for (auto it = m->_possibleBindingsstencil.begin(); it !=
                m->_possibleBindingsstencil.end(); ) {

            if (get<0>(it->second) == cc && get<1>(it->second) == bindingSite)
                m->_possibleBindingsstencil.erase(it++);

            else ++it;
        }

        int oldNOther = m->_bindingSpecies->getN();
        int newNOther = m->numBindingSitesstencil();

        m->updateBindingReaction(oldNOther, newNOther);
    }
#endif
}
void MotorBindingManager::crosscheck(){
    cout<<"Motor NLORIGINAL size "<<_possibleBindings.size()<<" NLSTENCIL size "
        <<_possibleBindingsstencil.size()<<endl;
    if(_possibleBindings.size() != _possibleBindingsstencil.size()) {
        cout << "Motor.. The two methods compared do not yield the same number of "
        "binding sites" << endl;
        exit(EXIT_FAILURE);
    }
    short matches = 0;

    for(auto it1 = _possibleBindings.begin();it1!=_possibleBindings.end();it1++){
        short matchperelement = 0;

        auto cyl1_o = get<0>(it1->first)->getCylinder();
        auto bs1_o =  get<1>(it1->first);
        auto cyl2_o = get<0>(it1->second)->getCylinder();
        auto bs2_o =  get<1>(it1->second);
        //
        short sum = 0;
        for(auto it2 = _possibleBindingsstencil.begin();it2!=_possibleBindingsstencil.end();
            it2++){
            auto cyl1_s = get<0>(it2->first)->getCylinder();
            auto bs1_s =  get<1>(it2->first);
            auto cyl2_s = get<0>(it2->second)->getCylinder();
            auto bs2_s =  get<1>(it2->second);
            sum = 0;
            if(cyl1_o->getId() == cyl1_s->getId() )
                sum++;
            if(cyl2_o->getId() == cyl2_s->getId() )
                sum++;
            if(bs1_o == bs1_s )
                sum++;
            if(bs2_o == bs2_s )
                sum++;
            if(sum == 4) {
                //                cout << "match found" << endl;
                if(matchperelement == 0 ) {
                    matches++;
                    matchperelement++;
                } else
                cout<<"ERROR. Multiple matches for chosen binding site pair in "
                "STENCILLIST. Check stencillist.."<<endl;
            }
        }
    }
    std::cout<<"Motor possible bindings size "<<_possibleBindings.size()<<" Total matches"
            " "<<
             matches<<endl;
    if(_possibleBindings.size() != matches || _possibleBindings.size() !=
       _possibleBindingsstencil.size()){
        cout<<"Motor. All binding site pairs are not found in Stencil list"<<endl;
        exit(EXIT_FAILURE);
    }
}
vector<tuple<CCylinder*, short>> MotorBindingManager::chooseBindingSitesstencil() {
#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
    auto HManager = _compartment->getHybridBindingSearchManager();
    return HManager->chooseBindingSitesstencil(_Hbsmidvec);
#else
    assert((_possibleBindingsstencil.size() != 0)
           && "Major bug: Linker binding manager should not have zero binding \
                   sites when called to choose a binding site.");

    int randomIndex = Rand::randInteger(0, _possibleBindingsstencil.size() - 1);
    auto it = _possibleBindingsstencil.begin();

    advance(it, randomIndex);

    return vector<tuple<CCylinder*, short>>{it->first, it->second};
#endif
}

void MotorBindingManager::clearpossibleBindingsstencil() {
    #ifdef NLSTENCILLIST
    floatingpoint oldN=numBindingSitesstencil();
        _possibleBindingsstencil.clear();
        updateBindingReaction(oldN,0);
    #else
    auto HManager = _compartment->getHybridBindingSearchManager();
    HManager->clearPossibleBindingsstencil(_Hbsmidvec);

    #endif
}
int MotorBindingManager::numBindingSitesstencil() {
#ifdef NLSTENCILLSIT
    return _possibleBindingsstencil.size();
#else
    auto HManager = _compartment->getHybridBindingSearchManager();
    return HManager->numBindingSitesstencil(_Hbsmidvec);
#endif

}

void MotorBindingManager::printbindingsitesstencil() {
	#ifdef NLSTENCILLSIT
	cout<<"BINDINGSITES: CYL1(SIDX) CYL2(SIDX) SITE1 SITE2"<<endl;
	for(auto it1 = _possibleBindingsstencil.begin();it1!=_possibleBindingsstencil.end();
	    it1++) {
		auto cyl1 = get<0>(it1->first)->getCylinder();
		auto bs1 = get<1>(it1->first);
		auto cyl2 = get<0>(it1->second)->getCylinder();
		auto bs2 = get<1>(it1->second);
		cout<<cyl1->getStableIndex()<<" "
		    <<cyl2->getStableIndex()<<" "<<bs1<<" "<<bs2<<endl;
	}
	#else
	auto HManager = _compartment->getHybridBindingSearchManager();
	HManager->printbindingsitesstencil(_Hbsmidvec);
	#endif
}
#endif
#ifdef CUDAACCL_NL
void MotorBindingManager::assigncudavars() {

//    if(gpu_numpairs == NULL) {
        CUDAcommon::handleerror(cudaMalloc((void **) &gpu_numpairs, sizeof(int)), "cuda data transfer", "BindingManager.cu");
//    int n[1];
//    n[0] = 0;
//    CUDAcommon::handleerror(cudaMemcpy(gpu_numpairs, n, sizeof(int), cudaMemcpyHostToDevice));
//    }

//    if(gpu_rminmax == NULL) {
        CUDAcommon::handleerror(cudaMalloc((void **) &gpu_rminmax, 2 * sizeof(floatingpoint)), "cuda data transfer", " "
                "BindingManager.cu");
        floatingpoint dist[2];
        dist[0] = _rMin;
        dist[1] = _rMax;
        CUDAcommon::handleerror(cudaMemcpy(gpu_rminmax, dist, 2 * sizeof(floatingpoint), cudaMemcpyHostToDevice));
//    }

//    delete dist;
}
void MotorBindingManager::freecudavars() {
    CUDAcommon::handleerror(cudaFree(gpu_rminmax),"cudaFree", "BindingManager");
    CUDAcommon::handleerror(cudaFree(gpu_numpairs),"cudaFree", "BindingManager");
}
#endif

SubSystem* FilamentBindingManager::_subSystem = 0;

vector<CylinderCylinderNL*> LinkerBindingManager::_neighborLists;
vector<CylinderCylinderNL*> MotorBindingManager::_neighborLists;

} // namespace medyan

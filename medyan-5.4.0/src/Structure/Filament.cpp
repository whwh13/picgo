
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

#include <time.h>
#include <random>
#include <math.h>
#include <algorithm>
#include <functional>
#include <numeric>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/range/numeric.hpp>
#include <fstream>

#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"

#include "SubSystem.h"

#include "SysParams.h"
#include "MathFunctions.h"
#include "Controller/GController.h"
#include "Rand.h"

#ifdef CROSSCHECK_CYLINDER
#include "Controller/CController.h"
#endif

namespace medyan {
using namespace mathfunc;

Histogram* Filament::_turnoverTimes;

Filament::Filament(SubSystem* s, short filamentType, const vector<floatingpoint>& position,
                   const vector<floatingpoint>& direction, bool nucleation, bool branch)

    : Trackable(), _subSystem(s), _filType(filamentType) {

    if(getId() >= SysParams::Chemistry().maxStableIndex) {
        LOG(ERROR) << "Filament ID assigned ("<< getId()<<
                   ") equals/exceeds the maximum permissible value "
                   "("<<SysParams::Chemistry().maxStableIndex<<"). Exiting."<<endl;
        throw std::logic_error("Max value reached");
    }
    //create beads
    Bead* b1 = _subSystem->addTrackable<Bead>(position, this, 0);
    b1->monomerSerial = 0;
    
    //choose length
    floatingpoint length = 0.0;
    int numberOfMonomers = 0;
    
    if(branch) {
        length = SysParams::Geometry().monomerSize[_filType];
        numberOfMonomers = 1;
    }
    else if(nucleation) {
        length = SysParams::Geometry().minCylinderSize[_filType];
        numberOfMonomers = SysParams::Geometry().minCylinderNumMon;
    }
    
    auto pos2 = nextPointProjection(position, length, direction);
        
    Bead* b2 = _subSystem->addTrackable<Bead>(pos2, this, 1);
    b2->monomerSerial = numberOfMonomers;
    
    //create cylindera
    Cylinder* c0 = _subSystem->addTrackable<Cylinder>(this, b1, b2, _filType, 0);
    
    c0->setPlusEnd(true);
    c0->setMinusEnd(true);
    _cylinderVector.push_back(c0);
        
    // set cylinder's filID
    c0->setFilID(getId());


    //set plus end marker
    _plusEndPosition = 1;

    // update reaction rates
    if(const auto& loadForceFunc = _subSystem->getCylinderLoadForceFunc()) {
        loadForceFunc(c0, ForceFieldTypes::LoadForceEnd::Plus);
        loadForceFunc(c0, ForceFieldTypes::LoadForceEnd::Minus);
        c0->updateReactionRates();
    }
}


Filament::Filament(SubSystem* s, short filamentType, const std::vector<Vec<3, FP>>& position,
                   int numBeads, const SimulConfig& conf, string projectionType)

    : Trackable(), _subSystem(s), _filType(filamentType) {

    
    //create a projection of beads
    std::vector<Vec<3, FP>> tmpBeadsCoord;

    const auto cylLength = conf.geoParams.cylinderSize[filamentType];
    
    //straight projection
    if(projectionType == "STRAIGHT")
        tmpBeadsCoord = straightFilamentProjection(position[0], position[1], numBeads, cylLength);
    //zigzag projection
    else if(projectionType == "ZIGZAG")
        tmpBeadsCoord = zigZagFilamentProjection(position[0], position[1], numBeads, cylLength);
    //arc projection
    else if(projectionType == "ARC")
        tmpBeadsCoord = arcFilamentProjection(position[0], position[1], numBeads);
    // Predefined coordinates.
    else if(projectionType == "PREDEFINED") {
        tmpBeadsCoord = position;
    }
   
    //create beads
    auto direction = normalizedVector(tmpBeadsCoord[1] - tmpBeadsCoord[0]);
        
    Bead* b1 = _subSystem->addTrackable<Bead>(vec2Vector(tmpBeadsCoord[0]), this, 0);
    b1->monomerSerial = 0;
    Bead* b2 = _subSystem->addTrackable<Bead>(vec2Vector(tmpBeadsCoord[1]), this, 1);
    b2->monomerSerial = conf.geoParams.cylinderNumMon[_filType];

    Cylinder* c0 = _subSystem->addTrackable<Cylinder>(this, b1, b2, _filType, 0,
                                                      false, false, true);
        
    c0->setPlusEnd(true);
    c0->setMinusEnd(true);
    _cylinderVector.push_back(c0);

    // set cylinder's filID
    c0->setFilID(getId());


    for (int i = 2; i < numBeads; ++i) {
        extendPlusEnd(tmpBeadsCoord[i]);
    }
        
    //set plus end marker
    _plusEndPosition = numBeads - 1;

    // update reaction rates
    if(const auto& loadForceFunc = _subSystem->getCylinderLoadForceFunc()) {
        loadForceFunc(_cylinderVector.back(), ForceFieldTypes::LoadForceEnd::Plus);
        _cylinderVector.back()->updateReactionRates();
        loadForceFunc(_cylinderVector.front(), ForceFieldTypes::LoadForceEnd::Minus);
        _cylinderVector.front()->updateReactionRates();
    }
}



Filament::~Filament() {
    
    //remove cylinders, beads from system
    for(auto &c : _cylinderVector) {

        _subSystem->removeTrackable<Bead>(c->getFirstBead());
        
        if(c->isPlusEnd())
            _subSystem->removeTrackable<Bead>(c->getSecondBead());

        #ifdef CROSSCHECK_CYLINDER
        cout<<"RemoveTrackable Cylinder "<<c->getId()<<" "<<c->getStableIndex() <<endl;
        #endif
        _subSystem->removeTrackable<Cylinder>(c);
    }
}


//Extend front for initialization
void Filament::extendPlusEnd(const Vec<3, FP>& coordinates) {
    
    Cylinder* cBack = _cylinderVector.back();
    cBack->setPlusEnd(false);
    
    int lpf = cBack->getPosition();
    
    Bead* b2 = cBack->getSecondBead();
    
    //create a new bead
    auto newBeadCoords = vec2Vector(coordinates);
    //create
    Bead* bNew = _subSystem->addTrackable<Bead>(newBeadCoords, this, b2->getPosition() + 1);
    // Currently, this function extends a full size cylinder
    bNew->monomerSerial = b2->monomerSerial + SysParams::Geometry().cylinderNumMon[_filType];

    Cylinder* c0 = _subSystem->addTrackable<Cylinder> (this, b2, bNew, _filType,
                                                       lpf + 1, false, false, true);
    c0->setPlusEnd(true);
    _cylinderVector.push_back(c0);
    
    // set cylinder's filID

    c0->setFilID(getId());


}

//Extend back for initialization
void Filament::extendMinusEnd(const Vec<3, FP>& coordinates) {

    Cylinder* cFront = _cylinderVector.front();
    cFront->setMinusEnd(false);
    
    int lpf = cFront->getPosition();
    
    Bead* b2 = cFront->getFirstBead();
    
    //create a new bead
    auto direction = twoPointDirection(b2->vcoordinate(), vec2Vector(coordinates));
    auto newBeadCoords = nextPointProjection(b2->vcoordinate(),
    SysParams::Geometry().cylinderSize[_filType], direction);
    
    //create
    Bead* bNew = _subSystem->addTrackable<Bead>(newBeadCoords, this, b2->getPosition() - 1);
    // Currently, this function extends a full size cylinder
    bNew->monomerSerial = b2->monomerSerial - SysParams::Geometry().cylinderNumMon[_filType];

    Cylinder* c0 = _subSystem->addTrackable<Cylinder>(this, bNew, b2, _filType,
                                                  lpf - 1, false, false, true);
    c0->setMinusEnd(true);
    _cylinderVector.push_front(c0);
    

    // set cylinder's filID
    c0->setFilID(getId());

}

//Initialize for restart
void Filament::initializerestart(vector<Cylinder*> cylinderVector,
        vector<restartCylData>& _rCDatavec) {

    if(SysParams::RUNSTATE){
        LOG(ERROR) << "initializerestart Function from Filament class can only be called "
                      "during restart phase. Exiting.";
        throw std::logic_error("Illegal function call pattern");
    }
    if(_cylinderVector.size())
        cout<<_cylinderVector.size()<<endl;
    for(auto cyl:cylinderVector) {
        auto rcdata = _rCDatavec[cyl->getStableIndex()];
        cyl->initializerestart( rcdata.totalmonomers, rcdata.endmonomerpos[0],
                                rcdata.endmonomerpos[1], rcdata.endstatusvec[0],
                                rcdata.endstatusvec[1], rcdata.endtypevec[0],
                                rcdata.endtypevec[1]);
	    _cylinderVector.push_back(cyl);
    }

    //set plus end marker
    _plusEndPosition = getPlusEndCylinder()->getSecondBead()->getPosition();

	// update reaction rates
	if(const auto& loadForceFunc = _subSystem->getCylinderLoadForceFunc()) {
		loadForceFunc(_cylinderVector.back(), ForceFieldTypes::LoadForceEnd::Plus);
		_cylinderVector.back()->updateReactionRates();
		loadForceFunc(_cylinderVector.front(), ForceFieldTypes::LoadForceEnd::Minus);
		_cylinderVector.front()->updateReactionRates();
	}

}

//extend front at runtime
void Filament::extendPlusEnd(short plusEnd) {

    chrono::high_resolution_clock::time_point mins, mine;

    mins = chrono::high_resolution_clock::now();
    Cylinder* cBack = _cylinderVector.back();
    
    int lpf = cBack->getPosition();
    
    Bead* b1 = cBack->getFirstBead();
    Bead* b2 = cBack->getSecondBead();
    
    //move last bead of last cylinder forward
    auto direction1 = twoPointDirection(b1->vcoordinate(), b2->vcoordinate());
    
    auto npp = nextPointProjection(b2->vcoordinate(),
    SysParams::Geometry().monomerSize[_filType], direction1);

    mine = chrono::high_resolution_clock::now();
	chrono::duration<floatingpoint> elapsed_time1(mine - mins);
	FilextendPlusendtimer1 += elapsed_time1.count();

    //create a new bead in same place as b2
    Bead* bNew = _subSystem->addTrackable<Bead>(npp, this, b2->getPosition() + 1);
    bNew->monomerSerial = b2->monomerSerial + 1;
    
    //transfer the same load force to new bead
    //(approximation until next minimization)
    bNew->loadForcesP = b2->loadForcesP;
    bNew->lfip = b2->lfip + 1;
    
    Cylinder* c0 = _subSystem->addTrackable<Cylinder>(this, b2, bNew, _filType,
                                                      lpf + 1, true);

    mins = chrono::high_resolution_clock::now();
    _cylinderVector.back()->setPlusEnd(false);
    _cylinderVector.push_back(c0);
    _cylinderVector.back()->setPlusEnd(true);

    #ifdef CROSSCHECK_CYLINDER
    CController::_crosscheckdumpFilechem <<"PlusEnd cylinder set"<<endl;
    #endif
    
    // set cylinder's filID

    c0->setFilID(getId());

    #ifdef CROSSCHECK_CYLINDER
    CController::_crosscheckdumpFilechem <<"Cyl FilID set"<<endl;
    #endif


    //get last cylinder, mark species
    CMonomer* m = _cylinderVector.back()->getCCylinder()->getCMonomer(0);
    m->speciesPlusEnd(plusEnd)->up();

    #ifdef CROSSCHECK_CYLINDER
    CController::_crosscheckdumpFilechem <<"Species PlusEnd set"<<endl;
    #endif

    //update reaction rates
    _cylinderVector.back()->updateReactionRates();

    #ifdef CROSSCHECK_CYLINDER
    CController::_crosscheckdumpFilechem <<"Reaction Rates updated"<<endl;
    #endif
    
    _deltaPlusEnd++;

/*    cout<<"Extend minus End Cylinder ID = "<<cBack->getId()<<endl;
    cBack->printSelf();*/

    mine = chrono::high_resolution_clock::now();
	chrono::duration<floatingpoint> elapsed_time2(mine - mins);
	FilextendPlusendtimer2 += elapsed_time2.count();
    #ifdef CROSSCHECK_CYLINDER
    CController::_crosscheckdumpFilechem <<"extendPlusend"<<endl;
    #endif
}

//extend back at runtime
void Filament::extendMinusEnd(short minusEnd) {
    
    Cylinder* cFront = _cylinderVector.front();
    int lpf = cFront->getPosition();
    
    Bead* b2 = cFront->getFirstBead();
    Bead* b1 = cFront->getSecondBead();
    
    //move last bead of last cylinder forward
    auto direction1 = twoPointDirection(b1->vcoordinate(), b2->vcoordinate());
    
    auto npp = nextPointProjection(b2->vcoordinate(),
    SysParams::Geometry().monomerSize[_filType], direction1);
    
    //create a new bead in same place as b2
    Bead* bNew = _subSystem->addTrackable<Bead>(npp, this, b2->getPosition() - 1);
    bNew->monomerSerial = b2->monomerSerial - 1;

    //transfer the same load force to new bead
    //(approximation until next minimization)
    bNew->loadForcesM = b2->loadForcesM;
    bNew->lfim = b2->lfim + 1;
    
    Cylinder* c0 = _subSystem->addTrackable<Cylinder>(this, bNew, b2, _filType,
                                                      lpf - 1, false, true);
    _cylinderVector.front()->setMinusEnd(false);
    _cylinderVector.push_front(c0);
    _cylinderVector.front()->setMinusEnd(true);

    #ifdef CROSSCHECK_CYLINDER
    CController::_crosscheckdumpFilechem <<"MinusEnd cylinder set"<<endl;
    #endif
    
    // set cylinder's filID

    c0->setFilID(getId());

    //get first cylinder, mark species
    auto newCCylinder = getCylinderVector().front()->getCCylinder();
    CMonomer* m = newCCylinder->getCMonomer(newCCylinder->getSize() - 1);
    
    m->speciesMinusEnd(minusEnd)->up();
#ifdef CROSSCHECK_CYLINDER
    CController::_crosscheckdumpFilechem <<"MinusEnd species set"<<endl;
#endif

    
    //update reaction rates
    _cylinderVector.front()->updateReactionRates();

    #ifdef CROSSCHECK_CYLINDER
    CController::_crosscheckdumpFilechem <<"Reaction rates updated"<<endl;
    #endif
    _deltaMinusEnd++;

/*    cout<<"Extend plus End Cylinder ID = "<<cFront->getId()<<endl;
    cFront->printSelf();*/
    #ifdef CROSSCHECK_CYLINDER
    CController::_crosscheckdumpFilechem <<"extendMinusend"<<endl;
    #endif
}

//Depolymerize front at runtime
void Filament::retractPlusEnd() {
    
    Cylinder* retCylinder = _cylinderVector.back();

/*    cout<<"Ret plus End Cylinder ID = "<<retCylinder->getId()<<endl;
    retCylinder->printSelf();*/

    _cylinderVector.pop_back();
    
    //transfer load forces
    Bead* bd = _cylinderVector.back()->getSecondBead();
    bd->loadForcesP = retCylinder->getSecondBead()->loadForcesP;
    bd->lfip = retCylinder->getSecondBead()->lfip - 1;
    
    _subSystem->removeTrackable<Bead>(retCylinder->getSecondBead());
    removeChild(retCylinder->getSecondBead());
    #ifdef CROSSCHECK_CYLINDER
    cout<<"RemoveTrackable Cylinder "<<retCylinder->getId()<<" "
                                                             ""<<retCylinder->getStableIndex()<<endl;
    #endif
    _subSystem->removeTrackable<Cylinder>(retCylinder);
    removeChild(retCylinder);
    
    _cylinderVector.back()->setPlusEnd(true);

    
    //update rates of new front
    _cylinderVector.back()->updateReactionRates();
    
    _deltaPlusEnd--;

    #ifdef CROSSCHECK_CYLINDER
    CController::_crosscheckdumpFilechem <<"retractPlusend"<<endl;
    #endif
}

void Filament::retractMinusEnd() {
    
    Cylinder* retCylinder = _cylinderVector.front();

/*    cout<<"Ret minus End Cylinder ID = "<<retCylinder->getId()<<endl;
    retCylinder->printSelf();*/

    _cylinderVector.pop_front();
    
    //transfer load forces
    Bead* bd = _cylinderVector.front()->getFirstBead();
    bd->loadForcesM = retCylinder->getFirstBead()->loadForcesM;
    bd->lfim = retCylinder->getFirstBead()->lfim - 1;
    
    _subSystem->removeTrackable<Bead>(retCylinder->getFirstBead());
    removeChild(retCylinder->getFirstBead());
    #ifdef CROSSCHECK_CYLINDER
    cout<<"RemoveTrackable Cylinder "<<retCylinder->getId()<<" "<<retCylinder->getStableIndex() <<endl;
    #endif
    _subSystem->removeTrackable<Cylinder>(retCylinder);
    removeChild(retCylinder);
    
    _cylinderVector.front()->setMinusEnd(true);
    
    //update rates of new back
    _cylinderVector.front()->updateReactionRates();
    
    _deltaMinusEnd--;
    
    ///If filament has turned over, mark as such
    if(_plusEndPosition == getMinusEndCylinder()->getFirstBead()->getPosition()) {
        
        //reset
        _plusEndPosition = getPlusEndCylinder()->getSecondBead()->getPosition();
        _turnoverTime = tau();
    }
    #ifdef CROSSCHECK_CYLINDER
    CController::_crosscheckdumpFilechem <<"retractMinusEnd"<<endl;
    #endif
}

void Filament::polymerizePlusEnd() {
    
    Cylinder* cBack = _cylinderVector.back();
    
    Bead* b1 = cBack->getFirstBead();
    Bead* b2 = cBack->getSecondBead();
    
    auto direction = twoPointDirection(b1->vcoordinate(), b2->vcoordinate());
    
    b2->coordinate() = vector2Vec<3, floatingpoint>(nextPointProjection(b2->vcoordinate(),
    SysParams::Geometry().monomerSize[_filType], direction));
    ++b2->monomerSerial;

    // Update cylinder data
    cBack->getCoordinate() = (floatingpoint)0.5 * (b1->coordinate() + b2->coordinate());
    
    //increment load
    b2->lfip++;
    
    //increase eq length, update
    floatingpoint newEqLen = cBack->getMCylinder()->getEqLength() +
                      SysParams::Geometry().monomerSize[_filType];
    cBack->getMCylinder()->setEqLength(_filType, newEqLen);
    
    //update rates of new back
    _cylinderVector.back()->updateReactionRates();

    _polyPlusEnd++;

/*    cout<<"Poly plus End Cylinder ID = "<<cBack->getId()<<endl;
    cBack->printSelf();*/
    #ifdef CROSSCHECK_CYLINDER
    CController::_crosscheckdumpFilechem <<"polymerizePlusend"<<endl;
    #endif

}

void Filament::polymerizeMinusEnd() {
    
    Cylinder* cFront = _cylinderVector.front();
    
    Bead* b1 = cFront->getFirstBead();
    Bead* b2 = cFront->getSecondBead();

    auto direction = twoPointDirection(b2->vcoordinate(), b1->vcoordinate());
    
    b1->coordinate() = vector2Vec<3, floatingpoint>(nextPointProjection(b1->vcoordinate(),
    SysParams::Geometry().monomerSize[_filType], direction));
    --b1->monomerSerial;

    // Update cylinder data
    cFront->getCoordinate() = (floatingpoint)0.5 * (b1->coordinate() + b2->coordinate());

    
    //increment load
    b1->lfim++;
    
    //increase eq length, update
    floatingpoint newEqLen = cFront->getMCylinder()->getEqLength() +
                      SysParams::Geometry().monomerSize[_filType];
    cFront->getMCylinder()->setEqLength(_filType, newEqLen);

    //update rates of new back
    _cylinderVector.front()->updateReactionRates();

    _polyMinusEnd++;

/*    cout<<"Poly minus End Cylinder ID = "<<cFront->getId()<<endl;
    cFront->printSelf();*/
    #ifdef CROSSCHECK_CYLINDER
    CController::_crosscheckdumpFilechem <<"polymerizeMinusend"<<endl;
    #endif

}

void Filament::depolymerizePlusEnd() {
    
    Cylinder* cBack = _cylinderVector.back();
    
    Bead* b1 = cBack->getFirstBead();
    Bead* b2 = cBack->getSecondBead();

    auto direction = twoPointDirection(b2->vcoordinate(), b1->vcoordinate());
    
    b2->coordinate() = vector2Vec<3, floatingpoint>(nextPointProjection(b2->vcoordinate(),
    SysParams::Geometry().monomerSize[_filType], direction));
    --b2->monomerSerial;

    // Update cylinder data
    cBack->getCoordinate() = (floatingpoint)0.5 * (b1->coordinate() + b2->coordinate());

    
    //increment load
    b2->lfip--;
    
    //decrease eq length, update
    floatingpoint newEqLen = cBack->getMCylinder()->getEqLength() -
                      SysParams::Geometry().monomerSize[_filType];
    cBack->getMCylinder()->setEqLength(_filType, newEqLen);

    //update rates of new back
    _cylinderVector.front()->updateReactionRates();
    
    _depolyPlusEnd++;

/*    cout<<"DePoly plus End Cylinder ID = "<<cBack->getId()<<endl;
    cBack->printSelf();*/
    #ifdef CROSSCHECK_CYLINDER
    CController::_crosscheckdumpFilechem <<"depolymerizePlusend"<<endl;
    #endif

}

void Filament::depolymerizeMinusEnd() {
    
    Cylinder* cFront = _cylinderVector.front();
    
    Bead* b1 = cFront->getFirstBead();
    Bead* b2 = cFront->getSecondBead();
    
    auto direction = twoPointDirection(b1->vcoordinate(), b2->vcoordinate());
    
    b1->coordinate() = vector2Vec<3, floatingpoint>(nextPointProjection(b1->vcoordinate(),
    SysParams::Geometry().monomerSize[_filType], direction));
    ++b1->monomerSerial;

    // Update cylinder data
    cFront->getCoordinate() = (floatingpoint)0.5 * (b1->coordinate() + b2->coordinate());

    
    b1->lfim--;
    
    //decrease eq length, update
    floatingpoint newEqLen = cFront->getMCylinder()->getEqLength() -
                      SysParams::Geometry().monomerSize[_filType];
    cFront->getMCylinder()->setEqLength(_filType, newEqLen);

    //update rates of new back
    _cylinderVector.front()->updateReactionRates();

    _depolyMinusEnd++;
/*    cout<<"DePoly minus End Cylinder ID = "<<cFront->getId()<<endl;
    cFront->printSelf();*/
    #ifdef CROSSCHECK_CYLINDER
    CController::_crosscheckdumpFilechem <<"depolymerizeMinusend"<<endl;
    #endif
}


void Filament::nucleate(short plusEnd, short filament, short minusEnd) {
    
    //chemically initialize species
    CCylinder* cc = _cylinderVector[0]->getCCylinder();
    int monomerPosition = SysParams::Geometry().cylinderNumMon[_filType] / 2 + 1;
    
    CMonomer* m1 = cc->getCMonomer(monomerPosition - 1);
    CMonomer* m2 = cc->getCMonomer(monomerPosition);
    CMonomer* m3 = cc->getCMonomer(monomerPosition + 1);
    
    //minus end
    m1->speciesMinusEnd(minusEnd)->up();
    
    //filament
    m2->speciesFilament(filament)->up();
    
    for(auto j : SysParams::Chemistry().bindingIndices[_filType])
        m2->speciesBound(j)->up();
    
    //plus end
    m3->speciesPlusEnd(plusEnd)->up();

    _nucleationReaction++;

}


Filament* Filament::sever(int cylinderPosition) {
    
    const int vectorPosition = cylinderPosition - _cylinderVector.front()->getPosition();
    
    //if vector position is zero, we can't sever. return null
    if(vectorPosition == 0) return nullptr;
    
    //if any of the cylinders are only one monomer long, we can't sever
    CCylinder* ccBack  = _cylinderVector[vectorPosition - 1]->getCCylinder();
    CMonomer* cm = ccBack->getCMonomer(ccBack->getSize() - 1);
    
    if(cm->activeSpeciesMinusEnd() != -1) return nullptr;

    // If the plus side segment has only the plus end monomer, do not sever.
    {
        CCylinder* pcc = _cylinderVector[vectorPosition]->getCCylinder();
        CMonomer* pcm = pcc->getCMonomer(0);
        if(pcm->activeSpeciesPlusEnd() != -1) return nullptr;
    }
    
    //create a new filament
    Filament* newFilament = _subSystem->addTrackable<Filament>(_subSystem, _filType);
    
    //Split the cylinder vector at position, transfer cylinders to new filament
    for(int i = 0; i < vectorPosition; ++i) {
        
        Cylinder* c = _cylinderVector.front();
        _cylinderVector.pop_front();

        newFilament->_cylinderVector.push_back(c);
        
        // Transfer the cylinder and its minus end bead to new parent.
        transferChild(c, *newFilament);
        transferChild(c->getFirstBead(), *newFilament);
    }
    // Get plus end cylinder of new filament, and minus end cylinder of current filament.
    auto c1 = newFilament->_cylinderVector.back();
    auto c2 = _cylinderVector.front();
    
    ///copy bead at severing point, attach to new filament
    Bead* oldB = c2->getFirstBead();
    Bead* newB = _subSystem->addTrackable<Bead>(*oldB);
    newB->monomerSerial = oldB->monomerSerial;
    
    //offset these beads by a little for safety
    auto msize = SysParams::Geometry().monomerSize[_filType];

    medyan::Vec< 3, floatingpoint > offsetCoord {
        (Rand::randInteger(0,1) ? -1 : +1) * Rand::randfloatingpoint(msize, 2 * msize),
        (Rand::randInteger(0,1) ? -1 : +1) * Rand::randfloatingpoint(msize, 2 * msize),
        (Rand::randInteger(0,1) ? -1 : +1) * Rand::randfloatingpoint(msize, 2 * msize),
    };
    
    oldB->coord += offsetCoord;
    newB->coord -= offsetCoord;
    
    //add bead
    c1->setSecondBead(newB);
    newFilament->addChild(unique_ptr<Component>(newB));
    
    //set plus and minus ends
    c1->setPlusEnd(true);
    c2->setMinusEnd(true);
    
    //mark the plus and minus ends of the new and old filament
    CCylinder* cc1 = c1->getCCylinder();
    CCylinder* cc2 = c2->getCCylinder();
    
    CMonomer* m1 = cc1->getCMonomer(cc1->getSize() - 1);
    CMonomer* m2 = cc2->getCMonomer(0);
    
    short filamentInt1 = m1->activeSpeciesFilament();
    short filamentInt2 = m2->activeSpeciesFilament();
    
    //plus end
    m1->speciesFilament(filamentInt1)->down();
    m1->speciesPlusEnd(filamentInt1)->up();
    
    for(auto j : SysParams::Chemistry().bindingIndices[_filType])
        m1->speciesBound(j)->down();
    
    //minus end
    m2->speciesFilament(filamentInt2)->down();
    m2->speciesMinusEnd(filamentInt2)->up();
    
    for(auto j : SysParams::Chemistry().bindingIndices[_filType])
        m2->speciesBound(j)->down();
    
    //remove any cross-cylinder rxns between these two cylinders
    cc1->removeCrossCylinderReactions(cc2);
    cc2->removeCrossCylinderReactions(cc1);

    _severingReaction++;
    _severingID.push_back(newFilament->getId());
    return newFilament;
}

std::vector<Vec<3, FP>> Filament::straightFilamentProjection(const Vec<3, FP>& v1, const Vec<3, FP>& v2, int numBeads, FP cylLength) {
    
    std::vector<Vec<3, FP>> coordinate;
    coordinate.reserve(numBeads);
    const auto dir = (v2 - v1) / distance(v1, v2);
    
    for (int i = 0; i < numBeads; ++i) {
        coordinate.push_back(v1 + dir * (i * cylLength));
    }
    return coordinate;
}

std::vector<Vec<3, FP>> Filament::zigZagFilamentProjection(const Vec<3, FP>& v1, const Vec<3, FP>& v2, int numBeads, FP cylLength) {
    
    std::vector<Vec<3, FP>> coordinate;
    coordinate.reserve(numBeads);
    const auto dir = (v2 - v1) / distance(v1, v2);
    const Vec<3, FP> altDir { -dir[1], dir[0], dir[2] };

    for (int i = 0; i < numBeads; ++i) {
        auto tmp = (i % 2 == 0)
            ? (i * cylLength) * dir + v1
            : (i * cylLength) * altDir + v1;
        coordinate.push_back(tmp);
    }
    return coordinate;
}

/// Create a projection
/// @note - created by Aravind 12/2014
void marsagila(vector<floatingpoint>&v) {
    
    floatingpoint d1,d2,d3;
    floatingpoint *x=new floatingpoint[3];
    d1=2*Rand::randfloatingpoint(0,1)-1;
    d2=2*Rand::randfloatingpoint(0,1)-1;
    d3=pow(d1,2)+pow(d2,2);
    
    while(d3>=1) {
        d1=2*Rand::randfloatingpoint(0,1)-1;
        d2=2*Rand::randfloatingpoint(0,1)-1;
        d3=pow(d1,2)+pow(d2,2);
    }
    
    x[0]=2.0*d1*pow((1.0-d3),0.5);
    x[1]=2.0*d2*pow(1.0-d3,0.5);
    x[2]=1.0-2.0*d3;
    v[0]=2.0*d1*pow((1.0-d3),0.5);
    v[1]=2.0*d2*pow(1.0-d3,0.5);
    v[2]=1.0-2.0*d3;
}

/// Matrix multiply
/// @note - created by Aravind 12/2014
void matrix_mul(boost::numeric::ublas::matrix<floatingpoint>&X,
                boost::numeric::ublas::matrix<floatingpoint>&Y,
                boost::numeric::ublas::matrix<floatingpoint>&Z,
                vector<floatingpoint>&x,vector<floatingpoint>&y,
                vector<floatingpoint>&z,int nbeads,
                std::vector<Vec<3, FP>> &coordinate) {
    
    int t,i;
    floatingpoint dt,length,cyl_length,sum;
    vector<int> id;
    vector<floatingpoint> dx,dy,dz,dx2,dy2,dz2,length2,dxdy2;
    using namespace boost::numeric::ublas;
    matrix<floatingpoint> B(4,4),B2(1,4),temp1(4,1),dummy(1,1),temp2(4,1),temp3(4,1);
    
    // B
    B(0,0)=1; B(0,1)=0; B(0,2)=0; B(0,3)=0;
    B(1,0)=-3; B(1,1)=3; B(1,2)=0; B(1,3)=0;
    B(2,0)=3; B(2,1)=-6; B(2,2)=3; B(2,3)=0;
    B(3,0)=-1; B(3,1)=3; B(3,2)=-3; B(3,3)=1;
    
    axpy_prod(B,X,temp1);
    axpy_prod(B,Y,temp2);
    axpy_prod(B,Z,temp3);
    B2(0,0)=1;
    
    for(t=0;t<=4000;t++) {
        dt=0.00025*t;
        B2(0,1)=dt;
        B2(0,2)=dt*dt;
        B2(0,3)=dt*dt*dt;
        axpy_prod(B2,temp1,dummy);
        x.push_back(dummy(0,0));
        axpy_prod(B2,temp2,dummy);
        y.push_back(dummy(0,0));
        axpy_prod(B2,temp3,dummy);
        z.push_back(dummy(0,0));
    }
    
    adjacent_difference(x.begin(),x.end(),back_inserter(dx));//dx
    adjacent_difference(y.begin(),y.end(),back_inserter(dy));//dy
    adjacent_difference(z.begin(),z.end(),back_inserter(dz));//dz
    
    transform(dx.begin(), dx.end(),dx.begin(),
              back_inserter(dx2), multiplies<floatingpoint>());
    transform(dy.begin(), dy.end(),dy.begin(),
              back_inserter(dy2), multiplies<floatingpoint>());
    transform(dz.begin(), dz.end(),dz.begin(),
              back_inserter(dz2), multiplies<floatingpoint>());
    
    //array of sum(dx^2+dy^2)
    transform(dx2.begin(),dx2.end(),dy2.begin(),
              back_inserter(dxdy2),plus<floatingpoint>());
    //array of sum(dx^2+dy^2+dz^2)
    transform(dxdy2.begin(),dxdy2.end(),dz2.begin(),
              back_inserter(length2),plus<floatingpoint>());
    
    std::vector<floatingpoint> tempLength;
    for(auto x: length2) tempLength.push_back(sqrt(x));
    length2 = tempLength; length2[0]=0.0;
    
    length = boost::accumulate(length2, 0.0);//arc length.
    
    //making equal divisions.
    i=0;sum=0.0;id.push_back(0.0);
    cyl_length=length/(nbeads-1);
    
    while(i<=4000) {
        sum+=length2[i];
        if(sum>=cyl_length||i==4000) {
            id.push_back(i);
            sum=0.0;
        }
        i++;
    }
    
    for(i=0;i<id.size();i++) {
        coordinate.push_back({
            x[id[i]],
            y[id[i]],
            z[id[i]],
        });
    }
}

void arcOutward(vector<floatingpoint>&v1,vector<floatingpoint>&v2, const Vec<3, FP>& v1in, const Vec<3, FP>& v2in) {
    
    vector<floatingpoint> center,tempv1,tempv2,temp2,temp3(3),
                   temp4(3),mid,mid2(3),mid3(3),temp5;
    
    center = GController::getCenter();
    
    // point v[0]
    temp3[0]=v1in[0];
    temp3[1]=v1in[1];
    temp3[2]=v1in[2];
    
    //point v[1]
    temp4[0]=v2in[0];
    temp4[1]=v2in[1];
    temp4[2]=v2in[2];
    
    mid=midPointCoordinate(temp3, temp4, 0.5);
    mid2=midPointCoordinate(mid, temp4, 0.5);
    mid3=midPointCoordinate(mid, temp3, 0.5);
    
    //vector between v[1] and center stored in tempv1
    std::transform(v2in.begin(), v2in.end(), center.begin(),
                   std::back_inserter(tempv1), std::minus<floatingpoint>());
    
    floatingpoint dist=twoPointDistance(center, temp4);
    dist=300/dist;
    
    std::transform(tempv1.begin(),tempv1.end(),tempv1.begin(),
                   [&](auto x) { return x * dist; });
    std::transform(mid.begin(), mid.end(), tempv1.begin(),
                   std::back_inserter(v1), std::plus<floatingpoint>());
    
    //vector between v[0] and center stored in tempv2
    std::transform(v1in.begin(), v1in.end(), center.begin(),
                   std::back_inserter(tempv2), std::minus<floatingpoint>());
    
    dist=twoPointDistance(center, temp3);
    dist=100/dist;
    
    std::transform(tempv2.begin(),tempv2.end(),tempv2.begin(),
                   [&](auto x) { return x * dist; });
    std::transform(mid3.begin(), mid3.end(), tempv2.begin(),
                   std::back_inserter(v2), std::plus<floatingpoint>());
}

std::vector<Vec<3, FP>> Filament::arcFilamentProjection(const Vec<3, FP>& v1, const Vec<3, FP>& v2, int numBeads) {
    
    using namespace boost::numeric::ublas;

    std::vector<floatingpoint> X3,x3(3),x4(3),X4,x,y,z;
    matrix<floatingpoint> C(3,3),B(4,4),X(4,1),Y(4,1),Z(4,1);
    std::vector<Vec<3, FP>> coordinates;

    arcOutward(X3,X4, v1, v2);
    X(0,0)=v1[0]; X(1,0)=X3[0]; X(2,0)=X4[0]; X(3,0)=v2[0];
    Y(0,0)=v1[1]; Y(1,0)=X3[1]; Y(2,0)=X4[1]; Y(3,0)=v2[1];
    Z(0,0)=v1[2]; Z(1,0)=X3[2]; Z(2,0)=X4[2]; Z(3,0)=v2[2];
    //
    matrix_mul(X,Y,Z,x,y,z,numBeads,coordinates);
    return coordinates;
}


void Filament::printSelf()const {
    
    cout << endl;
    
    cout << "Filament: ptr = " << this << endl;
    cout << "Filament ID = " << getId() << endl;
    cout << "Filament type = " << _filType << endl;
    
    cout << endl;
    cout << "Cylinder information..." << endl;
    
    for(auto c : _cylinderVector)
        c->printSelf();
    
    cout << endl;
    
}

bool Filament::isConsistent() {
    
    //check consistency of each individual cylinder
    for(auto &c : _cylinderVector) {
        if(!c->getCCylinder()->isConsistent()) {
         
            cout << "Cylinder at position " << c->getPosition()
                 << " is chemically inconsistent" << endl;
            return false;
        }
    }

    //check that it only has one plus end and one minus end
    int numPlusEnd = 0;
    int numMinusEnd = 0;
        
    for(auto &c : _cylinderVector) {
        
        for(int i = 0; i < c->getCCylinder()->getSize(); i++) {
            auto m = c->getCCylinder()->getCMonomer(i);
            
            if(m->activeSpeciesPlusEnd() != -1) numPlusEnd++;
            if(m->activeSpeciesMinusEnd() != -1) numMinusEnd++;
        }
    }
    if(numPlusEnd != 1) {
        cout << "This filament has more than one plus end species." << endl;
        return false;
    }
    if(numMinusEnd != 1) {
        cout << "This filament has more than one minus end species." << endl;
        return false;
    }
     
    return true;
}


species_copy_t Filament::countSpecies(short filamentType, const string& name) {
    
    species_copy_t copyNum = 0;
    
    for(auto f : getElements()) {
        
        if(f->getType() != filamentType) continue;
        
        //loop through the filament
        for(auto c : f->_cylinderVector) {
            
            for(int i = 0; i < c->getCCylinder()->getSize(); i++) {
                auto m = c->getCCylinder()->getCMonomer(i);
                
                //filament species
                int activeIndex = m->activeSpeciesFilament();
                
                if(activeIndex != -1) {
                    
                    auto s = m->speciesFilament(activeIndex);
                    string sname = SpeciesNamesDB::removeUniqueFilName(s->getName());
                    
                    if(sname == name)
                        copyNum += s->getN();
                    
                    continue;
                }
                
                //plus end species
                activeIndex = m->activeSpeciesPlusEnd();
                
                if(activeIndex != -1) {
                    
                    auto s = m->speciesPlusEnd(activeIndex);
                    string sname = SpeciesNamesDB::removeUniqueFilName(s->getName());
                    
                    if(sname == name)
                        copyNum += s->getN();
                    
                    continue;
                }
                
                //minus end species
                activeIndex = m->activeSpeciesMinusEnd();
                
                if(activeIndex != -1) {
                    
                    auto s = m->speciesMinusEnd(activeIndex);
                    string sname = SpeciesNamesDB::removeUniqueFilName(s->getName());
                    
                    if(sname == name)
                        copyNum += s->getN();
                    
                    continue;
                }
                
            }
        }
    }
    return copyNum;
}

floatingpoint Filament::FilextendPlusendtimer1 = 0.0;
floatingpoint Filament::FilextendPlusendtimer2 = 0.0;
floatingpoint Filament::FilextendPlusendtimer3 = 0.0;
floatingpoint Filament::FilextendMinusendtimer1 = 0.0;
floatingpoint Filament::FilextendMinusendtimer2 = 0.0;
floatingpoint Filament::FilextendMinusendtimer3 = 0.0;

} // namespace medyan

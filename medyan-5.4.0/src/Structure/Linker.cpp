
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

#include "Linker.h"

#include "Bead.h"
#include "Cylinder.h"
#include "ChemRNode.h"

#include "Controller/GController.h"
#include "SysParams.h"
#include "MathFunctions.h"
#include "Mechanics/CUDAcommon.h"

namespace medyan {
using namespace mathfunc;


void Linker::updateCoordinate() {
    
    auto x1 = _c1->getFirstBead()->vcoordinate();
    auto x2 = _c1->getSecondBead()->vcoordinate();
    auto x3 = _c2->getFirstBead()->vcoordinate();
    auto x4 = _c2->getSecondBead()->vcoordinate();
    
    auto m1 = midPointCoordinate(x1, x2, _c1->adjustedrelativeposition(_position1));
    auto m2 = midPointCoordinate(x3, x4, _c2->adjustedrelativeposition(_position2));
    
    coordinate = midPointCoordinate(m1, m2, 0.5);
}

Linker::Linker(
    Cylinder* c1, Cylinder* c2,
    short linkerType,
    int linkerSpeciesIndex1,
    int linkerSpeciesIndex2,
    floatingpoint position1, floatingpoint position2)

    : Trackable(true, true), _c1(c1), _c2(c2),
      _position1(position1), _position2(position2),
      _linkerType(linkerType), _birthTime(tau())
{
    using namespace std;

    // Initialize linker mechanics.
    //---------------------------------
    // Set stretching constant.
    if(!SysParams::Mechanics().LStretchingK.empty())
        mLinker_.kStretch = SysParams::Mechanics().LStretchingK[linkerType];

    // Set equilibrium length.
    auto x1 = _c1->getFirstBead()->vcoordinate();
    auto x2 = _c1->getSecondBead()->vcoordinate();
    auto x3 = _c2->getFirstBead()->vcoordinate();
    auto x4 = _c2->getSecondBead()->vcoordinate();

    auto m1 = midPointCoordinate(x1, x2, c1->adjustedrelativeposition(position1));
    auto m2 = midPointCoordinate(x3, x4, c2->adjustedrelativeposition(position2));
    mLinker_.eqLength = twoPointDistance(m1, m2);


    // Initialize coordinates and registration in compartments.
    //---------------------------------
    updateCoordinate();

    try {_compartment = &GController::getCompartment(coordinate);}
    catch (exception& e) {
        cout << e.what();
        
        printSelf();
        
        exit(EXIT_FAILURE);
    }

    // Initialize linker chemistry.
    //---------------------------------
    const int pos1 = int(position1 * SysParams::Geometry().cylinderNumMon[c1->getType()]);
    const int pos2 = int(position2 * SysParams::Geometry().cylinderNumMon[c2->getType()]);
  
    _cLinker = make_unique<CLinker>(linkerSpeciesIndex1, linkerSpeciesIndex2, _compartment, _c1->getCCylinder(), _c2->getCCylinder(), pos1, pos2);
    _cLinker->setLinker(this);
        
    
}

///@note - tracks lifetime data here
Linker::~Linker() noexcept {

//    floatingpoint lifetime = tau() - _birthTime;
//    
//    if(_lifetimes->getMax() > lifetime &&
//       _lifetimes->getMin() < lifetime)
//        _lifetimes->addValue(lifetime);

}


void Linker::updatePosition() {
    
    //update ccylinders
    _cLinker->setFirstCCylinder(_c1->getCCylinder());
    _cLinker->setSecondCCylinder(_c2->getCCylinder());
    
    updateCoordinate();
    
    Compartment* c;
    
    try {c = &GController::getCompartment(coordinate);}
    catch (exception& e) {
        cout << e.what();
        
        printSelf();
        
        exit(EXIT_FAILURE);
    }
    
    if(c != _compartment) {
	    chrono::high_resolution_clock::time_point mins, mine;
        _compartment = c;

        SpeciesBound* firstSpecies = _cLinker->getFirstSpecies();
        SpeciesBound* secondSpecies = _cLinker->getSecondSpecies();
        
        CLinker* clone = _cLinker->clone(c);
        setCLinker(clone);
        
        _cLinker->setFirstSpecies(firstSpecies);
        _cLinker->setSecondSpecies(secondSpecies);

	    mine = chrono::high_resolution_clock::now();
	    chrono::duration<floatingpoint> compartment_update(mine - mins);
	    CUDAcommon::tmin.timelinkerupdate += compartment_update.count();
	    CUDAcommon::tmin.callslinkerupdate++;
    }
    
    if(SysParams::RUNSTATE) {
        auto x1 = _c1->getFirstBead()->vcoordinate();
        auto x2 = _c1->getSecondBead()->vcoordinate();
        auto x3 = _c2->getFirstBead()->vcoordinate();
        auto x4 = _c2->getSecondBead()->vcoordinate();

        auto m1 = midPointCoordinate(x1, x2, _c1->adjustedrelativeposition(_position1));
        auto m2 = midPointCoordinate(x3, x4, _c2->adjustedrelativeposition(_position2));

    }
}

/// @note - The function uses the linker's stretching force at
/// the current state to change this rate. Does not consider
/// compression forces, only stretching.

void Linker::updateReactionRates() {

    //if no rate changers were defined, skip
    if(_unbindingChangers.empty()) return;
    
    //current force on linker
    floatingpoint force = max<floatingpoint>((floatingpoint)0.0, mLinker_.stretchForce);
    
    //get the unbinding reaction
    ReactionBase* offRxn = _cLinker->getOffReaction();

    //change the rate

    if(SysParams::RUNSTATE==false)
        offRxn->setRateMulFactor(0.0f, ReactionBase::RESTARTPHASESWITCH);
    else
        offRxn->setRateMulFactor(1.0f, ReactionBase::RESTARTPHASESWITCH);

#ifdef DETAILEDOUTPUT
    std::cout<<"Linker UB f "<<force<<" Rate "<<newRate<<" "<<coordinate[0]<<" "
            ""<<coordinate[1]<<" "
            ""<<coordinate[2]<<endl;
#endif
    if(_unbindingChangers.size() > 0) {
        float factor = _unbindingChangers[_linkerType]->getRateChangeFactor(force);
        offRxn->setRateMulFactor(factor, ReactionBase::mechanochemical);
        offRxn->updatePropensity();
    }
}


void Linker::printSelf()const {
    
    cout << endl;
    
    cout << "Linker: ptr = " << this << endl;
    cout << "Linker type = " << _linkerType << ", Linker ID = " << getId() << endl;
    cout << "Coordinates = " << coordinate[0] << ", " << coordinate[1] << ", " << coordinate[2] << endl;
    
    cout << "Position on first cylinder (floatingpoint) = " << _position1 << endl;
    cout << "Position on second cylinder (floatingpoint) = " << _position2 << endl;
    
    cout << "Birth time = " << _birthTime << endl;
    
    cout << endl;
    
    cout << "Associated species 1 = " << _cLinker->getFirstSpecies()->getName()
         << " , copy number = " << _cLinker->getFirstSpecies()->getN()
         << " , position on first cylinder (int) = " << _cLinker->getFirstPosition() << endl;
    
    cout << "Associated species 2 = " << _cLinker->getSecondSpecies()->getName()
         << " , copy number = " << _cLinker->getSecondSpecies()->getN()
         << " , position on second cylinder (int) = " << _cLinker->getSecondPosition() << endl;
    
    cout << endl;
    
    cout << "Associated cylinders (one and two): " << endl;
    _c1->printSelf();
    _c2->printSelf();
    
    cout << endl;
}

species_copy_t Linker::countSpecies(const string& name) {
    
    species_copy_t copyNum = 0;
    
    for(auto l : getElements()) {
        
        auto s = l->getCLinker()->getFirstSpecies();
        string sname = SpeciesNamesDB::removeUniqueFilName(s->getName());
        
        if(sname == name)
            copyNum += s->getN();
    }
    return copyNum;
}

vector<LinkerRateChanger*> Linker::_unbindingChangers;

Histogram* Linker::_lifetimes;

} // namespace medyan

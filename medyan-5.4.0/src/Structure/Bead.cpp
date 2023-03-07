
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

#include "Bead.h"

#include "Compartment.h"
#include "Composite.h"

#include "SysParams.h"
#include "Controller/GController.h"
#include "MathFunctions.h"

namespace medyan {
using namespace mathfunc;

std::vector<Bead*> Bead::_pinnedBeads;

Bead::Bead (vector<floatingpoint> v, Composite* parent, int position)
//add brforce, pinforce
    : Trackable(true, false, true, false),
      coord(vector2Vec<3>(v)),
      force{},
      coordinateP(v),
      brforce(3, 0), pinforce(3,0),
      _position(position), _birthTime(tau()) {

    parent->addChild(unique_ptr<Component>(this));
          
    loadForcesP = vector<floatingpoint>(SysParams::Geometry().cylinderNumMon[getType()], 0.0);
    loadForcesM = vector<floatingpoint>(SysParams::Geometry().cylinderNumMon[getType()], 0.0);

    //Find compartment
    try {_compartment = &GController::getCompartment(v);}
    catch (exception& e) {
        
        cout << e.what() << endl;
        
        printSelf();
        
        //also print parent info
        getParent()->printSelf();
        
        throw;
    }

}

Bead::Bead(Composite* parent, int position)
//add brforce, pinforce
    : Trackable(true, false, true, false),
    coord{},
    force{},
    coordinateP(3, 0),
    brforce(3, 0), pinforce(3,0), _position(position) {
    
    parent->addChild(unique_ptr<Component>(this));

}

void Bead::updatePosition() {
    
    try {GController::getCompartment(vec2Vector(coordinate()));}
    catch (exception& e) {
        
        //print exception
        cout << e.what() << endl;
        
        printSelf();
        
        //also print parent info
        getParent()->printSelf();
        
        throw;
    }
}

void Bead::printSelf()const {
    
    cout << endl;
    
    cout << "Bead: ptr = " << this << endl;
    cout << "Coordinates = " << coordinate()[0] << ", " << coordinate()[1] << ", " << coordinate()[2] << endl;
    cout << "Previous coordinates before minimization = " << coordinateP[0] << ", " << coordinateP[1] << ", " << coordinateP[2] << endl;
    cout << "Forces = " << force[0] << ", " << force[1] << ", " << force[2] << endl;

    cout << "Position on structure = " << _position << endl;
    cout << "Birth time = " << _birthTime << endl;
    
    
    cout << endl;
}

floatingpoint Bead::getLoadForcesP() {
    
    if (lfip < 0)
        return loadForcesP[0];
        
    if (lfip >= loadForcesP.size())
        return loadForcesP.back();
    
    else return loadForcesP[lfip];
}

floatingpoint Bead::getLoadForcesM() {
    
    if (lfim < 0)
        return loadForcesM[0];
    
    if (lfim >= loadForcesM.size())
        return loadForcesM.back();
    
    else return loadForcesM[lfim];
}

} // namespace medyan

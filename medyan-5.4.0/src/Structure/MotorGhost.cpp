
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

#include "MotorGhost.h"

#include "Filament.h"
#include "Cylinder.h"
#include "Bead.h"
#include "ChemRNode.h"

#include "Controller/GController.h"
#include "SysParams.h"
#include "MathFunctions.h"
#include "Mechanics/CUDAcommon.h"
#include "Rand.h"

namespace medyan {
using namespace mathfunc;

void MotorGhost::updateCoordinate() {
    
    auto x1 = _c1->getFirstBead()->vcoordinate();
    auto x2 = _c1->getSecondBead()->vcoordinate();
    auto x3 = _c2->getFirstBead()->vcoordinate();
    auto x4 = _c2->getSecondBead()->vcoordinate();
    
    auto m1 = midPointCoordinate(x1, x2, _c1->adjustedrelativeposition(_position1));
    auto m2 = midPointCoordinate(x3, x4, _c2->adjustedrelativeposition(_position2));
    
    coordinate = midPointCoordinate(m1, m2, 0.5);
}


MotorGhost::MotorGhost(
    Cylinder* c1, Cylinder* c2, short motorType,
    int motorSpeciesIndex1, int motorSpeciesIndex2,
    floatingpoint position1, floatingpoint position2,
    floatingpoint onRate, floatingpoint offRate)

    : Trackable(true, true),
      _c1(c1), _c2(c2),
      _position1(position1), _position2(position2),
      _motorType(motorType), _birthTime(tau()),
      _onRate(onRate), _offRate(offRate) {

    using namespace std;

    // Initialize motor heads.
    //---------------------------------
    //set number of heads by picking random int between maxheads and minheads
    _numHeads = Rand::randInteger(SysParams::Chemistry().motorNumHeadsMin[_motorType],
                                  SysParams::Chemistry().motorNumHeadsMax[_motorType]);
    
    if(!_unbindingChangers.empty())
        _numBoundHeads = _unbindingChangers[_motorType]->numBoundHeads(_onRate, _offRate, 0, _numHeads);
    else
        _numBoundHeads = _numHeads;


    // Initialize motor mechanics.
    //---------------------------------
    // Set stretching constant.
#ifdef PLOSFEEDBACK
    mMotorGhost_.setStretchingConstant(motorType, _numHeads);
#else
    mMotorGhost_.setStretchingConstant(motorType, _numBoundHeads);
#endif

    // Set equilibrium length.
    auto x1 = _c1->getFirstBead()->vcoordinate();
    auto x2 = _c1->getSecondBead()->vcoordinate();
    auto x3 = _c2->getFirstBead()->vcoordinate();
    auto x4 = _c2->getSecondBead()->vcoordinate();

    auto m1 = midPointCoordinate(x1, x2, c1->adjustedrelativeposition(position1));
    auto m2 = midPointCoordinate(x3, x4, c2->adjustedrelativeposition(position2));
    mMotorGhost_.eqLength = twoPointDistance(m1, m2);


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
    const auto filType1 = c1->getType();
    const auto filType2 = c2->getType();
          
    const int pos1 = int(position1 * SysParams::Geometry().cylinderNumMon[filType1]);
    const int pos2 = int(position2 * SysParams::Geometry().cylinderNumMon[filType2]);
          
    _cMotorGhost = make_unique<CMotorGhost>(motorSpeciesIndex1, motorSpeciesIndex2, _compartment, _c1->getCCylinder(), _c2->getCCylinder(), pos1, pos2);
    _cMotorGhost->setMotorGhost(this);

    
}

///@note - record lifetime data here
MotorGhost::~MotorGhost() noexcept {

//    floatingpoint lifetime = tau() - _birthTime;
//    
//    if(_lifetimes->getMax() > lifetime &&
//       _lifetimes->getMin() < lifetime)
//        _lifetimes->addValue(lifetime);
//        
//    
//    if(_walkLengths->getMax() > _walkLength &&
//       _walkLengths->getMin() < _walkLength)
//        _walkLengths->addValue(_walkLength);
    
}

void MotorGhost::updatePosition() {
    //update ccylinders
    _cMotorGhost->setFirstCCylinder(_c1->getCCylinder());
    _cMotorGhost->setSecondCCylinder(_c2->getCCylinder());
    
    //check if in same compartment
    updateCoordinate();
    Compartment* c;
    
    try {c = &GController::getCompartment(coordinate);}
    catch (exception& e) {
        cout << e.what();
        
        printSelf();
        
        exit(EXIT_FAILURE);
    }
    
    if(c != _compartment) {

        mins = chrono::high_resolution_clock::now();
        
        _compartment = c;

        SpeciesBound* firstSpecies = _cMotorGhost->getFirstSpecies();
        SpeciesBound* secondSpecies = _cMotorGhost->getSecondSpecies();
        CMotorGhost* clone = _cMotorGhost->clone(c);
        setCMotorGhost(clone);
        
        _cMotorGhost->setFirstSpecies(firstSpecies);
        _cMotorGhost->setSecondSpecies(secondSpecies);

        mine = chrono::high_resolution_clock::now();
        chrono::duration<floatingpoint> compartment_update(mine - mins);
        CUDAcommon::tmin.timemotorupdate += compartment_update.count();
        CUDAcommon::tmin.callsmotorupdate++;
    }
    
    if(SysParams::RUNSTATE) {
        auto x1 = _c1->getFirstBead()->vcoordinate();
        auto x2 = _c1->getSecondBead()->vcoordinate();
        auto x3 = _c2->getFirstBead()->vcoordinate();
        auto x4 = _c2->getSecondBead()->vcoordinate();

        auto m1 = midPointCoordinate(x1, x2, _c1->adjustedrelativeposition(_position1));
        auto m2 = midPointCoordinate(x3, x4, _c2->adjustedrelativeposition(_position2));

        //update the spring constant, based on numboundheads
        //current force
        floatingpoint force = max((floatingpoint) 0.0, mMotorGhost_.stretchForce);

        //update number of bound heads
        if (!_unbindingChangers.empty())
            _numBoundHeads = _unbindingChangers[_motorType]->numBoundHeads(_onRate, _offRate,
                                                                        force, _numHeads);
        else
            _numBoundHeads = _numHeads;

        #ifdef PLOSFEEDBACK
            mMotorGhost_.setStretchingConstant(_motorType, _numHeads);
        #else
            mMotorGhost_.setStretchingConstant(_motorType, _numBoundHeads);
        #endif
    }
}

/// @note - This function updates forward walking rates using the
/// stetching force in the opposite direction of the motor walk.
/// Does not consider negative forces in this direction.

/// Updates unbinding rates based on the stretch force. Does not
/// consider compression forces, only stretching.
void MotorGhost::updateReactionRates() {

    //current force
    floatingpoint force = max<floatingpoint>((floatingpoint)0.0,
            mMotorGhost_.stretchForce);
    
    //update number of bound heads
    if(!_unbindingChangers.empty())
        _numBoundHeads = _unbindingChangers[_motorType]->numBoundHeads(_onRate, _offRate, force, _numHeads);
    else
        _numBoundHeads = _numHeads;
    
    //walking rate changer
    if(!_walkingChangers.empty()) {
        auto x1 = _c1->getFirstBead()->vcoordinate();
        auto x2 = _c1->getSecondBead()->vcoordinate();
        auto x3 = _c2->getFirstBead()->vcoordinate();
        auto x4 = _c2->getSecondBead()->vcoordinate();

	    const auto& cylinderInfoData = Cylinder::getDbData();

	    auto c1struct = cylinderInfoData[_c1->getStableIndex()];
	    auto c2struct = cylinderInfoData[_c2->getStableIndex()];
	    auto fType1 = c1struct.type;
	    auto fType2 = c2struct.type;

	    bool consider_passivation = false;
	    bool isc1leftofc2 = false;
	    /* isc1leftofc2 tells which of the two cylinders is close to minus end. c1 is
	     * connected to leg1 and c2 is connected to leg2 of motor.
	     * If isc1leftofc2= true, in MotorWalkingForward reactions, the reaction
	     * involving leg1 should be passivated.
	     * If isc1leftofc2 = false, in MotorWalkingBackward reactions, the reaction
	     * involving leg2 should be passivated.*/

	    if((c1struct.filamentId == c2struct.filamentId)) {
	    	auto c1posonFil = c1struct.positionOnFilament;
		    auto c2posonFil = c2struct.positionOnFilament;
            consider_passivation = abs(c1posonFil - c2posonFil) <= ChemParams::minCylinderDistanceSameFilament + 1;
		    //A distance of 3 or lesser between two cylinders on the same filament is not
		    // acceptable.
		    isc1leftofc2 = c1posonFil < c2posonFil;
	    }

        auto mp1 = midPointCoordinate(x1, x2, _c1->adjustedrelativeposition(_position1));
        auto mp2 = midPointCoordinate(x3, x4, _c2->adjustedrelativeposition(_position2));

        //get component of force in direction of forward walk for C1, C2
        vector<floatingpoint> motorC1Direction =
        twoPointDirection(mp1, mp2);
        
        vector<floatingpoint> motorC2Direction =
        twoPointDirection(mp2, mp1);
        
        vector<floatingpoint> c1Direction = twoPointDirection(x2,x1);
        vector<floatingpoint> c2Direction = twoPointDirection(x4,x3);
        
        floatingpoint forceDotDirectionC1 = force * dotProduct(motorC1Direction, c1Direction);
        floatingpoint forceDotDirectionC2 = force * dotProduct(motorC2Direction, c2Direction);
        
        //WALKING REACTIONS
        Species* s1 = _cMotorGhost->getFirstSpecies();
        Species* s2 = _cMotorGhost->getSecondSpecies();
        
        for(auto r : s1->getRSpecies().reactantReactions()) {
            
            if(r->getReactionType() == ReactionType::MOTORWALKINGFORWARD) {
                if (SysParams::RUNSTATE == false)
                    r->setRateMulFactor(0.0f, ReactionBase::RESTARTPHASESWITCH);
                else
                    r->setRateMulFactor(1.0f, ReactionBase::RESTARTPHASESWITCH);

            	if(consider_passivation && isc1leftofc2){

		            float c1lastbindingsite = float(*(SysParams::Chemistry()
				            .bindingSites[fType1].end() -1))/float(SysParams::Geometry()
				            		.cylinderNumMon[fType1]);
		            if(areEqual(c1lastbindingsite,_position1))
			            r->setRateMulFactor(0.0f, ReactionBase::MOTORWALKCONSTRAINTFACTOR);
            	}
            	else{
		            r->setRateMulFactor(1.0f, ReactionBase::MOTORWALKCONSTRAINTFACTOR);
		            float newRate =
				            _walkingChangers[_motorType]->
						            changeRate(_cMotorGhost->getOnRate(),
						                       _cMotorGhost->getOffRate(),
						                       _numHeads, max<floatingpoint>((floatingpoint)0.0, forceDotDirectionC1));
/*		            if(SysParams::RUNSTATE==false){
			            newRate=0.0;}*/
#ifdef DETAILEDOUTPUT
                    std::cout<<"Motor WF1 f "<<force<<" Rate "<<newRate<<" "<<coordinate[0]<<" "
                        ""<<coordinate[1]<<" "<<coordinate[2]<<" Fdirn "<<
                         forceDotDirectionC2<<" NH "<<_numHeads<<endl;

#endif
		            r->setBareRate(newRate);
	            }/*if(consider_passivation && isc1leftofc2)*/
	            r->updatePropensity();

            }/*MOTORWALKINGFORWARD*/
            else if(r->getReactionType() == ReactionType::MOTORWALKINGBACKWARD) {
                if (SysParams::RUNSTATE == false) {
                    r->setRateMulFactor(0.0f, ReactionBase::RESTARTPHASESWITCH);
                }else {
                    r->setRateMulFactor(1.0f, ReactionBase::RESTARTPHASESWITCH);
                }
	            if(consider_passivation && !isc1leftofc2) {

		            auto c1firstbindingsite = float(*(SysParams::Chemistry()
				            .bindingSites[fType1].begin()))/
		                                      float(SysParams::Geometry()
		                                      .cylinderNumMon[fType1]);
		            if(areEqual(c1firstbindingsite,_position1))
			            r->setRateMulFactor(0.0f, ReactionBase::MOTORWALKCONSTRAINTFACTOR);

	            }else{
		            r->setRateMulFactor(1.0f, ReactionBase::MOTORWALKCONSTRAINTFACTOR);
		            float newRate =
				            _walkingChangers[_motorType]->
						            changeRate(_cMotorGhost->getOnRate(),
						                       _cMotorGhost->getOffRate(),
						                       _numHeads, max<floatingpoint>((floatingpoint)0.0, -forceDotDirectionC1));

/*		            if(SysParams::RUNSTATE==false){
		                newRate=0.0;}*/
#ifdef DETAILEDOUTPUT
		            std::cout<<"Motor WB1 f "<<force<<" Rate "<<newRate<<" "<<coordinate[0]<<" "
                        ""<<coordinate[1]<<" "<<coordinate[2]<<" Fdirn "<<
                         forceDotDirectionC2<<" NH "<<_numHeads<<endl;
#endif
		            r->setBareRate(newRate);
	            }
	            r->updatePropensity();
            }
        }
        for(auto r : s2->getRSpecies().reactantReactions()) {
            
            if(r->getReactionType() == ReactionType::MOTORWALKINGFORWARD) {
                if (SysParams::RUNSTATE == false) {
                    r->setRateMulFactor(0.0f, ReactionBase::RESTARTPHASESWITCH);
                }else {
                    r->setRateMulFactor(1.0f, ReactionBase::RESTARTPHASESWITCH);
                }
	            if(consider_passivation && !isc1leftofc2) {

		            float c2lastbindingsite = float(*(SysParams::Chemistry()
				            .bindingSites[fType2].end()-1))
				            		/float(SysParams::Geometry().cylinderNumMon[fType2]);
		            if(areEqual(c2lastbindingsite,_position2))
			            r->setRateMulFactor(0.0f, ReactionBase::MOTORWALKCONSTRAINTFACTOR);

	            }else {
		            r->setRateMulFactor(1.0f, ReactionBase::MOTORWALKCONSTRAINTFACTOR);
		            if(r->isPassivated()) {r->activateReaction();}
		            float newRate =
				            _walkingChangers[_motorType]->
						            changeRate(_cMotorGhost->getOnRate(),
						                       _cMotorGhost->getOffRate(),
						                       _numHeads,
						                       max<floatingpoint>((floatingpoint) 0.0,
						                                          forceDotDirectionC2));
/*		            if (SysParams::RUNSTATE == false) {
		                newRate = 0.0; }*/

		            r->setBareRate(newRate);
	            }
	            r->updatePropensity();

            }/*MOTORWALKINGFORWARD*/
            else if(r->getReactionType() == ReactionType::MOTORWALKINGBACKWARD) {
                if (SysParams::RUNSTATE == false) {
                    r->setRateMulFactor(0.0f, ReactionBase::RESTARTPHASESWITCH);
                }else {
                    r->setRateMulFactor(1.0f, ReactionBase::RESTARTPHASESWITCH);
                }
	            if(consider_passivation && isc1leftofc2) {
		            auto c2firstbindingsite = float(*(SysParams::Chemistry()
				            .bindingSites[fType2].begin()))/float(SysParams::Geometry()
				                                      .cylinderNumMon[fType2]);
		            if(areEqual(c2firstbindingsite,_position2))
			            r->setRateMulFactor(0.0f, ReactionBase::MOTORWALKCONSTRAINTFACTOR);
	            }else {
		            r->setRateMulFactor(1.0f, ReactionBase::MOTORWALKCONSTRAINTFACTOR);
	            	if(r->isPassivated()) {r->activateReaction();}
		            float newRate =
				            _walkingChangers[_motorType]->
						            changeRate(_cMotorGhost->getOnRate(),
						                       _cMotorGhost->getOffRate(),
						                       _numHeads, max<floatingpoint>((floatingpoint)0.0, -forceDotDirectionC2));
/*		            if(SysParams::RUNSTATE==false)
		            { newRate=0.0;}*/

		            r->setBareRate(newRate);
	            }
	            r->updatePropensity();

            }
        }
    }
    
    //unbinding rate changer
    if(!_unbindingChangers.empty()) {

        //get the unbinding reaction
        ReactionBase* offRxn = _cMotorGhost->getOffReaction();

        if (SysParams::RUNSTATE == false)
            offRxn->setRateMulFactor(0.0f, ReactionBase::RESTARTPHASESWITCH);
        else
            offRxn->setRateMulFactor(1.0f, ReactionBase::RESTARTPHASESWITCH);
        
        //change the rate
        float newRate =
        _unbindingChangers[_motorType]->
        changeRate(_cMotorGhost->getOnRate(), _cMotorGhost->getOffRate(), _numHeads, force);
#ifdef DETAILEDOUTPUT
        std::cout<<"Motor UB f "<<force<<" Rate "<<newRate<<" "<<coordinate[0]<<" "
                ""<<coordinate[1]<<" "
                ""<<coordinate[2]<<endl;
#endif
        offRxn->setBareRate(newRate);
        offRxn->activateReaction();
    }
}

void MotorGhost::moveMotorHead(Cylinder* c,
                               floatingpoint oldPosition, floatingpoint newPosition,
    int speciesMotorIndex, short boundType, SubSystem* ps
) {

    //shift the position of one side of the motor
    floatingpoint shift =  newPosition - oldPosition;
    
    //shift the head
    if(c == _c1) {
        _position1 += shift;
    }
    else {
        _position2 += shift;
    }
    short filType = c->getType();
    
    //record walk length
    _walkLength += shift * SysParams::Geometry().cylinderSize[filType];
    
    short oldpos = int (oldPosition * SysParams::Geometry().cylinderNumMon[filType]);
    short newpos = int (newPosition * SysParams::Geometry().cylinderNumMon[filType]);
    
    _cMotorGhost->moveMotorHead(c->getCCylinder(), oldpos, newpos,
                                speciesMotorIndex, boundType, ps);

    
}

void MotorGhost::moveMotorHead(Cylinder* oldC, Cylinder* newC,
                               floatingpoint oldPosition, floatingpoint newPosition,
    int speciesMotorIndex, short boundType, SubSystem* ps
) {
    //shift the head
    if(oldC == _c1) {
        _position1 = newPosition;
        _c1 = newC;
    }
    else {
        _position2 = newPosition;
        _c2 = newC;
    }
    const auto filType = newC->getType();
    
    //record walk length
    _walkLength += (1-oldPosition + newPosition) * SysParams::Geometry().cylinderSize[filType];
    
    short oldpos = int (oldPosition * SysParams::Geometry().cylinderNumMon[filType]);
    short newpos = int (newPosition * SysParams::Geometry().cylinderNumMon[filType]);
    
    _cMotorGhost->moveMotorHead(oldC->getCCylinder(), newC->getCCylinder(),
                                oldpos, newpos, speciesMotorIndex, boundType, ps);
}


void MotorGhost::printSelf()const {
    
    cout << endl;
    
    cout << "MotorGhost: ptr = " << this << endl;
    cout << "Motor type = " << _motorType << ", Motor ID = " << getId() << endl;
    cout << "Coordinates = " << coordinate[0] << ", " << coordinate[1] << ", " << coordinate[2] << endl;
    
    cout << "Position on first cylinder (floatingpoint) = " << _position1 << endl;
    cout << "Position on second cylinder (floatingpoint) = " << _position2 << endl;
    
    cout << "Number of heads = " << _numHeads << endl;
    cout << "Birth time = " << _birthTime << endl;
    
    cout << endl;
    
    cout << "Associated species 1 = " << _cMotorGhost->getFirstSpecies()->getName()
    << " , copy number = " << _cMotorGhost->getFirstSpecies()->getN()
    << " , position on first cylinder (int) = " << _cMotorGhost->getFirstPosition() << endl;
    
    cout << "Associated species 2 = " << _cMotorGhost->getSecondSpecies()->getName()
    << " , copy number = " << _cMotorGhost->getSecondSpecies()->getN()
    << " , position on second cylinder (int) = " << _cMotorGhost->getSecondPosition() << endl;
    
    cout << endl;
    
    cout << "Associated cylinders (one and two): " << endl;
    _c1->printSelf();
    _c2->printSelf();
    
    cout << endl;
}

species_copy_t MotorGhost::countSpecies(const string& name) {
    
    species_copy_t copyNum = 0;
    
    for(auto m : getElements()) {
        
        auto s = m->getCMotorGhost()->getFirstSpecies();
        string sname = SpeciesNamesDB::removeUniqueFilName(s->getName());
        
        if(sname == name)
            copyNum += s->getN();
    }
    return copyNum;
}

vector<MotorRateChanger*> MotorGhost::_unbindingChangers;
vector<MotorRateChanger*> MotorGhost::_walkingChangers;

Histogram* MotorGhost::_lifetimes;
Histogram* MotorGhost::_walkLengths;

} // namespace medyan

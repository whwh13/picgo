
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

#include "Structure/Cylinder.h"

#include "SubSystem.h"
#include "Controller/CController.h"
#include "ChemManager.h"
#include "ChemRNode.h"

#include "Filament.h"
#include "Bead.h"

#include "Controller/GController.h"
#include "MathFunctions.h"

namespace medyan {
using namespace mathfunc;

void Cylinder::updateData() {
    auto& data = getDbData()[getStableIndex()];

    data.filamentId = static_cast<Filament*>(getParent())->getId();
    data.positionOnFilament = _position;
    data.compartmentId = getCompartment()->getId();
    data.beadIndices[0] = _b1->getStableIndex();
    data.beadIndices[1] = _b2->getStableIndex();
    data.coord = vector2Vec<3, floatingpoint>(coordinate);
    data.type = getType();
    data.id = getId();

    data.chemCylinder = getCCylinder();
}

void Cylinder::updateCoordinate() {
    coordinate = midPointCoordinate(_b1->vcoordinate(), _b2->vcoordinate(), 0.5);
    //update the coordiante in cylinder structure.
    Cylinder::getDbData()[getStableIndex()].coord = 0.5 * (_b1->coordinate() + _b2->coordinate());
}

Cylinder::Cylinder(Composite* parent, Bead* b1, Bead* b2, short type, int position,
                   bool extensionFront, bool extensionBack, bool initialization,
                   floatingpoint eqLength)

        : Trackable(true, true, true, false),
          _b1(b1), _b2(b2), _type(type), _position(position)
{

    #if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
    if(getStableIndex() >= SysParams::Chemistry().maxStableIndex) {

        LOG(ERROR) << "Total number of cylinders initialized("<< getStableIndex()<<
                    ") equals/exceeds the maximum ("<<SysParams::Chemistry()
                    .maxStableIndex<<")."
                    "Check number of binding sites to continue to use "
                    "HYBRID_NLSTENCILLIST or SIMDBINDINGSEARCH. If not, shift to "
                    "other binding search algorithms. Exiting.";
        throw std::logic_error("Max value reached");
    }
    #endif

    parent->addChild(unique_ptr<Component>(this));
    //setID
    _filID = static_cast<Filament*>(parent)->getId();
    //@{

    //Set coordinate
    updateCoordinate();

    Compartment* compartment;
    try {compartment = &GController::getCompartment(coordinate);}
    catch (exception& e) {
        cout << e.what() << endl;
        throw;
    }

    //add to compartment
    _cellElement.manager = compartment->cylinderCell.manager;
    _cellElement.manager->addElement(this, _cellElement, compartment->cylinderCell);

    //@}

    _cCylinder = unique_ptr<CCylinder>(new CCylinder(compartment, this));
    _cCylinder->setCylinder(this);

    if(SysParams::RUNSTATE) {
        eqLength = twoPointDistance(b1->vcoordinate(), b2->vcoordinate());

        //init using chem manager
        _chemManager->initializeCCylinder(_cCylinder.get(), extensionFront,
                                          extensionBack, initialization);
    }

    _mCylinder = unique_ptr<MCylinder>(new MCylinder(_type, eqLength));
    _mCylinder->setCylinder(this);

    // Update the stored data
    updateData();
    #ifdef CROSSCHECK_CYLINDER
    cout<<"Cylinder created "<<getId()<<" "<<getStableIndex()<<endl;
    #endif
}

void Cylinder::initializerestart(int nummonomers, int firstmonomer, int lastmonomer, bool
                                    minusendstatus, bool plusendstatus, short
                                    minusendtype, short plusendtype){
    if(SysParams::RUNSTATE){
        LOG(ERROR) << "initializerestart Function from Cylinder class can only be called "
                      "during restart phase. Exiting.";
        throw std::logic_error("Illegal function call pattern");
    }
    _chemManager->initializeCCylinder(_cCylinder.get(), false,
                                      false, true, nummonomers, firstmonomer,
                                      lastmonomer, minusendstatus, plusendstatus,
                                      minusendtype, plusendtype);
}

Cylinder::~Cylinder() noexcept {
    #ifdef CROSSCHECK_CYLINDER
    cout<<"Cylinder deleting "<<getId()<<" "<<getStableIndex()<<endl;
    #endif

    //remove from compartment
    _cellElement.manager->removeElement(_cellElement);
    #ifdef CROSSCHECK_CYLINDER
    cout<<"Cylinder deleted "<<getId()<<endl;
    #endif
}

/// Get filament type
int Cylinder::getType() {return _type;}

void Cylinder::updatePosition() {
    if(!setpositionupdatedstate) {

        //check if Cylinder is still in same compartment, set new position
        updateCoordinate();
        #ifdef CROSSCHECK_CYLINDER
        _crosscheckdumpFile <<"Coord updated "<<getId()<<endl;
        #endif
        Compartment *c;
        try { c = &GController::getCompartment(coordinate); }
        catch (exception &e) {
            cout << e.what();

            printSelf();

            throw;
        }

        Compartment* curCompartment = getCompartment();
        if (c != curCompartment) {
            #ifdef CROSSCHECK_CYLINDER
            _crosscheckdumpFile <<"Attempt Move cmp "<<getId()<<endl;
            #endif

            mins = chrono::high_resolution_clock::now();

            //remove from old compartment, add to new
            _cellElement.manager->updateElement(_cellElement, c->cylinderCell);
            #ifdef CROSSCHECK_CYLINDER
            _crosscheckdumpFile <<"Complete Move cmp "<<getId()<<endl;
            #endif

//			auto oldCCylinder = _cCylinder.get();

            //Remove old ccylinder from binding managers
            //Removed March 8, 2019 Aravind. Unnecessary as all UpdatePosition calls are
            // immediately followed by UpdateNeighborLists call in Controller.cpp/.cu
/*        for(auto &manager : oldCompartment->getFilamentBindingManagers()) {
#ifdef NLORIGINAL
            manager->removePossibleBindings(oldCCylinder);
#endif
#ifdef NLSTENCILLIST
            manager->removePossibleBindingsstencil(oldCCylinder);
#endif
        }*/

            //clone and set new ccylinder
            CCylinder *clone = _cCylinder->clone(c);
            setCCylinder(clone);
#ifdef CROSSCHECK_CYLINDER
            _crosscheckdumpFile <<"Clone CCyl "<<getId()<<endl;
#endif

//			auto newCCylinder = _cCylinder.get();

            //change both CCylinder and Compartment ID in the vector
            auto& data = getDbData()[getStableIndex()];
            data.compartmentId = c->getId();
            data.chemCylinder = _cCylinder.get();
#ifdef CROSSCHECK_CYLINDER
            _crosscheckdumpFile <<"Update CylinderData "<<getId()<<endl;
#endif

            mine = chrono::high_resolution_clock::now();
            chrono::duration<floatingpoint> compartment_update(mine - mins);
            CUDAcommon::tmin.timecylinderupdate += compartment_update.count();
            CUDAcommon::tmin.callscylinderupdate++;

        }

        //update length
        _mCylinder->setLength(twoPointDistance(_b1->vcoordinate(), _b2->vcoordinate()));

        #ifdef CROSSCHECK_CYLINDER
        _crosscheckdumpFile <<"MCylinder updated "<<getId()<<endl;
        #endif
    }
}

/// @note -  The function uses the bead load force to calculate this changed rate.
/// If there is no force on the beads the reaction rates are set to the bare.

void Cylinder::updateReactionRates() {

    floatingpoint force;

    //if no rate changer was defined, skip
    if(_polyChanger.empty()) return;

    //load force from front (affects plus end polymerization)
    if(_plusEnd) {

        //get force of front bead
        force = _b2->getLoadForcesP();

        //change all plus end polymerization rates
        for(auto &r : _cCylinder->getInternalReactions()) {
            floatingpoint factor;
            if(r->getReactionType() == ReactionType::POLYMERIZATIONPLUSEND) {

                //If reaching a threshold time for manual treadmilling rate changer
                if(tau() > SysParams::DRParams.manualCharStartTime){
                    //all bare rate will be change by a threshold ratio
                    factor = _polyChanger[_type]->getRateChangeFactor(force)*
                             SysParams::DRParams.manualPlusPolyRate;
                }
                else{
                    factor = _polyChanger[_type]->getRateChangeFactor(force);
                }

                r->setRateMulFactor(factor, ReactionBase::mechanochemical);
                r->updatePropensity();

            }

            //change all plus end depolymerization rates, not force dependent
            //If reaching a threshold time for manual treadmilling rate changer
            if(tau() > SysParams::DRParams.manualCharStartTime){
                if(r->getReactionType() == ReactionType::DEPOLYMERIZATIONPLUSEND) {
                    r->setRateMulFactor(SysParams::DRParams.manualPlusDepolyRate,
                                        ReactionBase::MANUALRATECHANGEFACTOR1);
                    r->updatePropensity();
                }
            }
        }
    }

    //load force from back (affects minus end polymerization)
    if(_minusEnd) {

        //get force of front bead
        force = _b1->getLoadForcesM();

        //change all plus end polymerization rates
        for(auto &r : _cCylinder->getInternalReactions()) {
            floatingpoint factor;
            if(r->getReactionType() == ReactionType::POLYMERIZATIONMINUSEND) {

                //If reaching a threshold time for manual treadmilling rate changer
                if(tau() > SysParams::DRParams.manualCharStartTime){
                    //all bare rate will be change by a threshold ratio
                    factor = _polyChanger[_type]->getRateChangeFactor(force) *
                            SysParams::DRParams.manualMinusPolyRate;
                }
                else{

                    factor = _polyChanger[_type]->getRateChangeFactor(force);
                }

                r->setRateMulFactor(factor, ReactionBase::mechanochemical);
                r->updatePropensity();
            }

            //change all minus end depolymerization rates, not force dependent
            //If reaching a threshold time for manual treadmilling rate changer
            if(tau() > SysParams::DRParams.manualCharStartTime){

                if(r->getReactionType() == ReactionType::DEPOLYMERIZATIONMINUSEND) {
                    r->setRateMulFactor(SysParams::DRParams.manualPlusDepolyRate,
                                        ReactionBase::MANUALRATECHANGEFACTOR1);
                    r->updatePropensity();
                }
            }
        }
    }
}

bool Cylinder::isFullLength() {

    return areEqual(_mCylinder->getEqLength(), SysParams::Geometry().cylinderSize[_type]);
}

void Cylinder::printSelf()const {

    cout << endl;

    cout << "Cylinder: ptr = " << this << endl;
    cout << "Cylinder ID = " << getId() << endl;
    cout << "Stable Index = "<< getStableIndex() << endl;
    cout << "Parent ptr = " << getParent() << endl;
    cout << "Coordinates = " << coordinate[0] << ", " << coordinate[1] << ", " << coordinate[2] << endl;

    if(_plusEnd) cout << "Is a plus end." << endl;
    if(_minusEnd) cout << "Is a minus end." << endl;

    if(_branchingCylinder != nullptr) cout << "Has a branching cylinder." << endl;

    cout << "Position = " << _position << endl;

    cout<< "Length "<<_mCylinder->getLength()<<endl;
    cout<< "Eq Length "<<_mCylinder->getEqLength()<<endl;
    cout<< "Eq Theta "<<_mCylinder->getEqTheta()<<endl;
    cout<<" Stretching constant "<<_mCylinder->getStretchingConst()<<endl;
    cout<<" Bending constant "<<_mCylinder->getBendingConst()<<endl;

    cout << endl;

    cout << "Chemical composition of cylinder:" << endl;
    _cCylinder->printCCylinder();

    cout << endl;

    cout << "Bead information..." << endl;

    _b1->printSelf();
    _b2->printSelf();

    cout << endl;
}
//Ask Qin when this is used
bool Cylinder::within(Cylinder* other, floatingpoint dist) {

    //check midpoints
    if(twoPointDistancesquared(coordinate, other->coordinate) <= (dist * dist))
        return true;

    //briefly check endpoints of other
    if(twoPointDistancesquared(coordinate, other->_b1->vcoordinate()) <= (dist * dist) ||
       twoPointDistancesquared(coordinate, other->_b2->vcoordinate()) <= (dist * dist))
        return true;

    return false;
}

//adjust the position variable according to the length of cylinder
//Refer Docs/Design/PartialCylinderAlpha.pdf
floatingpoint Cylinder::adjustedrelativeposition(floatingpoint _alpha, bool verbose){
#ifdef ADJUSTRELPOS
    //Full Length Cylinder
    if(isFullLength())
        return _alpha;
    floatingpoint _alphacorr = (floatingpoint)0.0;
    auto x1 = _b1->vcoordinate();
    auto x2 = _b2->vcoordinate();
    floatingpoint L = twoPointDistance(x1, x2);
    short filamentType = _type;
    short minusendmonomer = 0;
    floatingpoint fullcylinderSize = SysParams::Geometry().cylinderSize[filamentType];
    //Partial Plus End cylinder
    if(_plusEnd == true){
        floatingpoint Lm = (floatingpoint)0.0;//Distance of minus end from 0th monomer.
        // Both Minus and Plus End at the same time (Filament is one cylinder long)
        if(_minusEnd == true){
            int numMonomers = SysParams::Geometry().cylinderNumMon[filamentType];
            auto monomersize = SysParams::Geometry().monomerSize[filamentType];
            for(int midx = 0; midx<numMonomers; midx++){
                short m = _cCylinder->getCMonomer(midx)->activeSpeciesMinusEnd();
                short p = _cCylinder->getCMonomer(midx)->activeSpeciesPlusEnd();
                if(m != -1) {
                    minusendmonomer = midx;
                    break;
                }
            }
            Lm = minusendmonomer*monomersize;//Distance of minus end from 0th monomer.
        }
        _alphacorr = (fullcylinderSize*_alpha - Lm)/L;
    }
    //Parital Minus End Cylinder
    else if(_minusEnd == true){
        _alphacorr = 1-(1-_alpha)*fullcylinderSize/L;
    }

    if(verbose){
        cout<<"cyl ID "<<getId()<<" minsendstatus "<<_minusEnd<<" <<plusendstatus "
            <<_plusEnd<<" minusendmonomer "<<minusendmonomer<<" alpha "<<_alpha
            <<" _alphacorr "<<_alphacorr<<endl;
    }

    if(_alphacorr < (floatingpoint)0.0)
        return (floatingpoint)0.0;
    else if(_alphacorr > (floatingpoint)1.0)
        return (floatingpoint)1.0;
    else
        return _alphacorr;
#else
    return _alpha;
#endif
}

vector<FilamentRateChanger*> Cylinder::_polyChanger;

bool Cylinder::setpositionupdatedstate = false;
floatingpoint Cylinder::timecylinder1 = 0.0;
floatingpoint Cylinder::timecylinder2= 0.0;
floatingpoint Cylinder::timecylinderchem= 0.0;
floatingpoint Cylinder::timecylindermech= 0.0;
ofstream Cylinder::_crosscheckdumpFile;

} // namespace medyan

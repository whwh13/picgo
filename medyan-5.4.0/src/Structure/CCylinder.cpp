
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

#include "CCylinder.h"

#include "Cylinder.h"

#include "CBound.h"
#include "ChemManager.h"
#include "Structure/Compartment.h"

namespace medyan {
/// Default constructor, sets compartment and cylinder
CCylinder::CCylinder(Compartment* C, Cylinder* c)
    : _compartment(C), _pCylinder(c) {
    //set size based on parent cylinder
    _size = SysParams::Geometry().cylinderSize[c->getType()] /
            SysParams::Geometry().monomerSize[c->getType()];
}


CCylinder::CCylinder(const CCylinder& rhs, Compartment* c)
    : _compartment(c), _pCylinder(rhs._pCylinder), _size(rhs._size) {
        
    CCylinder* rhsPtr = const_cast<CCylinder*>(&rhs);
    #ifdef CROSSCHECK_CYLINDER
    Cylinder::_crosscheckdumpFile <<"Clone begin rhs cylinder set "<<getId()<<endl;
    Cylinder::_crosscheckdumpFile <<"Total monomers "<<rhs._monomers.size()<<getId()<<endl;
    #endif
    //copy all monomers, bounds
    for(auto &m : rhs._monomers)
        _monomers.emplace_back(m->clone(c));

    #ifdef CROSSCHECK_CYLINDER
    Cylinder::_crosscheckdumpFile <<"Monomers cloned"<<endl;
    #endif
    
    //copy all internal reactions
    #ifdef CROSSCHECK_CYLINDER
    Cylinder::_crosscheckdumpFile <<"Internal rxn count "<<rhs
        ._internalReactions.size()<<endl;
    #endif
#ifdef OPTIMOUT
    mins = chrono::high_resolution_clock::now();
#endif
    unsigned int count = 0;
    for(auto &r: rhs._internalReactions) {
        #ifdef CROSSCHECK_CYLINDER
        Cylinder::_crosscheckdumpFile <<"Internal ReactionType "<<r->getReactionType()
        <<endl;
        if(r->getCBound() != nullptr)
            Cylinder::_crosscheckdumpFile <<"CBound exists"<<endl;
        else
            Cylinder::_crosscheckdumpFile <<"CBound DOESNOT exist"<<endl;
        #endif
#ifdef OPTIMOUT
        minsi = chrono::high_resolution_clock::now();
#endif
        ReactionBase* rxnClone = r->clone(c->getSpeciesContainer());
#ifdef OPTIMOUT
        minei = chrono::high_resolution_clock::now();
        chrono::duration<floatingpoint> irxnclone(minei - minsi);
        CUDAcommon::cdetails.internalrxnclone += irxnclone.count();
#endif
        rxnClone->setVolumeFrac(c->getVolumeFrac());
        #ifdef CROSSCHECK_CYLINDER
        Cylinder::_crosscheckdumpFile <<"Internal Reaction cloned."<<endl;
        #endif
        
        if(r->getCBound() != nullptr)
            r->getCBound()->setOffReaction(rxnClone);

        #ifdef CROSSCHECK_CYLINDER
        if(r->getCBound() != nullptr)
            Cylinder::_crosscheckdumpFile <<"OffReaction cloned + set."<<endl;
        #endif
#ifdef OPTIMOUT
        minsi = chrono::high_resolution_clock::now();
#endif
        addInternalReaction(rxnClone);
#ifdef OPTIMOUT
        minei = chrono::high_resolution_clock::now();
        chrono::duration<floatingpoint> irxnadd(minei - minsi);
        CUDAcommon::cdetails.internalrxnadd += irxnadd.count();
#endif
        #ifdef CROSSCHECK_CYLINDER
        Cylinder::_crosscheckdumpFile <<"Internal ReactionAdded."<<endl;
        #endif
        count++;
    }
#ifdef OPTIMOUT
    mine = chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> internalrxn(mine - mins);
    CUDAcommon::cdetails.ccylclonecounter[0]++;
    CUDAcommon::cdetails.ccylclonetimer[0] += internalrxn.count();
    CUDAcommon::cdetails.ccylclonerxncounter[0] += rhs._internalReactions.size();
    mins = chrono::high_resolution_clock::now();
#endif
    //copy all cross-cylinder reactions
    #ifdef CROSSCHECK_CYLINDER
    Cylinder::_crosscheckdumpFile <<"Cross cylinder rxn count "<<rhs
    ._crossCylinderReactions.size()<<endl;
    #endif
    for(auto it = rhs._crossCylinderReactions.begin();
             it != rhs._crossCylinderReactions.end(); it++) {
        for(auto &r : it->second) {
            #ifdef CROSSCHECK_CYLINDER
            Cylinder::_crosscheckdumpFile <<"CrossCylinder ReactionType "
                                            ""<<r->getReactionType()
                                          <<endl;
            if(r->getCBound() != nullptr)
                Cylinder::_crosscheckdumpFile <<"CBound exists"<<endl;
            else
                Cylinder::_crosscheckdumpFile <<"CBound DOESNOT exist"<<endl;
            #endif

            //copy cbound if any
            ReactionBase* rxnClone = r->clone(c->getSpeciesContainer());
            rxnClone->setVolumeFrac(c->getVolumeFrac());
            #ifdef CROSSCHECK_CYLINDER
            Cylinder::_crosscheckdumpFile <<"CrossCylinder Reaction cloned."<<endl;
            #endif
            
            if(r->getCBound() != nullptr) {
	            r->getCBound()->setOffReaction(rxnClone);
            }
            #ifdef CROSSCHECK_CYLINDER
            if(r->getCBound() != nullptr)
                Cylinder::_crosscheckdumpFile <<"OffReaction cloned + set."<<endl;
            #endif
            
            addCrossCylinderReaction(it->first, rxnClone);
            #ifdef CROSSCHECK_CYLINDER
            Cylinder::_crosscheckdumpFile <<"CrossCylinder ReactionAdded."<<endl;
            #endif
        }
    }
    //Copy reacting cylinders, Clone reactions where this cylinder is involved
    #ifdef CROSSCHECK_CYLINDER
    Cylinder::_crosscheckdumpFile <<"Reacting cylinder count "<<rhs
        ._reactingCylinders.size()<<endl;
    #endif
#ifdef OPTIMOUT
    mine = chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> crosscyl(mine - mins);
    CUDAcommon::cdetails.ccylclonecounter[1]++;
    CUDAcommon::cdetails.ccylclonetimer[1] += crosscyl.count();
    CUDAcommon::cdetails.ccylclonerxncounter[1] += rhs._crossCylinderReactions.size();
    mins = chrono::high_resolution_clock::now();
#endif
    for(auto &ccyl : rhs._reactingCylinders) {

        //clone reactions
        for(auto &r: ccyl->getCrossCylinderReactions()[rhsPtr]) {
            #ifdef CROSSCHECK_CYLINDER
            Cylinder::_crosscheckdumpFile <<"Reacting Cylinder ReactionType "
                                            ""<<r->getReactionType()
                                          <<endl;
            if(r->getCBound() != nullptr)
                Cylinder::_crosscheckdumpFile <<"CBound exists"<<endl;
            else
                Cylinder::_crosscheckdumpFile <<"CBound DOESNOT exist"<<endl;
            #endif
            
            //copy cbound if any
            ReactionBase* rxnClone = r->clone(c->getSpeciesContainer());
            rxnClone->setVolumeFrac(c->getVolumeFrac());
            #ifdef CROSSCHECK_CYLINDER
            Cylinder::_crosscheckdumpFile <<"ReactingCylinder Reaction cloned."<<endl;
            #endif

            if(r->getCBound() != nullptr) {
                r->getCBound()->setOffReaction(rxnClone);
            }
            #ifdef CROSSCHECK_CYLINDER
            if(r->getCBound() != nullptr)
                Cylinder::_crosscheckdumpFile <<"OffReaction cloned + set."<<endl;
            #endif
            
            ccyl->addCrossCylinderReaction(this, rxnClone);
            #ifdef CROSSCHECK_CYLINDER
            Cylinder::_crosscheckdumpFile <<"CrossCylinder ReactionAdded."<<endl;
            #endif
        }
    }
#ifdef OPTIMOUT
    mine = chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> reactingcyl(mine - mins);
    CUDAcommon::cdetails.ccylclonecounter[2]++;
    CUDAcommon::cdetails.ccylclonetimer[2] += reactingcyl.count();
    CUDAcommon::cdetails.ccylclonerxncounter[2] += rhs._reactingCylinders.size();
#endif
}

void CCylinder::addInternalReaction(ReactionBase* r) {
    
    //add to compartment and chemsim
    _compartment->addInternalReaction(r);
    _chemSim->addReaction(r);
    
    //add to local reaction list
    _internalReactions.insert(r);
    
    //activate reaction
    r->activateReaction();
}


void CCylinder::removeInternalReaction(ReactionBase* r) {
    
    //remove from internal reaction list
    if (_internalReactions.find(r) != _internalReactions.end()) {
//        std::cout<<"passivate removeInternalReaction"<<endl;
        //passivate reaction, removing from dependents
        r->passivateReaction();
        
        //remove from compartment and chemsim
        _chemSim->removeReaction(r);
        _compartment->removeInternalReaction(r);
        
        _internalReactions.erase(r);
    }
}

void CCylinder::addCrossCylinderReaction(CCylinder* other,
                                         ReactionBase* r) {
    
    //add to compartment and chemsim
    _compartment->addInternalReaction(r);
    _chemSim->addReaction(r);

    //add to this reaction map
    _crossCylinderReactions[other].insert(r);
    other->addReactingCylinder(this);
    
    //activate reaction
    r->activateReaction();
}

void CCylinder::addReactingCylinder(CCylinder* other) {
    _reactingCylinders.insert(other);
}

void CCylinder:: removeAllInternalReactions() {
    
    auto tempReactions = _internalReactions;
    for (auto &r : tempReactions) removeInternalReaction(r);
}

void CCylinder::removeCrossCylinderReaction(CCylinder* other,
                                            ReactionBase* r) {
	if (r == nullptr) return;
    auto it = _crossCylinderReactions[other].find(r);
    if(it != _crossCylinderReactions[other].end()) {
       
        //erase the reaction
        _crossCylinderReactions[other].erase(it);
//        std::cout<<"passivate removeCrossCylinderReaction"<<endl;

        //passivate reaction, removing from dependents
        r->passivateReaction();
        
        //remove from compartment and chemsim
        _chemSim->removeReaction(r);
        _compartment->removeInternalReaction(r);
        
        //if number of reactions in cross-cylinder
        //has dropped to zero, delete it
        if(_crossCylinderReactions[other].empty()) {
            
            _crossCylinderReactions.erase(other);
            
            //also remove from reacting of other ccylinder
            auto it2 =other->_reactingCylinders.find(this);
            
            if(it2 != other->_reactingCylinders.end())
                other->_reactingCylinders.erase(it2);
        }
    }
}

void CCylinder::removeCrossCylinderReactions(CCylinder* other) {
    
    auto tempReactions = _crossCylinderReactions[other];
    
    for(auto &r : tempReactions)
        removeCrossCylinderReaction(other, r);
}

void CCylinder::removeAllCrossCylinderReactions() {
    
    auto tempMap = _crossCylinderReactions;
    
    for(auto it = tempMap.begin(); it != tempMap.end(); it++)
        removeCrossCylinderReactions(it->first);
}

void CCylinder::removeReactingCylinder(CCylinder* other) {
    
    other->removeCrossCylinderReactions(this);
}

void CCylinder::removeAllReactingCylinders() {
    
    auto tempReactingCylinders = _reactingCylinders;
    
    for(auto &cc : tempReactingCylinders)
        cc->removeCrossCylinderReactions(this);
}

CCylinder::~CCylinder() {
    #ifdef CROSSCHECK_CYLINDER
    cout<<"CCylinder deleting "<<endl;
    cout<<"Total monomers "<<_monomers.size()<<endl;
    #endif
    
    //Remove all reactions owned by this ccylinder
    removeAllInternalReactions();
    #ifdef CROSSCHECK_CYLINDER
    cout<<"CCylinder removed all InternalReactions "<<endl;
    #endif
    removeAllCrossCylinderReactions();
    #ifdef CROSSCHECK_CYLINDER
    cout<<"CCylinder removed all CrossCylinderReactions "<<endl;
    #endif
    
    //remove all reactions involving this ccylinder
    removeAllReactingCylinders();
    #ifdef CROSSCHECK_CYLINDER
    cout<<"CCylinder removed all ReactingCylinders "<<endl;
    #endif
    
    //Remove all species
    for(auto &m: _monomers) {
        #ifdef CROSSCHECK_CYLINDER
        cout<<"Total Filament Species #"<<(CMonomer::_numFSpecies[_pCylinder->getType()])<<endl;
        #endif
        for(int i = 0; i < CMonomer::_numFSpecies[_pCylinder->getType()]; i++) {
            Species* s = m->speciesFilament(i);
            #ifdef CROSSCHECK_CYLINDER
            cout<<"Removing Species Filament "<<s->getFullName()<<endl;
            #endif
            if(s != nullptr) _compartment->removeSpecies(s);
        }
        #ifdef CROSSCHECK_CYLINDER
        cout<<"CCylinder removed all FilamentSpecies "<<endl;
        cout<<"Total Bound Species #"<<CMonomer::_numBSpecies[_pCylinder->getType()]<<endl;
        #endif
        for(int i = 0; i < CMonomer::_numBSpecies[_pCylinder->getType()]; i++) {
            SpeciesBound* s = m->speciesBound(i);
            #ifdef CROSSCHECK_CYLINDER
            cout<<"Removing Species Bound "<<s->getFullName()<<endl;
            #endif
            if(s != nullptr) _compartment->removeSpecies(s);
        }
        #ifdef CROSSCHECK_CYLINDER
        cout<<"CCylinder removed all BoundSpecies "<<endl;
        #endif
    }
}

void CCylinder::passivatefilcrossreactions(){
    
    for (auto it2=_crossCylinderReactions.begin(); it2!=_crossCylinderReactions.end(); ++it2){
        auto mySet = it2->second;
        for (auto it: mySet) {
            if(it->getReactionType() ==ReactionType::POLYMERIZATIONMINUSEND
               ||it->getReactionType() ==ReactionType::POLYMERIZATIONPLUSEND
               ||it->getReactionType() ==ReactionType::DEPOLYMERIZATIONMINUSEND
               ||it->getReactionType() ==ReactionType::DEPOLYMERIZATIONPLUSEND
               ||it->getReactionType() ==ReactionType::SEVERING
               ||it->getReactionType() ==ReactionType::FILAMENTDESTRUCTION
               ||it->getReactionType() ==ReactionType::AGING)
            {it->passivateReaction();}
        
        }}
//    auto tempReactions = _crossCylinderReactions[this];
//    if(this->getCylinder()->isPlusEnd()){
//        for(auto &it : tempReactions){
//            std::cout<<it->getReactionType()<<endl;
//        }
//    }
//    for(auto &it : tempReactions){
//        if(it->getReactionType() ==ReactionType::POLYMERIZATIONMINUSEND
//           ||it->getReactionType() ==ReactionType::POLYMERIZATIONPLUSEND
//           ||it->getReactionType() ==ReactionType::DEPOLYMERIZATIONMINUSEND
//           ||it->getReactionType() ==ReactionType::DEPOLYMERIZATIONPLUSEND
//           ||it->getReactionType() ==ReactionType::SEVERING
//           ||it->getReactionType() ==ReactionType::FILAMENTDESTRUCTION
//           ||it->getReactionType() ==ReactionType::AGING)
//        {it->passivateReaction();}
//    }
}

void CCylinder::activatefilcrossreactions(){    
    for (auto it2=_crossCylinderReactions.begin(); it2!=_crossCylinderReactions.end(); ++it2){
        auto mySet = it2->second;
        for (auto it: mySet) {
            if(it->getReactionType() ==ReactionType::POLYMERIZATIONMINUSEND
               ||it->getReactionType() ==ReactionType::POLYMERIZATIONPLUSEND
               ||it->getReactionType() ==ReactionType::DEPOLYMERIZATIONMINUSEND
               ||it->getReactionType() ==ReactionType::DEPOLYMERIZATIONPLUSEND
               ||it->getReactionType() ==ReactionType::SEVERING
               ||it->getReactionType() ==ReactionType::FILAMENTDESTRUCTION
               ||it->getReactionType() ==ReactionType::AGING)
            {it->activateReaction();}
            
        }}}
void CCylinder::passivatefilreactions(){
    for(auto &it: _internalReactions){
        if(it->getReactionType() ==ReactionType::POLYMERIZATIONMINUSEND
           ||it->getReactionType() ==ReactionType::POLYMERIZATIONPLUSEND
           ||it->getReactionType() ==ReactionType::DEPOLYMERIZATIONMINUSEND
           ||it->getReactionType() ==ReactionType::DEPOLYMERIZATIONPLUSEND
           ||it->getReactionType() ==ReactionType::SEVERING
           ||it->getReactionType() ==ReactionType::FILAMENTDESTRUCTION
           ||it->getReactionType() ==ReactionType::AGING)
        {it->passivateReaction();}}}
void CCylinder::activatefilreactions(){
    for(auto &it: _internalReactions){
        if(it->getReactionType() ==ReactionType::POLYMERIZATIONPLUSEND
           ||it->getReactionType() ==ReactionType::POLYMERIZATIONMINUSEND
           ||it->getReactionType() ==ReactionType::DEPOLYMERIZATIONMINUSEND
           ||it->getReactionType() ==ReactionType::DEPOLYMERIZATIONPLUSEND
           ||it->getReactionType() ==ReactionType::SEVERING
           ||it->getReactionType() ==ReactionType::FILAMENTDESTRUCTION
           ||it->getReactionType() ==ReactionType::AGING)
        {it->activateReaction();}}}

vector<ReactionBase*> CCylinder::getAllReactions() {
    
    vector<ReactionBase*> reactions;
    
    //get internal rxns
    for(auto r : _internalReactions) reactions.push_back(r);
    
    //get cross cylinder rxns
    for(auto it = _crossCylinderReactions.begin();
            it != _crossCylinderReactions.end(); it++)
        
        for(auto it2 = it->second.begin();
                it2 != it->second.end(); it2++)
            
            reactions.push_back(*it2);
    
    return reactions;
}

void CCylinder::printCCylinder()
{
    cout << "Compartment:" << _compartment << endl;
    cout << "CCylinder: ptr = " << this << endl;
    
    cout << "Composition of CCylinder: " << endl;
    for (auto &m : _monomers){
        m->print();
        cout << ":";
    }
    cout << endl;
}

bool CCylinder::isConsistent() {

    int index = 0;
    for(auto &m : _monomers) {
        
        if(!m->isConsistent())
            cout << "CMonomer inconsistency is at index "
                 << index << "." << endl;
        
        index++;
    }
    return true;
}

short CCylinder::getType() {
    
    return _pCylinder->getType();
}

int CCylinder::getId(){ return _pCylinder->getId();}

} // namespace medyan

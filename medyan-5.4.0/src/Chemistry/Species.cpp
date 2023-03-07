
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

#include "Species.h"

#include "Reaction.h"
#include "Composite.h"

namespace medyan {
void Species::connect(std::function<void (RSpecies *, int)> callback)
{
    _rspecies->callbacks_.push_back(std::move(callback));
}

Composite* Species::getRoot()
{
    if(hasParent())
       return this->getParent()->getRoot();
    return nullptr;
}

ostream& operator<<(ostream& os, const Species& s){
    os << s.getFullName() << "[" << s.getN() << "]";
    return os;
}

void Species::updateReactantPropensities() {
    
    for(auto r : _rspecies->reactantReactions()){
        r->updatePropensity();
    }
}
//Only used for bulk reactions
void Species::activateReactantReactions() {
    
    for(auto r : _rspecies->reactantReactions()) {

        r->activateReaction();
    }
}


void Species::passivateReactantReactions() {
//    std::cout<<"passivate Species::passivateReactantReactions"<<endl;
    
    for(auto r : _rspecies->reactantReactions())
        r->passivateReaction();
}


unordered_map<string,int> SpeciesNamesDB::_map_string_int;
vector<string> SpeciesNamesDB::_vec_int_string;
unsigned long  SpeciesNamesDB::_num = 0;

} // namespace medyan

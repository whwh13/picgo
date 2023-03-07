
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

#include "ChemSimpleGillespieImpl.h"
#include "Rand.h"

namespace medyan {

void ChemSimpleGillespieImpl::initialize() {
    resetTime();
}

void ChemSimpleGillespieImpl::initializerestart(floatingpoint restarttime){

    if(SysParams::RUNSTATE){
        LOG(ERROR) << "initializerestart Function from ChemSimpleGillespieImpl class can "
                      "only be called "
                      "during restart phase. Exiting.";
        throw std::logic_error("Illegal function call pattern");
    }

    setTime(restarttime);
}


ChemSimpleGillespieImpl::~ChemSimpleGillespieImpl() noexcept {
    _reactions.clear();
}

floatingpoint ChemSimpleGillespieImpl::generateTau(floatingpoint a){
    #ifdef DEBUGCONSTANTSEED
    Rand::chemistrycounter++;
    #endif
    return medyan::rand::safeExpDist(_exp_distr, a, Rand::eng);
}

floatingpoint ChemSimpleGillespieImpl::generateUniform(){
    #ifdef DEBUGCONSTANTSEED
    Rand::chemistrycounter++;
    #endif
    return _uniform_distr(Rand::eng);
}

floatingpoint ChemSimpleGillespieImpl::computeTotalA(){
    floatingpoint rates_sum = 0;
    for (auto &r : _reactions){
        rates_sum+=r->computePropensity();
    }
    return rates_sum;
}

bool ChemSimpleGillespieImpl::makeStep(floatingpoint endTime) {
    
    floatingpoint a_total = computeTotalA();

    floatingpoint tau = generateTau(a_total);
    // Check if a reaction happened before endTime
    if (_t+tau>endTime){ 
        setTime(endTime);
        return true;
    }
    // this means that the network has come to a halt
    if(a_total<1e-15)
        return false;
    _t+=tau;
    syncGlobalTime();
    
    //Gillespie algorithm's second step: finding which reaction happened;
    floatingpoint mu = a_total*generateUniform();
    floatingpoint rates_sum = 0;
    ReactionBase* r_selected = nullptr;
    for (auto &r : _reactions){
        
        rates_sum+=r->computePropensity();
    
        if(rates_sum>mu){
            r_selected = r;
            break;
        }
    }
    if(r_selected==nullptr){
        cout << "ChemSimpleGillespieImpl::makeStep() for loop: rates_sum=" <<
            rates_sum << ", mu="
            << mu << ", a_total=" << a_total << endl;
        throw runtime_error( "ChemSimpleGillespieImpl::makeStep(): No Reaction was selected during the Gillespie step!");
    }
    r_selected->makeStep();

    // Send signal.
    r_selected->emitSignal();
    
    return true;
}

void ChemSimpleGillespieImpl::addReaction(ReactionBase *r) {
    auto vit = find(_reactions.begin(), _reactions.end(), r);
    if(vit==_reactions.end())
        _reactions.push_back(r);
}

void ChemSimpleGillespieImpl::removeReaction(ReactionBase *r) {
    auto vit = find(_reactions.begin(), _reactions.end(), r);
    if(vit!=_reactions.end())
        _reactions.erase(vit);
}

void ChemSimpleGillespieImpl::printReactions() const {
    for (auto &r : _reactions)
        cout << (*r);
}

} // namespace medyan

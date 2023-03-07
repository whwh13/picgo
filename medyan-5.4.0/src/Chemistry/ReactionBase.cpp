
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

#include "ReactionBase.h"

#include "Composite.h"

namespace medyan {
size_t ReactionBase::_Idcounter = 0;

ReactionBase::ReactionBase (FP rate, bool isProtoCompartment, FP volumeFrac, int rateVolumeDepExp)
    : _rnode(nullptr), _parent(nullptr), _rate(rate), 
      _rate_bare(rate), _isProtoCompartment(isProtoCompartment),
      _volumeFrac(volumeFrac), _rateVolumeDepExp(rateVolumeDepExp) {
    
	for(int i = 0; i < RateMulFactorType::RATEMULFACTSIZE; i++)
		_ratemulfactors[i] = 1.0;

    // Scale the rate
	recalcRateVolumeFactor();

    //All reactions are generated passivated.
#if defined TRACK_ZERO_COPY_N || defined TRACK_UPPER_COPY_N
    _passivated=true;
#endif
    _Id = _Idcounter;
    _Idcounter++;
}

Composite* ReactionBase::getRoot() {
    if(hasParent())
        return this->getParent()->getRoot();
    return nullptr;
}

void ReactionBase::registerNewDependent(ReactionBase *r){ _dependents.insert(r);}

void ReactionBase::unregisterDependent(ReactionBase *r){ _dependents.erase(r);}

void ReactionBase::clearSignaling () {
    callbacks_.clear();
}

void ReactionBase::connect(CallbackType callback) {
    callbacks_.push_back(std::move(callback));
}

void ReactionBase::printDependents()  {
    cout << "ReactionBase: ptr=" << this << "\n"
    << (*this) << "the following ReactionBase objects are dependents: ";
    if(_dependents.size()==0)
        cout << "NONE" << endl;
    else
        cout << endl;
    for(auto r : _dependents)
        cout << (*r) << endl;
}

bool afterchemsiminit = false;
void ReactionBase::activateReaction() {
#ifdef TRACK_ZERO_COPY_N
	if(areEqual(getProductOfReactants(), 0.0)) // One of the reactants is still at zero copy n,
		// no need to activate yet...
		return;
#endif
#ifdef TRACK_UPPER_COPY_N
	if(areEqual(getProductOfProducts(), 0.0)) // One of the products is at the maximum allowed
		//copy number, no need to activate yet...
		return;
#endif
	activateReactionUnconditional();
}

} // namespace medyan

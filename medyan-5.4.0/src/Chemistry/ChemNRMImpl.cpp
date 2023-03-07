
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

#include <iostream>
#include <chrono>
#include <cmath>

#ifdef BOOST_MEM_POOL
    #include <boost/pool/pool.hpp>
    #include <boost/pool/pool_alloc.hpp>
#endif

#include "ChemNRMImpl.h"
#include "Chemistry/DissipationTracker.h"
#include "Rand.h"
#include "Controller/CController.h"

#include "CUDAcommon.h"

namespace medyan {

#ifdef BOOST_MEM_POOL
#ifdef BOOST_POOL_MEM_RNODENRM
boost::pool<> allocator_rnodenrm(sizeof(RNodeNRM),BOOL_POOL_NSIZE);
#endif

#ifdef BOOST_POOL_MEM_PQNODE
boost::pool<> allocator_pqnode(sizeof(PQNode),BOOL_POOL_NSIZE);
#endif


#ifdef BOOST_POOL_MEM_PQNODE
void* PQNode::operator new(size_t size) {
    void *ptr = boost::fast_pool_allocator<PQNode>::allocate();
    return ptr;
}

void PQNode::operator delete(void* ptr) noexcept {
    boost::fast_pool_allocator<PQNode>::deallocate((PQNode*)ptr);
}
#endif

#ifdef BOOST_POOL_MEM_RNODENRM
void* RNodeNRM::operator new(size_t size) {
    void *ptr = boost::fast_pool_allocator<RNodeNRM>::allocate();
    return ptr;
}

void RNodeNRM::operator delete(void* ptr) noexcept {
    boost::fast_pool_allocator<RNodeNRM>::deallocate((RNodeNRM*)ptr);
}
#endif
#endif

RNodeNRM::RNodeNRM(ReactionBase *r, ChemNRMImpl &chem_nrm)
    : _chem_nrm (chem_nrm), _react(r) {
    _react->setRnode(this);
     boost_heap *heap = _chem_nrm.getHeap();
    _handle = heap->emplace(this);
    _a = _react->computePropensity();
}

RNodeNRM::~RNodeNRM() noexcept {
    boost_heap *heap = _chem_nrm.getHeap();
    heap->erase(_handle);
    _react->setRnode(nullptr);
}

void RNodeNRM::printSelf() const {
	std::cout.precision(10);
    cout << "RNodeNRM: ptr=" << this <<", tau=" << getTau() <<
        ", a=" << _a <<
        ", points to Reaction:\n";
    cout << (*_react);
}

void RNodeNRM::printDependents() const {
    cout << "RNodeNRM: ptr=" << this
    << ", the following RNodeNRM objects are dependents:\n\n";
    for(auto rit = _react->dependents().begin();
        rit!=_react->dependents().end(); ++rit){
        RNodeNRM *rn_other = (RNodeNRM*)((*rit)->getRnode());
        rn_other->printSelf();
    }
    cout << endl;
}

void RNodeNRM::updateHeap() {
    boost_heap *heap = _chem_nrm.getHeap();
    heap->update(_handle);
}

void RNodeNRM::generateNewRandTau() {
    floatingpoint newTau;
    reComputePropensity();//calculated new _a

#ifdef TRACK_ZERO_COPY_N
    auto t1 = _chem_nrm.generateTau(_a);
    auto t2 = _chem_nrm.getTime();
    newTau = t1 + t2;
#else
    if(_a<1.0e-10) // numeric_limits< floatingpoint >::min()
        newTau = numeric_limits<floatingpoint>::infinity();
    else
        newTau = _chem_nrm.generateTau(_a) + _chem_nrm.getTime();
#endif
//    cout<<"Propensity of rxn "<<_a<<" tau "<<newTau<<endl;
    setTau(newTau);

}

void RNodeNRM::activateReaction() {
    generateNewRandTau();
    updateHeap();
}

void RNodeNRM::passivateReaction() {
    _a=0;
    floatingpoint tau = numeric_limits<floatingpoint>::infinity();
    setTau(tau);
    updateHeap();
}

void ChemNRMImpl::initialize() {
    resetTime();
    for (auto &x : _map_rnodes){
        auto rn = x.second.get();
        rn->getReaction()->activateReaction();
    }
}

void ChemNRMImpl::initializerestart(floatingpoint restarttime){

    if(SysParams::RUNSTATE){
        log::error("initializerestart Function from ChemSimpleGillespieImpl class can only be called during restart phase. Exiting.");
        throw std::logic_error("Illegal function call pattern");
    }

    setTime(restarttime);
}


ChemNRMImpl::~ChemNRMImpl() {
    _map_rnodes.clear();
}

floatingpoint ChemNRMImpl::generateTau(floatingpoint a){

	#ifdef DEBUGCONSTANTSEED
	Rand::chemistrycounter++;
	#endif
    return medyan::rand::safeExpDist(_exp_distr, a, Rand::eng);
}

bool ChemNRMImpl::makeStep(floatingpoint endTime) {
    chrono::high_resolution_clock::time_point mins, mine, minsT, mineT, minses, mintes;
    minsT = chrono::high_resolution_clock::now();
    //try to get a reaction
    if(_heap.empty()) {
        if (endTime==std::numeric_limits<floatingpoint>::infinity()){
            cout << "There are no reactions to fire, returning..." << endl;
            return false;
        } else {
            setTime(endTime);
            return true;
        }
    }
    RNodeNRM *rn = _heap.top()._rn;
    const floatingpoint tau_top = rn->getTau();
    // Check if a reaction happened before endTime
    if (tau_top>endTime){ 
        setTime(endTime);
        return true;
    }
    if(tau_top==numeric_limits<floatingpoint>::infinity()){

        log::info("The heap has been exhausted - no more reactions to fire, returning...");
        return false;
    }
    ///DEBUG
    //assert heap ordering
    if(tau_top < _t) {
        log::error("ERROR: The heap is not correctly sorted, returning...");
        log::info("Tau top = {}", tau_top);
        log::info("Tau current = {}", _t);
        log::info("Reaction type = {}", text(rn->getReaction()->getReactionType()));
        rn->printSelf();
        return false;
    }

    floatingpoint t_prev = _t;

    setTime(tau_top);
    // if dissipation tracking is enabled and the reaction is supported, then compute the change in Gibbs free energy and store it
    if(SysParams::Chemistry().dissTracking){
        ReactionBase* react = rn->getReaction();
        string HRCDID = react->getHRCDID();
        string testString = "DNT";
        if((HRCDID != testString) || (react->getReactionType() == ReactionType::DIFFUSION)){
            dt->updateDelGChem(react);
        }
    }
    #ifdef CROSSCHECK_CYLINDER
    auto _react = rn->getReaction();
    if(_react->getReactionType()!= ReactionType::DIFFUSION){
        auto _a = rn->getPropensity();
        Compartment* c = static_cast<Compartment*>(_react->getParent());
        auto coord = c->coordinates();
        CController::_crosscheckdumpFilechem << "RNodeNRM: ptr=" << this <<", tau=" <<
                                             rn->getTau() <<
                                             ", a=" << _a <<" in Compartment "<<coord[0]<<" "<<coord[1]<<" "<<coord[2]<<
                                             ", points to Reaction Type "<< _react->getReactionType()<<endl;
//    CController::_crosscheckdumpFilechem << (*_react);

    }
    else{
        CController::_crosscheckdumpFilechem << "DIFFUSION "<<endl;
    }

    #endif
    rn->makeStep();

    #ifdef CROSSCHECK_CYLINDER
    CController::_crosscheckdumpFilechem <<"Update dependencies"<<endl;
    #endif

	#ifdef DEBUGCONSTANTSEED
    cout<<"tau "<<_t<<endl;
    #endif
    mins = chrono::high_resolution_clock::now();
#if defined TRACK_ZERO_COPY_N || defined TRACK_UPPER_COPY_N
    if(!rn->isPassivated()){
#endif
        //std::cout<<"Update R and Tau for fired reaction"<<endl;
        rn->generateNewRandTau();
        rn->updateHeap();

#if defined TRACK_ZERO_COPY_N || defined TRACK_UPPER_COPY_N
    }
#endif
    // Updating dependencies
    ReactionBase *r = rn->getReaction();
    if(r->updateDependencies()) {

        for(auto& prdep : r->dependents()) {

            RNodeNRM *rn_other = (RNodeNRM*)(prdep->getRnode());
            floatingpoint a_old = rn_other->getPropensity();

            //recompute propensity
            rn_other->reComputePropensity();

            floatingpoint tau_new;
            floatingpoint tau_old = rn_other->getTau();

            floatingpoint a_new = rn_other->getPropensity();

            // Recompute tau.
            // Unlikely exceptional cases for better branch prediction.
            if(a_new == 0.0 || a_new == numeric_limits<floatingpoint>::infinity() || a_old == 0.0 || a_old == numeric_limits<floatingpoint>::infinity()) [[unlikely]] {
                if(a_new == 0.0) {
                    tau_new = numeric_limits<floatingpoint>::infinity();
                }
                else if(a_new == numeric_limits<floatingpoint>::infinity()) {
                    // propensity can occasionally be infinite, indicating instant reaction.
                    tau_new = _t;
                }
                else {
                    rn_other->generateNewRandTau();
                    tau_new = rn_other->getTau();
                }
            }
            else {
                tau_new = (a_old/a_new)*(tau_old-_t) + _t;
            }

            // Debug.
            if(std::isnan(tau_new)) {
                log::error("tau_new is nan");
                log::info("a_old={}, a_new={}, tau_old={}, _t={}", a_old, a_new, tau_old, _t);
                throw std::runtime_error("tau_new is nan");
            }
            if(tau_new < _t) {

                log::warn("WARNING: Generated tau may be incorrect.");

                log::info("Tau new = {}", tau_new);
                log::info("Tau old = {}", tau_old);
                log::info("Current global t = {}", _t);
                log::info("Previous global t = {}", t_prev);
                log::info("a_old = {}", a_old);
                log::info("a_new = {}", a_new);

                log::info("Reaction type = {}", text(rn->getReaction()->getReactionType()));


                rn->printSelf();
                rn_other->printSelf();
            }

            rn_other->setTau(tau_new);
            rn_other->updateHeap();
        }
    }
    mine = chrono::high_resolution_clock::now();

#ifdef OPTIMOUT
    chrono::duration<floatingpoint> elapsed_time(mine - mins);
    auto rType = r->getReactionType();
    if(rType == 1){
        auto reactant = r->getReactantCopyNumbers();
        auto product = r->getProductCopyNumbers();
        if(reactant[0] == 0){
            CUDAcommon::cdetails.diffusion_passivate_count++;

        }else if (product[0] == 1){
            CUDAcommon::cdetails.diffusion_activate_count++;
        }
    }
    CUDAcommon::cdetails.reactioncount[rType]++;
    CUDAcommon::cdetails.dependencytime[rType]+= elapsed_time.count();
    CUDAcommon::cdetails.dependentrxncount[rType] += r->dependents().size();
#endif

    #ifdef CROSSCHECK_CYLINDER
    CController::_crosscheckdumpFilechem <<"emitSignal"<<endl;
    #endif
    minses = chrono::high_resolution_clock::now();
    // Send signal.
    r->emitSignal();

#ifdef OPTIMOUT
    mintes = chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_emitsignal(mintes - minses);
    CUDAcommon::cdetails.emitsignal[rType]+= elapsed_emitsignal.count();
#endif

    #ifdef CROSSCHECK_CYLINDER
    CController::_crosscheckdumpFilechem <<"----"<<endl;
    #endif
    mineT = chrono::high_resolution_clock::now();
#ifdef OPTIMOUT
    chrono::duration<floatingpoint> elapsed_timetotal(mineT - minsT);
    CUDAcommon::cdetails.totaltime[rType]+= elapsed_timetotal.count();
#endif
    return true;
}

void ChemNRMImpl::addReaction(ReactionBase *r) {
    _map_rnodes.emplace(r,make_unique<RNodeNRM>(r,*this));
}

void ChemNRMImpl::removeReaction(ReactionBase *r) {
    _map_rnodes.erase(r);
}

void ChemNRMImpl::printReactions() const {
    for (auto &x : _map_rnodes){
        auto rn = x.second.get();
        rn->printSelf();
    }
}

bool ChemNRMImpl::crosschecktau() const {
    bool status = true;
    for (auto &x : _map_rnodes){
        auto rn = x.second.get();
        if(rn->getTau() < tau()) {
            rn->printSelf();
            status = false;
            log::error("Tau in reaction is smaller than current time.");
            throw std::runtime_error("Tau in reaction is smaller than current time.");
        }
    }
    return status;
}

} // namespace medyan

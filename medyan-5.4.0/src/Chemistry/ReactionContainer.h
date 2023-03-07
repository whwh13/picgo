
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

#ifndef MEDYAN_ReactionContainer_h
#define MEDYAN_ReactionContainer_h

#include "common.h"

#include "Reaction.h"

namespace medyan {
/// An abstract interface for a container of pointers to reaction objects
class ReactionPtrContainerIFace {
public:
    ///Clear the reaction container
    virtual void clear() = 0;
    
    ///Add a unique reaction pointer to this container
    virtual ReactionBase* addReactionUnique(unique_ptr<ReactionBase> &&Reaction) = 0;
    
    ///Remove a reaction from this container
    virtual void removeReaction(ReactionBase* Reaction) = 0;
    
    ///Find a reaction in this container
    virtual ReactionBase* findReaction (size_t index) = 0;
    
    ///Print all reactions from this container
    virtual void printReaction() {}
};

/// A concrete class implementing the ReactionPtrContainerIFace,
/// using vector<unique_ptr<ReactionBase>> as the container implementation
class ReactionPtrContainerVector : public  ReactionPtrContainerIFace {
protected:
    vector<unique_ptr<ReactionBase>> _reactions;  ///< Reaction ptr container
public:
    
    /// Default constructor
    ReactionPtrContainerVector() : _reactions() {}
    
    /// Copying not allowed
    ReactionPtrContainerVector(const ReactionPtrContainerVector &) = delete;
    ReactionPtrContainerVector(ReactionPtrContainerVector &&)      = default;
    
    /// Assignment not allowed
    ReactionPtrContainerVector& operator=(ReactionPtrContainerVector &) = delete;

    friend void swap(ReactionPtrContainerVector& first,
                     ReactionPtrContainerVector& second) // nothrow
    {
        // enable ADL (not necessary in our case, but good practice)
        using std::swap;
        swap(first._reactions, second._reactions);
    }

    /// Clear the container
    virtual void clear() {_reactions.clear();}
    /// Add a unique reaction ptr to this container
    virtual ReactionBase* addReactionUnique (unique_ptr<ReactionBase> &&Reaction) {
        _reactions.push_back(move(Reaction));
        return _reactions.back().get();
    }
    
    /// Add a general reaction to this container
    template<unsigned short M, unsigned short N, typename ...Args>
    ReactionBase* addReaction( Args&& ...args )
    {
        _reactions.push_back(unique_ptr<ReactionBase>(
            new Reaction<M,N>( forward<Args>(args)...) ));
        return _reactions.back().get();
    }
    
    /// Add a general reaction class to this container
    template<template <unsigned short M, unsigned short N> class RXN,
                       unsigned short M, unsigned short N>
    ReactionBase* add(initializer_list<Species*> species, float rate)
    {
        _reactions.push_back(unique_ptr<ReactionBase>( new RXN<M,N>(species,rate)));
        return _reactions.back().get();
    }
    
    /// Remove a reaction from this container
    virtual void removeReaction (ReactionBase* R) {
        auto child_iter = find_if(_reactions.begin(),_reactions.end(),
                [R](const unique_ptr<ReactionBase> &element) {
                    return element.get()==R ? true : false; });
        if(child_iter!=_reactions.end()) {
            _reactions.erase(child_iter);
        }
    }
   //aravind June 24, 2016.
    virtual void updatePropensityComprtment()
    {        for(auto &it: _reactions)
        it->updatePropensity();    }
    
    /// Remove all reactions that contain a certain species from this container
    virtual void removeReactions (Species* s)
    {
        for(auto &r : _reactions)
        {
            if(r->containsSpecies(s)){
                removeReaction(r.get());
                return;
            }
        }
    }
    
    /// Find a reaction by index in this container
    /// @note no check on the index
    virtual ReactionBase* findReaction (size_t index) {
        return _reactions[index].get();
    }
    
    /// Get all reactions in vector form
    vector<unique_ptr<ReactionBase>>& reactions() {return _reactions;}
    const vector<unique_ptr<ReactionBase>>& reactions() const {return _reactions;}
    
    /// Print all reactions in this container
    virtual void printReactions()const {
        for(auto &r : _reactions)
            cout << (*r.get());
    }
    
    /// Find a similar reaction in this container (satistifes equality operator)
    /// @note returns the first similar reaction found
    virtual ReactionBase* findSimilarReaction (const ReactionBase &r) {
        auto it = find_if(_reactions.begin(),_reactions.end(),
                               [&r](const unique_ptr<ReactionBase> &element)
                               {return r==(*element);});
        if(it==_reactions.end()) throw out_of_range(
        "Reaction::findSimilarReaction(): The analogous Reaction was not found");
        return it->get();
    }
};

} // namespace medyan

#endif

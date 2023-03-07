
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

#ifndef MEDYAN_Reaction_h
#define MEDYAN_Reaction_h

#include "common.h"

#include "ReactionBase.h"

namespace medyan {
/// Represents a concrete chemical reaction, such as A + B -> C, where M is the number
/// of reactants and N is the number of products.

/*! Reaction<M,N> encodes a chemical reaction between M reactants and N products. It 
 *  follows the ReactionBase interface, where many methods are defined. Most of the 
 *  methods defined in Reaction<M,N> are specific implementations of virtual functions
 *  declared in ReactionBase. Hence, to a very large degree, Reaction<M,N> is an 
 *  implementation class, while ReactionBase provides the public interface for reaction 
 *  objects. The documentation of the latter should be mainly consulted when
 *  working with reactions.
 */
template <unsigned short M, unsigned short N>
class Reaction : public ReactionBase {
protected:
        array<RSpecies*, M+N> _rspecies; ///< An array of RSpecies objects
                                         ///< (reactants followed by products)
public:
        /// The main constructor:
        /// @param species - are reactants and products put together into a single list
        /// (starting from reactants)
        /// @param rate - the rate constant for this ReactionBase
        Reaction(initializer_list<Species*> species,
                FP rate = 0.0, bool isProtoCompartment = false, FP volumeFrac = 1.0,
                int rateVolumeDepExp = 0)
                : ReactionBase(rate, isProtoCompartment, volumeFrac, rateVolumeDepExp) {
            initializeSpecies(species);
        }
        
        /// The main constructor:
        /// @param it_begin - an iterator to the beginning of an RSpecies* container
        /// @param it_end - an iterator to the end of an RSpecies* container
        /// @param rate - the rate constant for this ReactionBase
        template <typename InputContainer>
        Reaction(const InputContainer &species,
                 FP rate = 0.0, bool isProtoCompartment = false, FP volumeFrac = 1.0,
                 int rateVolumeDepExp = 0)
                : ReactionBase(rate, isProtoCompartment, volumeFrac, rateVolumeDepExp) {
            initializeSpecies(species);
        }
        
        /// no copying (including all derived classes)
        Reaction (const Reaction &rb) = delete;
        
        /// no assignment (including all derived classes)
        Reaction& operator=(Reaction &rb) = delete;

#ifdef BOOST_MEM_POOL
        /// Advanced memory management
        void* operator new(size_t size);
        
        void operator delete(void* ptr) noexcept;
#endif
        /// Destructor
        /// Tell Rspecies to remove this Reaction from its internal lists of reactions
        /// @note noexcept is important here. Otherwise, gcc flags the constructor as
        /// potentially throwing, which in turn disables move operations by the STL
        /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will
        /// presumbaly be fixed in the future.
        virtual ~Reaction() noexcept
        {
            for(auto i=0U; i<M; ++i) _rspecies[i]->removeAsReactant(this);
            for(auto i=M; i<(M+N); ++i) _rspecies[i]->removeAsProduct(this);
        }
        
        /// Returns a pointer to the first element of array<RSpecies*, M+N>
        /// The type of the pointer is RSpecies**. In conjunction with getM() and
        /// getN(), this pointer can be used to iterate over RSpecies associated with
        /// this reaction.
        inline virtual RSpecies** rspecies() override {return &_rspecies[0];}

    virtual void setDependentReactions() override {
        _dependents.clear();
        for(int i = 0; i < M + N; i++) {
            auto s = _rspecies[i];
            for(auto it = s->beginReactantReactions();
                it != s->endReactantReactions(); it++) {
                ReactionBase* r = (*it);
                if(r!=this && !r->isPassivated())
                    _dependents.insert(r);
            }
        }
    }

    virtual void updatePropensityImpl() override;

protected:
        /// An implementation method used by the constructor.
        template <typename InputContainer>
        void initializeSpecies(const InputContainer &species){
            assert(species.size()==(M+N)
                   && "Reaction<M,N> Ctor: The species number does not match the template M+N");
            transform(species.begin(),species.end(),_rspecies.begin(),
                      [](Species *s){return &s->getRSpecies();});

            if(!_isProtoCompartment) {
#ifdef TRACK_DEPENDENTS
                //add dependents
                setDependentReactions();
#endif
                for(auto i=0U; i<M; ++i) _rspecies[i]->addAsReactant(this);
                for(auto i=M; i<(M+N); ++i) _rspecies[i]->addAsProduct(this);
            }
        }
        
        /// Implementation of getM()
        inline virtual unsigned short getMImpl() const override {return M;}

        /// Implementation of getN()
        inline virtual unsigned short getNImpl() const override {return N;}

        /// Implementation of size()
        inline virtual unsigned short sizeImpl() const override {return M+N;}

        /// Implementation of the operator==(...)
        virtual bool is_equal(const ReactionBase& other) const override
        {
            const Reaction<M,N> *a = this;
            const Reaction<M,N> *b = static_cast<const Reaction<M,N>*>(&other);
            auto it_pair = mismatch(a->_rspecies.begin(),a->_rspecies.end(),
                                    b->_rspecies.begin(),
            [](RSpecies* A, RSpecies* B){return A->getSpecies()==B->getSpecies();});
            if(it_pair.first==a->_rspecies.end())
                return true;
            return false;
        }
        
        /// Implementation of getProductOfReactants()
        inline virtual floatingpoint getProductOfReactantsImpl() const override
        {
            floatingpoint prod = 1;
            for(auto i=0U; i<M; ++i)
                prod*=_rspecies[i]->getN();
            return prod;
            
        }

        /// Implementation of computePropensity()
        inline virtual floatingpoint computePropensityImpl() const override
        {
            if(isPassivated()) return 0.0;
#ifdef TRACK_UPPER_COPY_N
            if(areEqual(this->Reaction<M,N>::getProductOfProductsImpl(),0.0)){
                return 0.0;
            }
#endif
            return _rate*Reaction<M,N>::getProductOfReactantsImpl();
        }
        
        /// Implementation of getProductOfProducts()
        inline virtual floatingpoint getProductOfProductsImpl() const override
        {
#ifdef TRACK_UPPER_COPY_N
            floatingpoint prod = 1;
            for(auto i=M; i<(M+N); ++i){
                prod*=_rspecies[i]->getN()-_rspecies[i]->getUpperLimitForN();
            }
            return prod;
#else
            return 1.0;
#endif
        }
    
        /// Implementation of containsSpecies()
        inline virtual bool containsSpeciesImpl(Species *s) const override
        {
            auto it = find_if(_rspecies.begin(), _rspecies.end(),
                [s](const RSpecies *rs){return (&rs->getSpecies())==s;});
            return it!=_rspecies.end();
            
        }
    
        /// Implementation of makeStep()
        inline virtual void makeStepImpl() override
        {
            for(auto i=0U; i<M; ++i) _rspecies[i]->down();
            for(auto i=M; i<(M+N); ++i) _rspecies[i]->up();
        }

        /// Implementation of activateReactionUnconditional()
        virtual void activateReactionUnconditionalImpl() override;

        /// Implementation of passivateReaction()
        virtual void passivateReactionImpl() override;
        
        /// Print information about this reaction to ostream
        virtual void printToStream(ostream& os) const override
        {
            unsigned char i=0;
            auto sit = _rspecies.cbegin();
            auto send = sit+M;
            for (; sit!=send; ++sit)
            {
                os << (*sit)->getFullName() << "{" << (*sit)->getN() << "}";
                if(i<M-1)
                os << " + ";
                ++i;
            }
            os << " ---> ";
            i=0;
            for (auto sit2 = send; sit2!=_rspecies.cend(); ++sit2)
            {
                os << (*sit2)->getFullName() << "{" << (*sit2)->getN() << "}";
                if(i<(N-1))
                    os << " + ";
                ++i;
            }
            os << ", " << "curr_rate = " << getRate() << ", a="
//                    << computePropensity() << "\n";
               << computePropensity() << ", ReactionBase ptr=" << this <<"\n";
        }
        
        /// Implementation of  clone()
        virtual Reaction<M,N>* cloneImpl(
            const SpeciesPtrContainerVector &spcv) override;
        
        virtual bool updateDependencies() override {return true;}
    };


/// A diffusive reaction in the system.
/*!
 * A DiffusionReaction is very similar to the Reaction class, the only difference
 * in being the marking of a few fields. The dependencies marker will be set to true
 * when a new average for the held RSpeciesAvg is calculated (if at all). 
 * The averaging marker specifies if this diffusion reaction is using RSpeciesAvg 
 * as its product and reactant. 
 */

class DiffusionReaction : public Reaction<1,1> {
    
private:
    bool _dependencies = true; ///< A marker to represent whether the dependents of
                               ///< this reaction should be updated. This will change
                               ///< based upon the diffusing species.
    
    bool _averaging = false;  ///< Using RSpeciesAvg as reactant and product
    
public:
    /// The main constructor
    DiffusionReaction(initializer_list<Species*> species,
                      FP rate = 0.0, bool isProtoCompartment = false, FP volumeFrac = 1.0)
        : Reaction(species, rate, isProtoCompartment, volumeFrac, -1) {
    
        //set averaging
        if(dynamic_cast<RSpeciesAvg*>(_rspecies[0]))
            _averaging = true;

        //set type
        _reactionType = ReactionType::DIFFUSION;
    }
    
    //Destructor does nothing new
    virtual ~DiffusionReaction() {}
    
#ifdef BOOST_MEM_POOL
    /// Advanced memory management
    void* operator new(size_t size);
    
    void operator delete(void* ptr) noexcept;
#endif
    
    /// Implementation of makeStep()
    inline virtual void makeStepImpl() override
    {
        _rspecies[0]->down();
        _rspecies[1]->up();
    
        //if averaging, update dependency marker
        if(_averaging) {
            
            bool newAvgR = ((RSpeciesAvg*)_rspecies[0])->newAverage();
            bool newAvgP = ((RSpeciesAvg*)_rspecies[1])->newAverage();
        
            _dependencies = newAvgR || newAvgP;
        }
    }

    ///This implementation returns the kept boolean value
    virtual bool updateDependencies() override {return _dependencies;}

    virtual void updatePropensityImpl() override;
};

} // namespace medyan

#endif


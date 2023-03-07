#ifndef MEDYAN_Chemistry_ReactionDy_hpp
#define MEDYAN_Chemistry_ReactionDy_hpp

#include <algorithm>
#include <stdexcept>

#ifdef BOOST_MEM_POOL
    #include <boost/pool/pool.hpp>
    #include <boost/pool/pool_alloc.hpp>
#endif

#include "Chemistry/ChemRNode.h"
#include "Chemistry/ReactionBase.h"
#include "Util/Io/Log.hpp"

namespace medyan {

class ReactionDy : public ReactionBase {
public:
    // Repetitive rspecies
    struct RepRSpecies {
        RSpecies*      prs    = nullptr;
        unsigned short stoich = 1;
    };

private:
    std::vector< RepRSpecies > repRSpecies_;
    unsigned short numReactants_ = 0;
    unsigned short numProducts_  = 0;
    // sum of stoich numbers
    unsigned short stoichReactants_ = 0;
    unsigned short stoichProducts_  = 0;

    void registerInRSpecies_() {
        for(unsigned short i = 0; i < numReactants_; ++i) {
            repRSpecies_[i].prs->addAsReactant(this);
        }
        for(unsigned short i = 0; i < numProducts_; ++i) {
            repRSpecies_[i + numReactants_].prs->addAsProduct(this);
        }
    }
    void unregisterInRSpecies_() {
        for(unsigned short i = 0; i < numReactants_; ++i) {
            repRSpecies_[i].prs->removeAsReactant(this);
        }
        for(unsigned short i = 0; i < numProducts_; ++i) {
            repRSpecies_[i + numReactants_].prs->removeAsProduct(this);
        }
    }

public:

    template< typename SpeciesPtrContainer >
    ReactionDy(
        const SpeciesPtrContainer& reactantSpecies,
        const SpeciesPtrContainer& productSpecies,
        ReactionType               reactionType,
        FP                         rate,
        FP                         volumeFrac = 1.0,
        int                        rateVolumeDepExp = 0
    ) : ReactionBase(rate, false, volumeFrac, rateVolumeDepExp) {
        _reactionType = reactionType;
        initializeSpecies(reactantSpecies, productSpecies);
    }
        
    ReactionDy(const ReactionDy &) = delete;
        
#ifdef BOOST_MEM_POOL
    /// Advanced memory management
    void* operator new(size_t size) {
        void *ptr = boost::fast_pool_allocator< ReactionDy >::allocate();
        return ptr;
    }

    void operator delete(void* ptr) noexcept {
        boost::fast_pool_allocator< ReactionDy >::deallocate((ReactionDy*)ptr);
    }
#endif

    virtual ~ReactionDy() noexcept override {
        unregisterInRSpecies_();
    }

    ReactionDy& operator=(const ReactionDy &) = delete;


    auto& getRepRSpecies() { return repRSpecies_; }
    auto& getRepRSpecies() const { return repRSpecies_; }

    // Add affected reactions to the dependency list.
    virtual void setDependentReactions() override {
        _dependents.clear();
        for(int i = 0; i < repRSpecies_.size(); i++) {
            for(auto r : repRSpecies_[i].prs->reactantReactions()) {
                if(r != this && !r->isPassivated()) {
                    _dependents.insert(r);
                }
            }
        }
    }

    virtual void updatePropensityImpl() override {
        if(_rnode && !_passivated) _rnode->activateReaction();
    }

    template <typename SpeciesPtrContainer>
    void initializeSpecies(
        const SpeciesPtrContainer& reactantSpecies,
        const SpeciesPtrContainer& productSpecies
    ) {
        using namespace std;

        // Unregister this in all RSpecies
        unregisterInRSpecies_();

        stoichReactants_ = reactantSpecies.size();
        stoichProducts_  = productSpecies.size();

        repRSpecies_.clear();

        // Set reactants
        for(Species* ps : reactantSpecies) {
            const auto prs = &ps->getRSpecies();
            if(
                auto it = find_if(
                    repRSpecies_.begin(), repRSpecies_.end(),
                    [&](const RepRSpecies& rrs) { return rrs.prs == prs; }
                );
                it == repRSpecies_.end()
            ) {
                repRSpecies_.push_back({ prs, 1 });
            } else {
                ++it->stoich;
            }
        }
        numReactants_ = repRSpecies_.size();

        // Set products
        for(Species* ps : productSpecies) {
            const auto prs = &ps->getRSpecies();
            if(
                auto it = find_if(
                    repRSpecies_.begin() + numReactants_, repRSpecies_.end(),
                    [&](const RepRSpecies& rrs) { return rrs.prs == prs; }
                );
                it == repRSpecies_.end()
            ) {
                repRSpecies_.push_back({ prs, 1 });
            } else {
                ++it->stoich;
            }
        }
        numProducts_ = repRSpecies_.size() - numReactants_;
        
        // Add dependents.
        setDependentReactions();

        // Register this in all RSpecies
        registerInRSpecies_();
    }
        
    virtual unsigned short getMImpl() const override { return stoichReactants_; }
    virtual unsigned short getNImpl() const override { return stoichProducts_; }
    virtual unsigned short sizeImpl() const override { return stoichReactants_ + stoichProducts_; }

    /// Implementation of getProductOfReactants()
    virtual floatingpoint getProductOfReactantsImpl() const override {
        floatingpoint prod = 1;
        for(unsigned short i = 0; i < numReactants_; ++i) {
            const auto& rrs = repRSpecies_[i];
            auto copyNumber = static_cast<int>(rrs.prs->getN());
            for(unsigned j = 0; j < rrs.stoich; ++j) {
                prod *= copyNumber;

                --copyNumber;
                if(copyNumber <= 0) break;
            }
        }
        return prod;
    }

    /// Implementation of computePropensity()
    virtual floatingpoint computePropensityImpl() const override
    {
        if(isPassivated()) return 0.0;
#ifdef TRACK_UPPER_COPY_N
        if(areEqual(getProductOfProductsImpl(), 0.0)){
            return 0.0;
        }
#endif
        return _rate * getProductOfReactantsImpl();
    }
        
    /// Implementation of getProductOfProducts()
    virtual floatingpoint getProductOfProductsImpl() const override
    {
#ifdef TRACK_UPPER_COPY_N
        floatingpoint prod = 1;
        for(unsigned short i = 0; i < numProducts_; ++i) {
            const auto& rrs = repRSpecies_[i + numReactants_];
            auto       copyNumber    = static_cast<int>(rrs.prs->getN());
            const auto copyNumberMax = static_cast<int>(rrs.prs->getUpperLimitForN());
            for(unsigned j = 0; j < rrs.stoich; ++j) {
                prod *= copyNumber - copyNumberMax;
                ++copyNumber;
            }
        }
        return prod;
#else
        return 1.0;
#endif
    }


    /// Implementation of makeStep()
    virtual void makeStepImpl() override
    {
        for(unsigned short i = 0; i < numReactants_; ++i) {
            const auto& rrs = repRSpecies_[i];
            for(unsigned j = 0; j < rrs.stoich; ++j) {
                rrs.prs->down();
            }
        }
        for(unsigned short i = 0; i < numProducts_; ++i) {
            const auto& rrs = repRSpecies_[i + numReactants_];
            for(unsigned j = 0; j < rrs.stoich; ++j) {
                rrs.prs->up();
            }
        }
    }

    /// Implementation of activateReactionUnconditional()
    virtual void activateReactionUnconditionalImpl() override {
#ifdef TRACK_DEPENDENTS
        for(unsigned short i = 0; i < numReactants_; ++i)
        {
            auto& rs = *repRSpecies_[i].prs;
            for(auto r : rs.reactantReactions()) {
                if(this != r) r->registerNewDependent(this);
            }
            for(auto r : rs.productReactions()) {
                if(this != r) r->registerNewDependent(this);
            }
        }
#endif
#if defined TRACK_ZERO_COPY_N || defined TRACK_UPPER_COPY_N
        _passivated = false;
#endif

        if(_rnode) _rnode->activateReaction();
    }

    /// Implementation of passivateReaction()
    virtual void passivateReactionImpl() override {
        if(isPassivated()) return;
#ifdef TRACK_DEPENDENTS
        for(unsigned short i = 0; i < numReactants_; ++i)
        {
            auto& rs = *repRSpecies_[i].prs;
            for(auto r : rs.reactantReactions()) {
                r->unregisterDependent(this);
            }
            for(auto r : rs.productReactions()) {
                r->unregisterDependent(this);
            }
        }
#endif
#if defined TRACK_ZERO_COPY_N || defined TRACK_UPPER_COPY_N
        _passivated = true;
#endif
        if(_rnode) _rnode->passivateReaction();
    }


    /// Print information about this reaction to ostream
    virtual void printToStream(ostream& os) const override
    {
        unsigned char i=0;
        auto sit = repRSpecies_.cbegin();
        auto send = sit + numReactants_;
        for (; sit!=send; ++sit)
        {
            os << (*sit).prs->getFullName() << "{" << (*sit).prs->getN() << "}";
            if(i < numReactants_ - 1)
                os << " + ";
            ++i;
        }
        os << " ---> ";
        i=0;
        for (auto sit2 = send; sit2 != repRSpecies_.cend(); ++sit2)
        {
            os << (*sit2).prs->getFullName() << "{" << (*sit2).prs->getN() << "}";
            if(i< numProducts_ - 1)
                os << " + ";
            ++i;
        }
        os << ", " << "curr_rate = " << getRate() << ", a="
            << computePropensity() << ", ReactionBase ptr=" << this <<"\n";
    }

    virtual bool updateDependencies() override { return true; }

    // Only for the purpose of overriding
    // If possible, these functions should be removed.
    //---------------------------------

    /// Returns a pointer to the first element of array<RSpecies*, M+N>
    /// The type of the pointer is RSpecies**. In conjunction with getM() and
    /// getN(), this pointer can be used to iterate over RSpecies associated with
    /// this reaction.
    virtual RSpecies** rspecies() override {
        log::error("RSpecies is not stored in a contiguous array in ReactionDy.");
        throw std::runtime_error("Invalid rspecies() function in ReactionDy");
    }

    virtual bool is_equal(const ReactionBase& other) const override {
        log::error("Virtual equality comparison of ReactionDy should not be used");
        throw std::runtime_error("Invalid comparsion of ReactionDy");
    }

    /// Implementation of containsSpecies()
    virtual bool containsSpeciesImpl(Species *s) const override
    {
        auto it = find_if(repRSpecies_.begin(), repRSpecies_.end(),
            [s](const RepRSpecies& rrs){ return (&rrs.prs->getSpecies()) == s; });
        return it != repRSpecies_.end();
    }

    /// Implementation of  clone()
    virtual ReactionDy* cloneImpl(const SpeciesPtrContainerVector &spcv) override {
        log::error("Virtual clone function in ReactionDy is not allowed.");
        throw std::runtime_error("Invalid clone in ReactionDy");
    }

};

} // namespace medyan

#endif

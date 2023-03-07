
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

#include <catch2/catch.hpp>

#include "Species.h"
#include "Reaction.h"
#include "Compartment.h"

namespace medyan {

void rspecies_callback (RSpecies *r, int delta){
    r->getSpecies().setN(33);
}

void reaction_callback (ReactionBase *r){
    r->setBareRate(0.05);
}

struct ReactionCallback {
    void operator() (ReactionBase *r){
        ++_count;
        r->setBareRate(1.05);
    }
    int _count;
};

TEST_CASE("RSpecies test", "[RSpecies]") {
    // Light testing RSpecies by itself, without Reactions
    Species A("A",  10, max_ulim, SpeciesType::BULK, RSpeciesType::REG);
    RSpecies& RA(A.getRSpecies());
    CHECK(10 == RA.getN());
    CHECK(&A == &RA.getSpecies());

    CHECK(0 == RA.reactantReactions().size());
    CHECK(0 == RA.productReactions().size());

    Species B("B",  10, max_ulim, SpeciesType::BULK, RSpeciesType::REG);
    RSpecies& RB(B.getRSpecies());
    
    // Using a reaction A->B to see if copy number n is correctly driven
    // (implicitly tests up() and down() methods
    Reaction<1,1> rxn1 = {{&A,&B}, 10.0 };
    rxn1.makeStep();
    rxn1.makeStep();
    CHECK(8 == RA.getN());
    CHECK(12 == RB.getN());
    
    // Introducing the back reaction, B->A, to bring the RSpecies back to
    // their original state
    Reaction<1,1> rxn2 = {{&B,&A}, 10.0 };
    rxn2.makeStep();
    rxn2.makeStep();
    CHECK(10 == RA.getN());
    CHECK(10 == RB.getN());    
}

//Testing const rspecies
TEST_CASE("RSpecies counting", "[RSpecies]") {
    
    Species A("A",  10, max_ulim, SpeciesType::BULK, RSpeciesType::CONST);
    RSpecies& RA(A.getRSpecies());
    
    Species B("B",  10, max_ulim, SpeciesType::BULK, RSpeciesType::CONST);
    RSpecies& RB(B.getRSpecies());
    
    // Using a reaction A->B to see if copy number n is correctly driven
    // (implicitly tests up() and down() methods
    Reaction<1,1> rxn1 = {{&A,&B}, 10.0 };
    rxn1.makeStep();
    rxn1.makeStep();
    CHECK(10 == RA.getN());
    CHECK(10 == RB.getN());
    
    Reaction<1,1> rxn2 = {{&B,&A}, 10.0 };
    rxn2.makeStep();
    rxn2.makeStep();
    CHECK(10 == RA.getN());
    CHECK(10 == RB.getN());
}

//Testing avg rspecies
TEST_CASE("Average RSpecies", "[RSpecies]") {
    
    global_time = 0.0;
    
    Species A("A", 100, max_ulim, SpeciesType::BULK, RSpeciesType::AVG);
    RSpecies& RA(A.getRSpecies());
    RSpeciesAvg* RAA = (RSpeciesAvg*)&RA;
    
    RAA->setNumEvents(5);
    
    Species B("B", 100, max_ulim, SpeciesType::BULK, RSpeciesType::AVG);
    RSpecies& RB(B.getRSpecies());
    RSpeciesAvg* RBA = (RSpeciesAvg*)&RB;
    
    RBA->setNumEvents(5);
    
    //check first averages
    CHECK(100 == RA.getN());
    CHECK(100 == RB.getN());
    
    //check true n too
    CHECK(100 == RA.getTrueN());
    CHECK(100 == RB.getTrueN());
    
    // Using a reaction A->B to see if copy number n is correctly driven
    // (implicitly tests up() and down() methods)
    // we also have to explicity change tau so the averages
    // are computed correctly
    Reaction<1,1> rxn1 = {{&A,&B}, 10.0 };
    global_time = 1.0;
    rxn1.makeStep();
    global_time = 2.0;
    rxn1.makeStep();
    global_time = 3.0;
    rxn1.makeStep();
    global_time = 4.0;
    rxn1.makeStep();
    global_time = 5.0;
    rxn1.makeStep();
    global_time = 6.0;
    rxn1.makeStep();
    
    //new avg for A = 97.5, B = 102.5
    //true A = 94, B = 106
    CHECK(97.5 == RA.getN());
    CHECK(102.5 == RB.getN());
    
    CHECK(94 == RA.getTrueN());
    CHECK(106 == RB.getTrueN());
    
    //now reverse
    Reaction<1,1> rxn2 = {{&B,&A}, 10.0 };
    global_time = 7.0;
    rxn2.makeStep();
    global_time = 8.0;
    rxn2.makeStep();
    global_time = 9.0;
    rxn2.makeStep();
    global_time = 10.0;
    rxn2.makeStep();
    global_time = 11.0;
    rxn2.makeStep();
    global_time = 12.0;
    rxn2.makeStep();
    
    //new avg for A = 97.5, B = 102.5
    //true A = 100, B = 100
    CHECK(96.5 == RA.getN());
    CHECK(103.5 == RB.getN());
    
    CHECK(100 == RA.getTrueN());
    CHECK(100 == RB.getTrueN());
}


TEST_CASE("Reaction constructors", "[Reaction]") {
    Species A("A",  10, max_ulim, SpeciesType::BULK, RSpeciesType::REG);
    Species B("B",  10, max_ulim, SpeciesType::BULK, RSpeciesType::REG);
    Species C("C",  10, max_ulim, SpeciesType::BULK, RSpeciesType::REG);
    Species D("D",  10, max_ulim, SpeciesType::BULK, RSpeciesType::REG);
    Reaction<1,1> rxn1 = { {&A,&B}, 10.0 }; // A -> B
    Reaction<2,1> rxn2 = { {&A,&B,&C}, 10.0 }; // A+B -> C
    Reaction<1,2> rxn3 = { {&A,&B,&C}, 10.0 }; // A -> B + C
    Reaction<2,2> rxn4 = { {&A,&B,&C,&D}, 10.0 }; // A + B -> C + D

    CHECK(1 == rxn1.getM()); CHECK(1 == rxn1.getN());
    CHECK(2 == rxn2.getM()); CHECK(1 == rxn2.getN());
    CHECK(1 == rxn3.getM()); CHECK(2 == rxn3.getN());
    CHECK(2 == rxn4.getM()); CHECK(2 == rxn4.getN());

    // rxn1
    auto it1 = rxn1.rspecies();
    CHECK("A{Bulk}" == (*it1)->getFullName());
    // rxn2
    auto it2 = rxn2.rspecies();
    CHECK("A{Bulk}" == (*it2)->getFullName());
    ++it2;
    CHECK("B{Bulk}" == (*it2)->getFullName());
    CHECK("C{Bulk}" == (*(rxn2.rspecies()+rxn2.getM()))->getFullName());
    // rxn3
    auto it3 = rxn3.rspecies()+rxn2.getM();
    CHECK("A{Bulk}" == (*rxn3.rspecies())->getFullName());
    CHECK("B{Bulk}" == (*it2)->getFullName());
    ++it3;
    CHECK("C{Bulk}" == (*(rxn2.rspecies()+rxn2.getM()))->getFullName());
    // rxn4
    auto r_it4 = rxn4.rspecies();
    auto p_it4 = rxn4.rspecies()+rxn4.getM();
    CHECK("A{Bulk}" == (*r_it4)->getFullName());
    ++r_it4;
    CHECK("B{Bulk}" == (*r_it4)->getFullName());
    CHECK("C{Bulk}" == (*p_it4)->getFullName());
    ++p_it4;
    CHECK("D{Bulk}" == (*p_it4)->getFullName());

} 

TEST_CASE("RSpecies signaling", "[RSpecies]") {
    Species A("A",  8, max_ulim, SpeciesType::BULK, RSpeciesType::REG);
    A.connect(rspecies_callback);
    
    A.getRSpecies().emitSignal(0);
    CHECK(33 == A.getN());
}

#ifdef TRACK_DEPENDENTS
#ifdef TRACK_ZERO_COPY_N
TEST_CASE("Reaction dependents 1", "[Reaction]") {
    Species A("A",  10, max_ulim, SpeciesType::BULK, RSpeciesType::REG);
    Species B("B",  10, max_ulim, SpeciesType::BULK, RSpeciesType::REG);
    Species C("C",  10, max_ulim, SpeciesType::BULK, RSpeciesType::REG);
    RSpecies& RA(A.getRSpecies());
    Reaction<1,1> rxn1 = { {&A,&B}, 10.0 };
    Reaction<1,1> rxn3 = { {&A,&C}, 10.0 };
    Reaction<1,1> rxn2 = { {&B,&A}, 10.0 };
    
    rxn1.activateReaction();
    rxn2.activateReaction();
    rxn3.activateReaction();
    
    // note: we have three reactions, (1) A->B, (2) B->A, (3) A->C
    // (1) affects (2) and (3)
    // (2) affects (1) and (3) 
    // (3) affects (1)
    CHECK(2U == rxn1.dependents().size());
    CHECK(2U == rxn2.dependents().size());
    CHECK(1U == rxn3.dependents().size());
    CHECK(&rxn1 == *rxn3.dependents().begin());// (3) affects (1)
    
    // Testing passivateAssocReacts() and activateAssocReactions()
    // We run A->C 10 times, so A's copy number drops to 0.
    // Hence, the reaction (1) & (3) should be passivated.
    // In the context of this test, where we do not run a Gillespie-like algorithm,
    // this should lead to reaction (2) not having (1) or (3) as dependents - the
    // dependent count of (2) should become zero
    
    for (int i=0; i<10; ++i){
        rxn3.makeStep();
    }
    REQUIRE(0 == A.getN());
    CHECK( rxn1.isPassivated());
    CHECK(!rxn2.isPassivated());
    CHECK( rxn3.isPassivated());

    CHECK(1 == rxn1.dependents().size());
    CHECK(0 == rxn2.dependents().size());
    CHECK(0 == rxn3.dependents().size());
    
    // But the registration in RSpecies should not change.
    CHECK(A.getRSpecies().reactantReactions().size() == 2);
    CHECK(B.getRSpecies().reactantReactions().size() == 1);
    CHECK(C.getRSpecies().reactantReactions().size() == 0);

    CHECK(A.getRSpecies().productReactions().size() == 1);
    CHECK(B.getRSpecies().productReactions().size() == 1);
    CHECK(C.getRSpecies().productReactions().size() == 1);

    // Now let's activate (1) and (3) by moving the copy number of A from 0 to 1
    rxn2.makeStep();
    CHECK(1 == RA.getN());
    CHECK(!rxn1.isPassivated());
    CHECK(!rxn2.isPassivated());
    CHECK(!rxn3.isPassivated());
    CHECK(2 == rxn1.dependents().size());
    CHECK(2 == rxn2.dependents().size());
    CHECK(1 == rxn3.dependents().size());
}
#endif // of TRACK_ZERO_COPY_N
#endif // of TRACK_UPPER_COPY_N

TEST_CASE("Reaction dependents 2", "[Reaction]") {
    Species A("A",  10, max_ulim, SpeciesType::BULK, RSpeciesType::REG);
    Species B("B",  10, max_ulim, SpeciesType::BULK, RSpeciesType::REG);
    Species C("C",  10, max_ulim, SpeciesType::BULK, RSpeciesType::REG);
    Species D("D",  10, max_ulim, SpeciesType::BULK, RSpeciesType::REG);
    Species E("E",  10, max_ulim, SpeciesType::BULK, RSpeciesType::REG);
    Species F("F",  10, max_ulim, SpeciesType::BULK, RSpeciesType::REG);
    Reaction<1,1> rxn1 = { {&A,&B}, 10.0 };
    Reaction<1,1> rxn3 = { {&C,&D}, 10.0 };
    Reaction<1,1> rxn2 = { {&E,&F}, 10.0 };

    rxn1.activateReaction();
    rxn2.activateReaction();
    rxn3.activateReaction();
    
    // Let's add artificial dependencies by hand
    rxn1.registerNewDependent(&rxn2);
    CHECK(rxn1.dependents().find(&rxn2) != rxn1.dependents().end());
    rxn1.registerNewDependent(&rxn3);
    CHECK(rxn1.dependents().find(&rxn3) != rxn1.dependents().end());
    CHECK(2 == rxn1.dependents().size());
    rxn1.unregisterDependent(&rxn3);
    CHECK(1 == rxn1.dependents().size());
}

TEST_CASE("Reaction propensities", "[Reaction]") {
    Species A("A",  8,  max_ulim, SpeciesType::BULK, RSpeciesType::REG);
    Species B("B",  12, max_ulim, SpeciesType::BULK, RSpeciesType::REG);
    Species C("C",  14, max_ulim, SpeciesType::BULK, RSpeciesType::REG);
    
    Reaction<2,1> rxn = { {&A,&B,&C}, 3.14 }; // A+B -> C
    rxn.activateReaction();
    CHECK(Approx(8*12*3.14) == rxn.computePropensity());
    CHECK(8*12 == rxn.getProductOfReactants());
    
    Reaction<1,1> rxn2 = { {&A,&B}, 3.14 }; // A->B
    rxn2.activateReaction();
    CHECK(Approx(8*3.14) == rxn2.computePropensity());
    CHECK(8 == rxn2.getProductOfReactants());

}

TEST_CASE("Reaction signaling test", "[Reaction]") {
    Species A("A",  8,  max_ulim, SpeciesType::BULK, RSpeciesType::REG);
    Species B("B",  12, max_ulim, SpeciesType::BULK, RSpeciesType::REG);
    Reaction<1,1> rxn = { {&A,&B}, 3.14 }; // A->B
    
    // There are at least ways to set up callbacks: 1) functions; 2) functors; 3) lambda
    // functions. In terms of connection management, the return of sm.connect(...) can
    // be captured and used later to temporarily block the callback or permanently
    // disconnect it. See the boost documentation for signals2.
    
    // function<void (Reaction *)> rcb(reaction_callback);

    rxn.connect(reaction_callback);
    rxn.emitSignal();
    CHECK(Approx(0.05) == rxn.getRate());
    
    rxn.connect(ReactionCallback());
    // The callbacks will be called at the order of addition.
    rxn.emitSignal();
    CHECK(Approx(1.05) == rxn.getRate());

    rxn.clearSignaling();
    rxn.connect([](ReactionBase *r){r->setBareRate(3.05);});
    rxn.emitSignal();
    CHECK(Approx(3.05) == rxn.getRate());
}

TEST_CASE("Reaction cloning", "[Reaction]") {
    
    Compartment* C1 = new Compartment;
    Compartment* C2 = new Compartment;
    
    Species* ADiff1 = C1->addSpeciesUnique(std::make_unique<Species>("ADiff", 10, max_ulim, SpeciesType::DIFFUSING, RSpeciesType::REG));
    Species* ADiff2 = C2->addSpeciesUnique(std::make_unique<Species>("ADiff", 10, max_ulim, SpeciesType::DIFFUSING, RSpeciesType::REG));
    
    Species* BDiff1 = C1->addSpeciesUnique(std::make_unique<Species>("BDiff", 10, max_ulim, SpeciesType::DIFFUSING, RSpeciesType::REG));
    Species* BDiff2 = C2->addSpeciesUnique(std::make_unique<Species>("BDiff", 10, max_ulim, SpeciesType::DIFFUSING, RSpeciesType::REG));
    
    ReactionBase* r1 = C1->addInternal<Reaction,1,1>({ADiff1,BDiff1}, 100.0);

    r1->connect([](ReactionBase *r){r->setBareRate(9.0);});
    
    ///Clone, check if valid
    ReactionBase* r2 = r1->clone(C2->getSpeciesContainer());
    
    CHECK(1 == r2->getM());
    CHECK(1 == r2->getN());
    
    CHECK(r2->containsSpecies(ADiff2));
    CHECK(r2->containsSpecies(BDiff2));
    
    ///Check signal cloning
    r2->emitSignal();
    CHECK(9.0 == r2->getRate());
    
    ///Clone a reaction where not all species are in compartment
    Compartment* C3 = new Compartment;
    Species* CDiff3 = C3->addSpeciesUnique(std::make_unique<Species>("CDiff", 10, max_ulim, SpeciesType::DIFFUSING, RSpeciesType::REG));
    Species* ADiff3 = C3->addSpeciesUnique(std::make_unique<Species>("ADiff", 10, max_ulim, SpeciesType::DIFFUSING, RSpeciesType::REG));
    
    ReactionBase* r3 = C3->addInternal<Reaction,1,1>({ADiff3, CDiff3}, 100.0);
    ReactionBase* r4 = r3->clone(C1->getSpeciesContainer());
    
    ///Should keep cdiff3 in reactants
    CHECK(r4->containsSpecies(CDiff3));
    CHECK(r4->containsSpecies(ADiff1));
    
    ///Check reaction equality
    CHECK(r3->is_equal(*r4));
}

} // namespace medyan

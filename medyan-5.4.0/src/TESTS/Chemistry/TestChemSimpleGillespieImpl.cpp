

#include "catch2/catch.hpp"

#include "Chemistry/ChemSimpleGillespieImpl.h"
#include "Chemistry/ReactionDy.hpp"

namespace medyan {

TEST_CASE("ChemSimpleGillespieImpl tests", "[ChemSim]") {

    Species a{"A", 0, 1000000, SpeciesType::unspecified, RSpeciesType::REG};
    Species b{"B", 0, 1000000, SpeciesType::unspecified, RSpeciesType::REG};
    Species c{"C", 0, 1000000, SpeciesType::unspecified, RSpeciesType::REG};
    vector<Species*> reactants{ &a, &b };
    vector<Species*> products{ &c };
    ReactionDy reaction{reactants, products, ReactionType::REGULAR, 1.0};
    ChemSimpleGillespieImpl sim{};//resets global time to 0
    sim.addReaction(&reaction);
    sim.initialize();
    //sim.printReactions();
    
    SECTION("Test runSteps") {
        REQUIRE(sim.computeTotalA()==0);
        a.up();
        REQUIRE(sim.computeTotalA()==0);
        b.up();
        REQUIRE(sim.computeTotalA()==1.0);
        b.up();
        REQUIRE(sim.computeTotalA()==2.0);
        REQUIRE(sim.runSteps(1));
        REQUIRE(a.getN()==0);
        REQUIRE(b.getN()==1);
        REQUIRE(c.getN()==1);
        REQUIRE(sim.computeTotalA()==0);
        REQUIRE(sim.getTime()>0);
    }

    SECTION("Test run") {
        REQUIRE(sim.run(1.0));
        REQUIRE(c.getN()==0);
        REQUIRE(a.getN()==0);
        REQUIRE(b.getN()==0);
        REQUIRE(sim.getTime() == Approx(1.0));
        a.up();
        b.up();
        REQUIRE(sim.run(100.0));
        REQUIRE(sim.getTime() == Approx(101.0));
        REQUIRE(a.getN()==0);
        REQUIRE(b.getN()==0);
        REQUIRE(c.getN()==1);
    }

}

} // namespace medyan

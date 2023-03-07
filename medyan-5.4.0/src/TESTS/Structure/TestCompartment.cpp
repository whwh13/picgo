
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

#include "Chemistry/Species.h"
#include "Controller/GController.h"
#include "Structure/CompartmentGrid.h"
#include "SysParams.h"

using namespace medyan;

//Testing neighbor capabilities
TEST_CASE("Compartment neighbors", "[Compartment]") {

    Compartment c1;     c1._ID = 100;
    Compartment c2(c1); c2._ID = 200;
    c1.setNeighbor(c2._ID, 0);
    c2.setNeighbor(c1._ID, 0);
    
    CHECK(1 == c1.numberOfNeighbours());
    CHECK(1 == c2.numberOfNeighbours());
    
    Compartment c3(c2); c3._ID = 300;
    c1.setNeighbor(c3._ID, 1);
    
    CHECK(2 == c1.numberOfNeighbours());
    CHECK(0 == c3.numberOfNeighbours());

    // Clear neighbor in direction 0 (which was c2).
    c1.resetNeighbor(0);
    
    CHECK(1 == c1.numberOfNeighbours());
}


#ifdef TESTING


//Testing species and reaction generation
TEST(Compartment, SpeciesAndReactions){
    
    Compartment C1;
    C1.setAsActive();
    
    Species *actin = C1.addSpeciesDiffusing("Actin",99U);
    C1.setDiffusionRate(actin,2000);
    Species *profilin = C1.addSpeciesDiffusing("Profilin",29U);
    C1.setDiffusionRate(profilin,2000);
    Species *arp23 = C1.addSpeciesDiffusing("Arp2/3",33U);
    C1.setDiffusionRate(arp23,2000);
    
    Compartment C2(C1);
    
    EXPECT_EQ(2000,C1.getDiffusionRate("Actin"));
    EXPECT_EQ(2000,C2.getDiffusionRate("Actin"));
    
    C1.setDiffusionRate(actin, 9000);
    EXPECT_EQ(9000, C1.getDiffusionRate("Actin"));
    EXPECT_EQ(2000, C2.getDiffusionRate("Actin"));
    
    
    C1.addNeighbour(&C2);
    C2.addNeighbour(&C1);
    C1.generateAllDiffusionReactions();
    C2.generateAllDiffusionReactions();
    
    EXPECT_EQ(3, C1.numberOfSpecies());
    EXPECT_EQ(3, C2.numberOfSpecies());
    
    EXPECT_EQ(0, C1.numberOfInternalReactions());
    EXPECT_EQ(3, C1.numberOfReactions());
    EXPECT_EQ(0, C2.numberOfInternalReactions());
    EXPECT_EQ(3, C2.numberOfReactions());
}


//Testing species transfers
TEST(Compartment, TransferSpecies) {
    
    //simple transfer
    Compartment C1;
    C1.setAsActive();
    
    Species *actin = C1.addSpeciesDiffusing("Actin",99U);
    C1.setDiffusionRate(actin,2000);
    Species *profilin = C1.addSpeciesDiffusing("Profilin",29U);
    C1.setDiffusionRate(profilin,2000);
    Species *arp23 = C1.addSpeciesDiffusing("Arp2/3",33U);
    C1.setDiffusionRate(arp23,2000);
    
    Compartment *C2 = C1.clone();
    
    C1.addNeighbour(C2);
    C2->addNeighbour(&C1);
    
    //transfer
    C1.transferSpecies();
    
    EXPECT_EQ(198, C2->findSpeciesByName("Actin")->getN());
    EXPECT_EQ(0, C1.findSpeciesByName("Actin")->getN());
    
    EXPECT_EQ(58, C2->findSpeciesByName("Profilin")->getN());
    EXPECT_EQ(0, C1.findSpeciesByName("Profilin")->getN());
    
    EXPECT_EQ(66, C2->findSpeciesByName("Arp2/3")->getN());
    EXPECT_EQ(0, C1.findSpeciesByName("Arp2/3")->getN());
    
    //now, add third
    Compartment* C3 = C1.clone();
    
    C2->addNeighbour(C3);
    C2->transferSpecies();
    
    EXPECT_EQ(198, C1.findSpeciesByName("Actin")->getN() +
                   C3->findSpeciesByName("Actin")->getN());
    EXPECT_EQ(0, C2->findSpeciesByName("Actin")->getN());
    
    EXPECT_EQ(58, C1.findSpeciesByName("Profilin")->getN() +
                  C3->findSpeciesByName("Profilin")->getN());
    EXPECT_EQ(0, C2->findSpeciesByName("Profilin")->getN());
    
    EXPECT_EQ(66, C1.findSpeciesByName("Arp2/3")->getN() +
                  C3->findSpeciesByName("Arp2/3")->getN());
    EXPECT_EQ(0, C2->findSpeciesByName("Arp2/3")->getN());
}

//Testing neighbors, species and reaciton generation
TEST(CompartmentContainerTest, Main) {

    GeoParams gparams;
    gparams.compartmentSizeX = 100.0;
    gparams.compartmentSizeY = 100.0;
    gparams.compartmentSizeZ = 100.0;
    
    gparams.NX = 50;
    gparams.NY = 50;
    gparams.NZ = 50;
    
    gparams.monomerSize = 2.7;
    gparams.cylinderSize = 27.0;
    
    int _numSpecies; ///for testing

    GController g;
    CompartmentGrid* grid = g.initializeGrid(gparams);
    
    ///activate all compartments for diffusion
    grid->setAllAsActive();
    
    Compartment &Cproto = grid->getProtoCompartment();
    
    Species *M1 = Cproto.addSpeciesDiffusing("Myosin",1U);
    Cproto.setDiffusionRate(M1,2000);
    Species *M2 = Cproto.addSpeciesDiffusing("Fascin",6U);
    Cproto.setDiffusionRate(M2,2000);
    
    ReactionBase* r1 = new Reaction<1,1>({M2, M1}, 90.9, true);
    ReactionBase* r2 = new Reaction<1,1>({M1, M2}, 40.2, true);
    
    Cproto.addInternalReaction(r1);
    Cproto.addInternalReaction(r2);

    ///initialize all compartments with species
    for(auto& C : grid->getCompartments())
        *C = Cproto;
    
    ///Generate all diffusion reactions
    for(auto& C : grid->getCompartments())
        C->generateAllDiffusionReactions();
    
    _numSpecies = grid->countSpecies();
    
    Compartment *C1 = &GController::getCompartment(vector<size_t>{10U,10U,10U});
    
    EXPECT_EQ(6, C1->numberOfNeighbours());
    EXPECT_EQ(2, C1->numberOfSpecies());
    EXPECT_EQ(14, C1->numberOfReactions());
    
    Compartment*C2 = &GController::getCompartment(vector<size_t>{10U,11U,10U});
    C1->removeNeighbour(C2);
    
    EXPECT_EQ(5, C1->numberOfNeighbours());
    
}

#endif //TESTING


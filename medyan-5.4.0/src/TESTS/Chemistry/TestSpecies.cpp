
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

#ifdef TESTING

// Note: This test omits many functions of Species that interact
// with Reaction objects. Separate tests weill cover those methods.

#include <catch2/catch.hpp>

#include "common.h"
#include "Species.h"
#include "SpeciesContainer.h"

using namespace medyan;

TEST_CASE("SpeciesNamesDB", "[Species]") {
    
    //basic test
    SpeciesNamesDB::clear();
    int y = SpeciesNamesDB::stringToInt("Arp2/3");
    EXPECT_EQ(0,y);
    string x = SpeciesNamesDB::intToString(y);
    EXPECT_EQ("Arp2/3",x);
    
    EXPECT_THROW(SpeciesNamesDB::intToString(1), out_of_range);

    y = SpeciesNamesDB::stringToInt("G-Actin");
    EXPECT_EQ(1,y);
    x = SpeciesNamesDB::intToString(y);
    EXPECT_EQ("G-Actin",x);
    EXPECT_NO_THROW(SpeciesNamesDB::intToString(1));
    
    //testing unique name generator
    string a1 = SpeciesNamesDB::genUniqueFilName("Actin");
    string a2 = SpeciesNamesDB::genUniqueFilName("Actin");
    string a3 = SpeciesNamesDB::genUniqueFilName("Actin");
    
    EXPECT_FALSE(a1 == a2);
    EXPECT_FALSE(a2 == a3);
    EXPECT_FALSE(a1 == a3);
}

TEST_CASE("Species constructors", "[Species]") {

    //Default CTor
    SpeciesBulk X;
    EXPECT_EQ(0,X.getN());
    EXPECT_EQ("",X.getName());

    //CTor 1
    SpeciesBulk A{"Arp2/3",25};
    EXPECT_EQ(25,A.getN());
    EXPECT_EQ("Arp2/3",A.getName());
    
    //Copy Ctor
    SpeciesBulk B(A);
    EXPECT_EQ(A,B);
    
    //Move Assignment
    SpeciesBulk D{"G-Actin",300};
    D = [] () {SpeciesBulk Y("G-Actin",39); return Y;}();
    EXPECT_EQ(39,D.getN());
    EXPECT_EQ("G-Actin",D.getName());
    
    //Assignment operator
    SpeciesBulk F{"Myosin-X",112};
    F=A;
    EXPECT_EQ(A,F);

    
}

TEST_CASE("Species copy number", "[Species]") {
    SpeciesBulk B{"Arp2/3",25};
    B.setN(32);
    EXPECT_EQ(32,B.getN());
}


TEST_CASE("Species formatting", "[Species]") {    
    SpeciesBulk B{"Arp2/3",25};

    ostringstream outStream;
    
    outStream << B;
    string res = outStream.str();
    string b_fullname = B.getFullName() + "[" + to_string(B.getN()) + "]";
    EXPECT_EQ(b_fullname,res);
}

TEST_CASE("Species vector", "[Species]") {
    vector<SpeciesBulk> vsb;
    vsb.push_back({"G-Actin",133});
    EXPECT_EQ(133,vsb.back().getN());
    EXPECT_EQ("G-Actin",vsb.back().getName());
    
    vsb.emplace_back("Arp2/3",244);
    EXPECT_EQ("Arp2/3",vsb.back().getName());

    EXPECT_NE(vsb[0],vsb[1]);
    
    vsb.push_back({"Motor",800});
    EXPECT_EQ("Motor",vsb[2].getName());
    RSpecies &rs_before = vsb[2].getRSpecies();
    
    for(int i=0; i<10; ++i){
        vsb.push_back({"Motor",900});
    }
    
    // To make sure that move should be enabled:
    EXPECT_TRUE(is_nothrow_move_constructible<SpeciesBulk>::value);
    EXPECT_TRUE(is_copy_constructible<SpeciesBulk>::value);
    
    // Now checking if Species have been moved and not
    // copied by the vector<Species> reallocators
    RSpecies &rs_after = vsb[2].getRSpecies();
    EXPECT_EQ(&rs_before,&rs_after);
    // This makes sure that the move constructors were involved during the internal
    // reallocation of the vector<Species>, and RSpecies pointer was conserved
}

TEST_CASE("Species pointer container", "[SpeciesContainer]") {
    SpeciesPtrContainerVector scv;
    
    Species *x0 = scv.addSpecies<SpeciesBulk>("Profilin",31);
    Species *x1 = scv.addSpecies<SpeciesDiffusing>("Actin",99);
    Species *x2 = scv.addSpecies<SpeciesDiffusing>("Arp2/3",11);
    Species *x3 = scv.addSpecies<SpeciesDiffusing>("Capping",22);
    
    EXPECT_EQ(x0, scv.findSpeciesByIndex(0));
    EXPECT_EQ(x1, scv.findSpeciesByIndex(1));
    EXPECT_EQ(x2, scv.findSpeciesByIndex(2));
    EXPECT_EQ(x3, scv.findSpeciesByIndex(3));

    EXPECT_EQ("Capping", scv.findSpeciesByIndex(3)->getName());
    
    scv.removeSpecies(x0);
    EXPECT_EQ(x3, scv.findSpeciesByIndex(2));
    
    scv.removeSpecies("Actin");
    EXPECT_EQ(x3, scv.findSpeciesByIndex(1));
    Species *y = scv.findSpeciesByName("Arp2/3");
    EXPECT_EQ(x2, y);



    SpeciesPtrContainerVector cont;
    Species *arp23 = cont.addSpecies(new SpeciesDiffusing("Arp2/3",55));
    Species *actin = cont.addSpecies<SpeciesDiffusing>("Actin",44);
    Species *profilin = cont.addSpeciesUnique(make_unique<SpeciesDiffusing>("Profilin",33));
    
    EXPECT_EQ(3U, cont.size());
    EXPECT_EQ(3U, cont.species().size());
    
    EXPECT_EQ(arp23,cont.findSpeciesByIndex(0));
    EXPECT_EQ(actin,cont.findSpeciesByIndex(1));
    EXPECT_EQ(profilin,cont.findSpeciesByIndex(2));
    
    EXPECT_EQ(arp23,cont.findSpeciesByName("Arp2/3"));
    EXPECT_EQ(actin,cont.findSpeciesByName("Actin"));
    EXPECT_EQ(profilin,cont.findSpeciesByName("Profilin"));
    
    SpeciesDiffusing another_profilin("Profilin",33);
    EXPECT_EQ(profilin,cont.findSimilarSpecies(another_profilin));
    
    auto molecule = another_profilin.getMolecule();
    EXPECT_EQ(profilin,cont.findSpeciesByMolecule(molecule));
    
    EXPECT_EQ(1,cont.areAllSpeciesUnique());
    
    cont.addSpecies(new SpeciesDiffusing("Arp2/3",66));
    
    EXPECT_EQ(2, cont.removeSpecies("Arp2/3"));
    EXPECT_EQ(2U, cont.size());
    
    EXPECT_EQ(1, cont.removeSpecies(actin));
    EXPECT_EQ(1U, cont.size());
}

TEST_CASE("Species vector container", "[SpeciesContainer]") {
    SpeciesContainerVector<SpeciesDiffusing> scv;
    
    size_t profilin = scv.addSpecies("Profilin",species_copy_t(31));
    size_t actin = scv.addSpecies("Actin",species_copy_t(99));
    size_t arp23 = scv.addSpecies("Arp2/3",species_copy_t(11));
    size_t capping = scv.addSpecies("Capping",species_copy_t(22));
    
    EXPECT_EQ(size_t(0), profilin);
    EXPECT_EQ(size_t(1), actin);
    EXPECT_EQ(size_t(2), arp23);
    EXPECT_EQ(size_t(3), capping);
    
    EXPECT_EQ("Profilin", scv.findSpecies(0).getName());
    EXPECT_EQ("Actin", scv.findSpecies(1).getName());
    EXPECT_EQ("Arp2/3", scv.findSpecies(2).getName());
    EXPECT_EQ("Capping", scv.findSpecies(3).getName());
    
    scv.removeSpecies("Profilin");
    EXPECT_EQ("Capping", scv.findSpecies(2).getName());
    
    scv.removeSpecies("Actin");
    EXPECT_EQ("Capping", scv.findSpecies(1).getName());
    
    Species &y = scv.findSpecies("Arp2/3");
    EXPECT_EQ("Arp2/3", y.getName());
    
    size_t index = scv.findSpeciesIndex("Arp2/3");
    EXPECT_EQ(0U,index);



    SpeciesContainerVector<SpeciesDiffusing> cont;
    size_t arp23 = cont.addSpecies("Arp2/3",55);
    size_t actin = cont.addSpecies("Actin",44);
    size_t profilin = cont.addSpecies("Profilin",33);
    
    EXPECT_EQ(3U, cont.size());
    EXPECT_EQ(3U, cont.species().size());
    
    EXPECT_EQ("Arp2/3",cont.findSpecies(arp23).getName());
    EXPECT_EQ("Actin",cont.findSpecies(actin).getName());
    EXPECT_EQ("Profilin",cont.findSpecies(profilin).getName());
    
    EXPECT_EQ("Arp2/3",cont.findSpecies("Arp2/3").getName());
    EXPECT_EQ("Actin",cont.findSpecies("Actin").getName());
    EXPECT_EQ("Profilin",cont.findSpecies("Profilin").getName());
    
    EXPECT_EQ(profilin,cont.findSpeciesIndex("Profilin"));
    
    
    SpeciesDiffusing another_profilin("Profilin",33);
    EXPECT_EQ("Profilin",cont.findSimilarSpecies(another_profilin).getName());
    
    auto molecule = another_profilin.getMolecule();
    EXPECT_EQ("Profilin",cont.findSpeciesByMolecule(molecule).getName());
    
    EXPECT_EQ(1,cont.areAllSpeciesUnique());
    
    cont.addSpecies("Arp2/3",66);
    
    EXPECT_EQ(2, cont.removeSpecies("Arp2/3"));
    EXPECT_EQ(2U, cont.size());
    
    EXPECT_EQ(1, cont.removeSpecies(actin));
    EXPECT_EQ(1U, cont.size());
    
}

#endif //TESTING

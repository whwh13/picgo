
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

#include "Controller/GController.h"
#include "Structure/SubSystem.h"
#include "SysParams.h"

TEST_CASE("Geometry controller", "[Geometry]") {

    using namespace medyan;

    GeoParams gparams;
    gparams.compartmentSizeX = 10.0;
    gparams.compartmentSizeY = 20.0;
    gparams.compartmentSizeZ = 100.0;
    gparams.NX = 15;
    gparams.NY = 5;
    gparams.NZ = 10;

    SubSystem s;
    GController g(&s);
    g.initializeGrid(gparams);

    CHECK_THROWS(GController::getCompartment(medyan::Vec<3, floatingpoint>{60.0,  50.0,  1050.0}));
    CHECK_THROWS(GController::getCompartment(medyan::Vec<3, floatingpoint>{200.0, 50.0,  900.0}));
    CHECK_THROWS(GController::getCompartment(medyan::Vec<3, floatingpoint>{100.0, 110.0, 900.0}));
    
    CHECK(GController::getCompartment(vector<size_t>{0,0,0}) == GController::getCompartment(medyan::Vec<3, floatingpoint>{5.0, 5.0, 5.0}));
    CHECK(GController::getCompartment(vector<size_t>{0,1,0}) == GController::getCompartment(medyan::Vec<3, floatingpoint>{5.0, 25.0, 5.0}));
    CHECK(GController::getCompartment(vector<size_t>{0,1,1}) == GController::getCompartment(medyan::Vec<3, floatingpoint>{5.0, 30.0, 190.0}));
}

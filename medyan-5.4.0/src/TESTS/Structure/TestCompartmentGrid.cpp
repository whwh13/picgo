#include <catch2/catch.hpp>

#include "Structure/CompartmentGrid.h"
#include "SysParams.h"

namespace medyan {

TEST_CASE("Compartment grid", "[Compartment]") {
    using namespace std;

    GeoParams geo;
    geo.compartmentSizeX = 10.0;
    geo.compartmentSizeY = 20.0;
    geo.compartmentSizeZ = 100.0;
    geo.NX = 15;
    geo.NY = 5;
    geo.NZ = 10;

    CompartmentGrid g(geo);

    // Check geometry parameters.
    CHECK(g.shape[0] == 15);
    CHECK(g.shape[1] == 5);
    CHECK(g.shape[2] == 10);

    CHECK(g.compartmentLengths[0] == 10.0);
    CHECK(g.compartmentLengths[1] == 20.0);
    CHECK(g.compartmentLengths[2] == 100.0);

    CHECK(g.compartmentAreas[0] == Approx(2000.0));
    CHECK(g.compartmentAreas[1] == Approx(1000.0));
    CHECK(g.compartmentAreas[2] == Approx(200.0));

    CHECK(g.compartmentVolume == Approx(20000.0));

    CHECK(g.gridLengths[0] == Approx(150.0));
    CHECK(g.gridLengths[1] == Approx(100.0));
    CHECK(g.gridLengths[2] == Approx(1000.0));

    CHECK(g.gridVolume == Approx(1.5e7));

    // Check compartment list and indices.
    REQUIRE(g.compartmentList.size() == geo.NX * geo.NY * geo.NZ);

    CHECK(g.getCompartmentIndex3({   2.0,  2.0,  40.0 }) == array<Index,3> {  0, 0, 0 });
    CHECK(g.getCompartmentIndex3({ 148.0, 98.0, 960.0 }) == array<Index,3> { 14, 4, 9 });
    CHECK_THROWS(g.getCompartmentIndex3({ -1.0, 2.0, 40.0 }));
    CHECK_THROWS(g.getCompartmentIndex3({ 2.0, 200.0, 40.0 }));

    CHECK(g.getCompartmentIndex(array<Index,3> { 0, 0, 0 }) == 0);
    CHECK(g.getCompartmentIndex(array<Index,3> { 1, 0, 0 }) == 1);
    CHECK(g.getCompartmentIndex(array<Index,3> { 14, 4, 9 }) == geo.NX * geo.NY * geo.NZ - 1);
}

} // namespace medyan

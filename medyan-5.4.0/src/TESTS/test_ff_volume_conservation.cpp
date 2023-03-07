
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2015-2016)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifdef TESTING

#  define DO_THIS_FF_VOLUME_CONSERVATION_TEST
#  ifdef DO_THIS_FF_VOLUME_CONSERVATION_TEST

#    include "gtest/gtest.h"
#    include "test_public.h"

#    include <random>

#    include "common.h"
#    include "MathFunctions.h"
using namespace mathfunc;
#    include "Rand.h"

#    include "Controller/GController.h"
#    include "SubSystem.h"

#    include "Membrane.hpp"
#    include "Vertex.hpp"

#    include "Mechanics/ForceField/VolumeConservation/VolConsrvInteractions.hpp"
#    include "Mechanics/ForceField/VolumeConservation/VolConsrvMembrane.hpp"
#    include "VolConsrvMembraneHarmonic.hpp"

namespace {
    using VertexData = tuple<array<double, 3>, vector<size_t>>;
    using MembraneData = vector<VertexData>;

    class VolumeConservationFFTest: public ::testing::Test {
    protected:
        double radius;
        SubSystem s;
        MembraneData memData;
        Membrane *m;

		VolumeConservationFFTest():
            radius(100),
            memData(test_public::membraneDataOctahedron({2*radius, 2*radius, 2*radius}, radius))
        {
            test_public::quickSetupPlayground(&s);

            SysParams::GParams.cylinderNumMon.resize(1, 3);

            SysParams::MParams.bulkModulus = 2150; // BM of water in pN/nm^2

            m = new Membrane(&s, 0, memData);
            m->addToSubSystem();
        }
        ~VolumeConservationFFTest() {
            m->removeFromSubSystem();
            delete m;
        }

    };

    void recordCoordinate(Membrane *m) {
        for(Vertex* it: m->getVertexVector()) it->coordinateP = it->coordinate;
    }
    void resetCoordinate(Membrane *m) {
        for(Vertex* it: m->getVertexVector()) it->coordinate = it->coordinateP;
    }
    void resetForce(Membrane *m) {
        for(Vertex* it: m->getVertexVector()) it->force.assign(3, 0);
    }
    void assignRandomForceAuxP(Membrane* m, double sigma) {
		normal_distribution<> nd(0, sigma);

		for(Vertex* it: m->getVertexVector()) {
			for(double& eachForce: it->forceAuxP) { // forceAuxP is used here simply because it is an empty container.
				eachForce = nd(Rand::engFixed);
			}
		}
    }
    void moveAlongForceAuxP(Membrane* m, double d) {
        for(Vertex* it: m->getVertexVector()) {
            for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
                it->coordinate[coordIdx] += it->forceAuxP[coordIdx] * d;
            }
        }
    }
    void resizeEqVolume(Membrane *m, double ratio) {
        m->mMembrane.eqVolume *= ratio;
    }
}

TEST_F(VolumeConservationFFTest, Force) {

    VolumeConservationMembrane<VolumeConservationMembraneHarmonic> vcmHarmonic;

    assignRandomForceAuxP(m, radius/500);
    resizeEqVolume(m, 0.95); // This is to shrink the equilibrium volume by a little bit to avoid zero forces.
    moveAlongForceAuxP(m, 2.5); // Distort the structure a little
    recordCoordinate(m);

    // Compare the results with force predictions
    // U(x+h) - U(x-h) = dotProduct(2h, dU/dx) = -dotProduct(2h, F)

    resetCoordinate(m);
    resetForce(m);
    m->updateGeometry(true);
    vcmHarmonic.computeForces();

    moveAlongForceAuxP(m, 1.0);
    m->updateGeometry(true);
    double U1 = vcmHarmonic.computeEnergy(0.0);

    resetCoordinate(m);
    moveAlongForceAuxP(m, -1.0);
    m->updateGeometry(true);
    double U2 = vcmHarmonic.computeEnergy(0.0);

    double exDiff = 0.0;
    for(Vertex* v: m->getVertexVector())
        exDiff -= 2 * dotProduct(v->forceAuxP, v->force);

    EXPECT_NEAR(U1 - U2, exDiff, abs(exDiff / 1000))
        << vcmHarmonic.getName() << " force not working properly.";

}


#  endif //DO_THIS_FF_VOLUME_CONSERVATION_TEST
#endif //TESTING


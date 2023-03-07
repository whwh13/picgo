
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

#  define DO_THIS_FF_TRIANGLE_CYLINDER_VOLUME_TEST
#  ifdef DO_THIS_FF_TRIANGLE_CYLINDER_VOLUME_TEST

#    include "gtest/gtest.h"

#    include "common.h"
#    include "test_public.h"
#    include "MathFunctions.h"
using namespace mathfunc;
#    include "Rand.h"

#    include "Controller/GController.h"
#    include "SubSystem.h"
#    include "Component.h"

#    include "Edge.hpp"
#    include "Triangle.hpp"
#    include "Cylinder.h"

#    include "TriangleBeadExclVolume.hpp"
#    include "TriangleBeadExclVolRepulsion.hpp"

namespace {

    class TriangleCylinderVolumeFFTest: public ::testing::Test {
    protected:
        double radius;
        double height;
        test_public::CompositeDummy dummyParent;
        SubSystem s;
        
        Triangle *t;
        array<Vertex*, 3> tv;
        array<Edge*, 3> te;
        Cylinder *c;
        array<Bead*, 2> cb;

		TriangleCylinderVolumeFFTest():
            dummyParent(0),
            radius(100), height(10) {
            
            test_public::quickSetupPlayground(&s);
            test_public::quickSetupChem(&s);

            SysParams::GParams.cylinderNumMon.resize(1, 3);
			SysParams::GParams.cylinderSize.resize(1, 1.0);
			SysParams::GParams.monomerSize.resize(1, 1.0);
            SysParams::CParams.bindingSites.resize(1); // For CCylinder use

            SysParams::MParams.triangleBeadVolume.k = 1725;
            SysParams::MParams.triangleBeadVolume.cutoffMech = 15;
            SysParams::MParams.triangleBeadVolume.cutoff = 15;

            double trans = radius + height; // To ensure that all coordinates of all the beads are greater than 0.

            // Add a triangle
            tv[0] = new Vertex({trans + radius, trans, trans}, &dummyParent, 0);
            tv[1] = new Vertex({trans - radius/2, trans + sqrt(3)*radius/2, trans}, &dummyParent, 0);
            tv[2] = new Vertex({trans - radius/2, trans - sqrt(3)*radius/2, trans}, &dummyParent, 0);
            for(Vertex* eachV: tv) eachV->addToSubSystem();
            te[0] = new Edge(&dummyParent, tv[0], tv[1]);
            te[1] = new Edge(&dummyParent, tv[1], tv[2]);
            te[2] = new Edge(&dummyParent, tv[2], tv[0]);
            for(Edge* eachE: te) eachE->addToSubSystem();
            t = new Triangle(&dummyParent, tv[0], tv[1], tv[2]);
            t->addToSubSystem();
            t->getEdges() = te;

            // Add a cylinder
            cb[0] = new Bead({trans + radius/20, trans + radius/20, trans + height}, &dummyParent, 0);
            cb[1] = new Bead({trans - radius/20, trans - radius/20, trans + height}, &dummyParent, 0);
            for(Bead* eachB: cb) eachB->addToSubSystem();
            c = new Cylinder(&dummyParent, cb[0], cb[1], 0, 0);
            c->addToSubSystem();
        }
        ~TriangleCylinderVolumeFFTest() {
            SysParams::GParams.cylinderNumMon.resize(0);
			SysParams::GParams.cylinderSize.resize(0);
			SysParams::GParams.monomerSize.resize(0);
            SysParams::CParams.bindingSites.resize(0);

            // Remove the triangle
            t->removeFromSubSystem();
            for(Vertex* eachV: tv) {
                eachV->removeFromSubSystem();
            }
            for(Edge* eachE: te) {
                eachE->removeFromSubSystem();
            }

            // Remove the cylinder
            c->removeFromSubSystem();
            for(Bead* eachB: cb) {
                eachB->removeFromSubSystem();
            }
        }

        void recordCoordinate() {
            for(Vertex* it: tv) it->coordinateP = it->coordinate;
            for(Bead* it: cb) it->coordinateP = it->coordinate;
        }
        void resetCoordinate() {
            for(Vertex* it: tv) it->coordinate = it->coordinateP;
            for(Bead* it: cb) it->coordinate = it->coordinateP;
        }
        void resetForce() {
            for(Vertex* it: tv) it->force.assign(3, 0);
            for(Bead* it: cb) it->force.assign(3, 0);
        }
        void assignRandomForceAuxP(double sigma) {
		normal_distribution<> nd(0, sigma);

            for(Vertex* it: tv) {
                for(double& eachForce: it->forceAuxP) { // forceAuxP is used here simply because it is an empty container.
                    eachForce = nd(Rand::engFixed);
                }
            }
            for(Bead* it: cb) {
                for(double& eachForce: it->forceAuxP) { // forceAuxP is used here simply because it is an empty container.
                    eachForce = nd(Rand::engFixed);
                }
            }
        }
        void moveAlongForceAuxP(double d) {
            for(Vertex* it: tv) {
                for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
                    it->coordinate[coordIdx] += it->forceAuxP[coordIdx] * d;
                }
            }
            for(Bead* it: cb) {
                for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
                    it->coordinate[coordIdx] += it->forceAuxP[coordIdx] * d;
                }
            }
        }

    };

}

TEST_F(TriangleCylinderVolumeFFTest, Force) {
    TriangleBeadExclVolume<TriangleBeadExclVolRepulsion> tcv;

    assignRandomForceAuxP(radius/200);
    recordCoordinate();

    // Compare the results with force predictions
    // U(x+h) - U(x-h) = dotProduct(2h, dU/dx) = -dotProduct(2h, F)

    resetCoordinate();
    resetForce();
    // Update position for neighbor list search
    t->updatePosition();
    c->updatePosition();
    // Update neighbor list
    tcv.getNeighborList()->addDynamicNeighbor(t);
	for (Edge* eachE : te) eachE->getGEdge()->calcLength();
	t->getGTriangle()->calcArea();
	tcv.computeForces();

    // Don't bother updating neighbor list here
    moveAlongForceAuxP(1.0);
    for(Edge* eachE: te) eachE->getGEdge()->calcLength();
    t->getGTriangle()->calcArea();
    double U1 = tcv.computeEnergy(0.0);

    resetCoordinate();
    moveAlongForceAuxP(-1.0);
	for (Edge* eachE : te) eachE->getGEdge()->calcLength();
	t->getGTriangle()->calcArea();
	double U2 = tcv.computeEnergy(0.0);

    double exDiff = 0.0;
    for(Vertex* v: tv)
        exDiff -= 2 * dotProduct(v->forceAuxP, v->force);
    for(Bead* b: cb)
        exDiff -= 2 * dotProduct(b->forceAuxP, b->force);

    EXPECT_NEAR(U1 - U2, exDiff, abs(exDiff / 1000))
        << tcv.getName() << " force not working properly.";

}


#  endif // DO_THIS_FF_TRIANGLE_CYLINDER_VOLUME_TEST
#endif // TESTING


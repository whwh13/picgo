
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

#  define DO_THIS_FF_MEMBRANE_TEST
#  ifdef DO_THIS_FF_MEMBRANE_TEST

#    include "gtest/gtest.h"

#    include <random>

#    include "common.h"
#    include "MathFunctions.h"
using namespace mathfunc;
#    include "Rand.h"

#    include "Controller/GController.h"
#    include "SubSystem.h"

#    include "Membrane.hpp"
#    include "Vertex.hpp"
#    include "MVoronoiCell.h"
#    include "Triangle.hpp"
#    include "MTriangle.hpp"

#    include "Interactions.hpp"
#    include "MembraneStretching.hpp"
#    include "MembraneStretchingImpl.hpp"
#    include "Bending.hpp"
#    include "MembraneBendingHelfrich.hpp"

namespace {
    using VertexData = tuple<array<double, 3>, vector<size_t>>;
    using MembraneData = vector<VertexData>;

    MembraneData membraneDataOctahedron(double radius) {
        double c = 2 * radius;

        return MembraneData{
            VertexData({c, c, c+radius}, {1, 2, 3, 4}),
            VertexData({c+radius, c, c}, {0, 4, 5, 2}),
            VertexData({c, c+radius, c}, {0, 1, 5, 3}),
            VertexData({c-radius, c, c}, {0, 2, 5, 4}),
            VertexData({c, c-radius, c}, {0, 3, 5, 1}),
            VertexData({c, c, c-radius}, {4, 3, 2, 1})
        };
    }

    class MembraneFFTest: public ::testing::Test {
    protected:
        double radius;
        SubSystem s;
        MembraneData memData;
        Membrane *m;

		MembraneFFTest(): radius(100), memData(membraneDataOctahedron(radius)) {
            SysParams::GParams.compartmentSizeX = 1e10;
            SysParams::GParams.compartmentSizeY = 1e10;
            SysParams::GParams.compartmentSizeZ = 1e10;
            
            SysParams::GParams.NX = 1;
            SysParams::GParams.NY = 1;
            SysParams::GParams.NZ = 1;
            
            SysParams::GParams.nDim = 3;

            GController g(&s); // Dummy variable to initialize the compartments
            g.initializeGrid();

            SysParams::GParams.cylinderNumMon.resize(1, 3);

            SysParams::MParams.memAreaK.resize(1, 400);
            SysParams::MParams.memEqAreaFactor.resize(1, 1.0);
            SysParams::MParams.MemBendingK.resize(1, 100);
            SysParams::MParams.MemEqCurv.resize(1, 0);

            m = new Membrane(&s, 0, memData);
            m->addToSubSystem();
        }
        ~MembraneFFTest() {
            SysParams::GParams.cylinderNumMon.resize(0);

            SysParams::MParams.memAreaK.resize(0);
            SysParams::MParams.memEqAreaFactor.resize(0);
            SysParams::MParams.MemBendingK.resize(0);
            SysParams::MParams.MemEqCurv.resize(0);

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
    void resizeEqArea(Membrane *m, double ratio) {
        for(Triangle* t: m->getTriangleVector()) {
            t->mTriangle.eqArea *= ratio;
        }
    }
}

TEST_F(MembraneFFTest, CompareStretchingEnergy) {

    // Check that the two different methods give the same stretching energy.
    MembraneStretching<MembraneStretchingHarmonic> mSTriangle;
    MembraneStretching<MembraneStretchingVoronoiHarmonic> mSVoronoi;

    resizeEqArea(m, 0.9); // This is to shrink the equilibrium area by a little bit to avoid zero forces.
    
    double mSTriangleU = mSTriangle.computeEnergy(0.0);
    double mSVoronoiU = mSVoronoi.computeEnergy(0.0);
    EXPECT_NEAR(mSTriangleU, mSVoronoiU, abs(mSTriangleU + mSVoronoiU) * 1e-9);

}

TEST_F(MembraneFFTest, Force) {
    MembraneStretching<MembraneStretchingHarmonic> mSTriangle;
    MembraneStretching<MembraneStretchingVoronoiHarmonic> mSVoronoi;
    MembraneBending<MembraneBendingHelfrich> mBVoronoi;
    vector<MembraneInteractions*> memInteraction = {&mSTriangle, &mSVoronoi, &mBVoronoi};

    assignRandomForceAuxP(m, radius/500);
    resizeEqArea(m, 0.9); // This is to shrink the equilibrium area by a little bit to avoid zero forces.
    moveAlongForceAuxP(m, 2.5); // Also, previous tests show that the right octahedron shape has the minimum
                                // bending energy, so we also need to distort the shape a little bit, to avoid
                                // zero force of bending energy.
    recordCoordinate(m);

    // Compare the results with force predictions
    // U(x+h) - U(x-h) = dotProduct(2h, dU/dx) = -dotProduct(2h, F)

    for(auto eachMemInteraction: memInteraction) {
        resetCoordinate(m);
        resetForce(m);
        m->updateGeometry(true);
        eachMemInteraction->computeForces();

        moveAlongForceAuxP(m, 1.0);
        m->updateGeometry(true);
        double U1 = eachMemInteraction->computeEnergy(0.0);

        resetCoordinate(m);
        moveAlongForceAuxP(m, -1.0);
        m->updateGeometry(true);
        double U2 = eachMemInteraction->computeEnergy(0.0);

        double exDiff = 0.0;
        for(Vertex* v: m->getVertexVector())
            exDiff -= 2 * dotProduct(v->forceAuxP, v->force);

        EXPECT_NEAR(U1 - U2, exDiff, abs(exDiff / 1000))
            << eachMemInteraction->getName() << " force not working properly.";
    
    }

}


#  endif //DO_THIS_FF_MEMBRANE_TEST
#endif //TESTING


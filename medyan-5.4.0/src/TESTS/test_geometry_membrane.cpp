
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

#  define DO_THIS_GEOMETRY_MEMBRANE_TEST
#  ifdef DO_THIS_GEOMETRY_MEMBRANE_TEST

#    include "gtest/gtest.h"

#    include <random>

#    include "common.h"
#    include "MathFunctions.h"
using namespace mathfunc;
#    include "Rand.h"

#    include "Controller/GController.h"
#    include "SubSystem.h"

#    include "Vertex.hpp"
#    include "Edge.hpp"
#    include "Triangle.hpp"
#    include "GVoronoiCell.h"
#    include "GEdge.h"
#    include "GTriangle.h"
#    include "Membrane.hpp"
#    include "MembraneHierarchy.hpp"

namespace {
    using VertexData = tuple<array<double, 3>, vector<size_t>>;
    using MembraneData = vector<VertexData>;

    MembraneData membraneDataOctahedron(const array<double, 3>& center, double radius) {

        return MembraneData{
            VertexData({center[0], center[1], center[2]+radius}, {1, 2, 3, 4}),
            VertexData({center[0]+radius, center[1], center[2]}, {0, 4, 5, 2}),
            VertexData({center[0], center[1]+radius, center[2]}, {0, 1, 5, 3}),
            VertexData({center[0]-radius, center[1], center[2]}, {0, 2, 5, 4}),
            VertexData({center[0], center[1]-radius, center[2]}, {0, 3, 5, 1}),
            VertexData({center[0], center[1], center[2]-radius}, {4, 3, 2, 1})
        };
    }

    class MembraneGeometryTest: public ::testing::Test {
    protected:
        double radius;
        SubSystem s;
        MembraneData memData;
        Membrane *m;

        MembraneGeometryTest(): radius(100), memData(membraneDataOctahedron({2*radius, 2*radius, 2*radius}, radius)) {
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
            m = new Membrane(&s, 0, memData);
        }
        ~MembraneGeometryTest() {
            SysParams::GParams.cylinderNumMon.resize(0);
            delete m;
        }

    };

    void recordCoordinate(Membrane *m) {
        for(Vertex* it: m->getVertexVector()) it->coordinateP = it->coordinate;
    }
    void resetCoordinate(Membrane *m) {
        for(Vertex* it: m->getVertexVector()) it->coordinate = it->coordinateP;
    }
    void assignRandomForce(Membrane* m, double sigma) {
        normal_distribution<> nd(0, sigma);

        for(Vertex* it: m->getVertexVector()) {
            for(double& eachForce: it->force) {
                eachForce = nd(Rand::engFixed);
            }
        }
    }
    void moveAlongForce(Membrane* m, double d) {
        for(Vertex* it: m->getVertexVector()) {
            for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
                it->coordinate[coordIdx] += it->force[coordIdx] * d;
            }
        }
    }

    class MembraneHierarchyFlattened {
    public:
        int id;
        vector<MembraneHierarchyFlattened> children;
        MembraneHierarchyFlattened(int newId): id(newId) {}
        MembraneHierarchyFlattened(int newId, const vector<MembraneHierarchyFlattened>& newHie):
            id(newId), children(newHie) {}
    };
    bool operator==(const MembraneHierarchyFlattened& o1, const MembraneHierarchyFlattened& o2) {
        return (o1.id == o2.id) && (o1.children == o2.children);
    }
    ostream& operator<<(ostream& os, const MembraneHierarchyFlattened& mhf) {
        if(mhf.id == -1) os << "root"; else os << mhf.id;

        size_t n = mhf.children.size();
        if(n) {
            os << '(';
            for(size_t idx = 0; idx < n; ++idx) {
                os << mhf.children[idx];
                if(idx < n - 1) os << ", ";
            }
            os << ')';
        }
        return os;
    }
    MembraneHierarchyFlattened FlattenMembraneHierarchy(const MembraneHierarchy< Membrane >& root) {
        MembraneHierarchyFlattened res(root.getMembrane()? root.getMembrane()->getId(): -1);

        for(auto& child: root.children()) {
            res.children.push_back(FlattenMembraneHierarchy(*static_cast<MembraneHierarchy< Membrane >*>(child.get())));
        }

        return res;
    }
}

TEST_F(MembraneGeometryTest, Topology) {

    // Check that the vertices, edges and triangles are correctly registered.
    EXPECT_EQ(m->getVertexVector().size(), 6);
    EXPECT_EQ(m->getEdgeVector().size(), 12);
    EXPECT_EQ(m->getTriangleVector().size(), 8);
    
}

TEST_F(MembraneGeometryTest, Geometry) {

    /**************************************************************************
        Check normal geometry
    **************************************************************************/
    m->updateGeometry(true);

    // Check edge length and pseudo normal (for edge 0)
    double exEdgeLen = radius * sqrt(2);
    for(Edge* it: m->getEdgeVector())
        EXPECT_DOUBLE_EQ(it->getGEdge()->getLength(), exEdgeLen);

    auto& edgePseudoUnitNormal0 = m->getEdgeVector()[0]->getGEdge()->getPseudoUnitNormal();
    array<double, 3> exEdgePseudoUnitNormal0 { 1/sqrt(2), 0, 1/sqrt(2) };
    for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx)
        EXPECT_DOUBLE_EQ(edgePseudoUnitNormal0[coordIdx], exEdgePseudoUnitNormal0[coordIdx])
            << "Edge pseudo unit normal doesn't match at coordinate index " << coordIdx;

    // Check triangle area, angle and unit normal (for triangle 0)
    double exTriangleArea = radius * radius * sqrt(3) / 2;
    double exTriangleAngle = M_PI / 3;
    for(Triangle* it: m->getTriangleVector()) {
        EXPECT_DOUBLE_EQ(it->getGTriangle()->getArea(), exTriangleArea);
        for(double eachTheta: it->getGTriangle()->getTheta())
            EXPECT_DOUBLE_EQ(eachTheta, exTriangleAngle);
    }

    auto& triangleUnitNormal0 = m->getTriangleVector()[0]->getGTriangle()->getUnitNormal();
    array<double, 3> exTriangleUnitNormal0 { 1/sqrt(3), 1/sqrt(3), 1/sqrt(3) };
    for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx)
        EXPECT_DOUBLE_EQ(triangleUnitNormal0[coordIdx], exTriangleUnitNormal0[coordIdx])
            << "Triangle unit normal doesn't match at coordinate index " << coordIdx;
    
    // Check Voronoi cell area, curvature and unit normal (for vertex 0)
    double exVCellArea = radius * radius * sqrt(3) * 2 / 3;
    double exVCellCurv = 1 / radius;
    for(Vertex* it: m->getVertexVector()) {
        EXPECT_DOUBLE_EQ(it->getGVoronoiCell()->getArea(), exVCellArea);
        EXPECT_DOUBLE_EQ(it->getGVoronoiCell()->getCurv(), exVCellCurv);
    }

    auto& vertexPseudoUnitNormal0 = m->getVertexVector()[0]->getGVoronoiCell()->getPseudoUnitNormal();
    array<double, 3> exVertexPseudoUnitNormal0 { 0, 0, 1 };
    for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx)
        EXPECT_DOUBLE_EQ(vertexPseudoUnitNormal0[coordIdx], exVertexPseudoUnitNormal0[coordIdx])
            << "Vertex pseudo unit normal doesn't match at coordinate index " << coordIdx;
    
    // Check total volume enclosed by the membrane
    double exVolume = radius * radius * radius * 4 / 3;
    EXPECT_DOUBLE_EQ(m->getGMembrane()->getVolume(), exVolume);

    /**************************************************************************
        Check stretched geometry
    **************************************************************************/
    Vertex* v0 = m->getVertexVector()[0];
    size_t numNeighborV0 = v0->getNeighborNum();
    ASSERT_EQ(numNeighborV0, 4) << "The 0th vertex does not have the right number of neighbors.";

    // Assign force to 0th vertex
    v0->force[2] = -radius;

    m->updateGeometry(false, 1.0); // Moves the 0th vertex to the center of the octahedron.

    // Check edge length (for edge around vertex 0) and pseudo normal (for edge 0)
    double exStretchedEdgeLen = radius;
    for(Edge* it: v0->getNeighborEdges())
        EXPECT_DOUBLE_EQ(it->getGEdge()->getStretchedLength(), exStretchedEdgeLen);
    
    auto& stretchedEdgePseudoUnitNormal0 = m->getEdgeVector()[0]->getGEdge()->getStretchedPseudoUnitNormal();
    array<double, 3> exStretchedEdgePseudoUnitNormal0 { 0, 0, 1 };
    for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx)
        EXPECT_DOUBLE_EQ(stretchedEdgePseudoUnitNormal0[coordIdx], exStretchedEdgePseudoUnitNormal0[coordIdx])
            << "Edge stretched pseudo unit normal doesn't match at coordinate index " << coordIdx;
    
    // Check triangle area (for triangle around v0), angle (for triangle around v0) and normal vector (for triangle 0)
    double exStretchedTriangleArea = radius * radius / 2;
    double exStretchedTriangleAngleIn = M_PI / 2;
    double exStretchedTriangleAngleOut = M_PI / 4;
    for(size_t nIdx = 0; nIdx < numNeighborV0; ++nIdx) {
        EXPECT_DOUBLE_EQ(v0->getNeighborTriangles()[nIdx]->getGTriangle()->getStretchedArea(), exStretchedTriangleArea);
        for(size_t angleIdx = 0; angleIdx < 3; ++angleIdx) {
            if((3 - v0->getTriangleHead()[nIdx]) % 3 == angleIdx)
                EXPECT_DOUBLE_EQ(
                    v0->getNeighborTriangles()[nIdx]->getGTriangle()->getStretchedTheta()[angleIdx],
                    exStretchedTriangleAngleIn
                );
            else
                EXPECT_DOUBLE_EQ(
                    v0->getNeighborTriangles()[nIdx]->getGTriangle()->getStretchedTheta()[angleIdx],
                    exStretchedTriangleAngleOut
                );
        }
    }

    auto& stretchedTriangleUnitNormal0 = m->getTriangleVector()[0]->getGTriangle()->getStretchedUnitNormal();
    array<double, 3> exStretchedTriangleUnitNormal0 { 0, 0, 1 };
    for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx)
        EXPECT_DOUBLE_EQ(stretchedTriangleUnitNormal0[coordIdx], exStretchedTriangleUnitNormal0[coordIdx])
            << "Triangle stretched unit normal doesn't match at coordinate index " << coordIdx;
    
    // Check Voronoi cell area (for vertex 0) and curvature (for vertex 0)
    double exStretchedVCellArea = radius * radius;
    double exStretchedVCellCurv = 0;
    EXPECT_DOUBLE_EQ(v0->getGVoronoiCell()->getStretchedArea(), exStretchedVCellArea);
    EXPECT_DOUBLE_EQ(v0->getGVoronoiCell()->getStretchedCurv(), exStretchedVCellCurv);

    auto& stretchedVertexPseudoUnitNormal0 = m->getVertexVector()[0]->getGVoronoiCell()->getStretchedPseudoUnitNormal();
    array<double, 3> exStretchedVertexPseudoUnitNormal0 { 0, 0, 1 };
    for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx)
        EXPECT_DOUBLE_EQ(stretchedVertexPseudoUnitNormal0[coordIdx], exStretchedVertexPseudoUnitNormal0[coordIdx])
            << "Vertex stretched pseudo unit normal doesn't match at coordinate index " << coordIdx;
    
    // Check the total volume enclosed by the membrane
    double exStretchedVolume = radius * radius * radius * 2 / 3;
    EXPECT_DOUBLE_EQ(m->getGMembrane()->getStretchedVolume(), exStretchedVolume);

}

TEST_F(MembraneGeometryTest, SignedDistance) {
    // Check the correctness of the signed distance field.

    // Prerequisites
    m->updateGeometry(true);

    // On triangle
    EXPECT_DOUBLE_EQ(m->signedDistance({2*radius, 2*radius, 2*radius}, true), -radius / sqrt(3));
    EXPECT_DOUBLE_EQ(m->signedDistance({2*radius, 2*radius, 2*radius}, false), -radius / sqrt(3));

    EXPECT_DOUBLE_EQ(m->signedDistance({3*radius, 3*radius, 3*radius}, true), radius * 2.0 / sqrt(3));
    EXPECT_DOUBLE_EQ(m->signedDistance({3*radius, 3*radius, 3*radius}, false), radius * 2.0 / sqrt(3));

    EXPECT_DOUBLE_EQ(m->signedDistance({0.0, 0.0, 0.0}, true), radius * 5.0 / sqrt(3));
    EXPECT_DOUBLE_EQ(m->signedDistance({0.0, 0.0, 0.0}, false), radius * 5.0 / sqrt(3));

    // On edge
    EXPECT_DOUBLE_EQ(m->signedDistance({2*radius, 3*radius, 3*radius}, true), radius / sqrt(2));
    EXPECT_DOUBLE_EQ(m->signedDistance({2*radius, 3*radius, 3*radius}, false), radius / sqrt(2));

    EXPECT_DOUBLE_EQ(m->signedDistance({0, 0, 2*radius}, true), radius * 3.0 / sqrt(2));
    EXPECT_DOUBLE_EQ(m->signedDistance({0, 0, 2*radius}, false), radius * 3.0 / sqrt(2));

    // On vertex
    EXPECT_DOUBLE_EQ(m->signedDistance({4*radius, 2*radius, 2*radius}, true), radius);
    EXPECT_DOUBLE_EQ(m->signedDistance({4*radius, 2*radius, 2*radius}, false), radius);

    EXPECT_DOUBLE_EQ(m->signedDistance({2*radius, 2*radius, 0}, true), radius);
    EXPECT_DOUBLE_EQ(m->signedDistance({2*radius, 2*radius, 0}, false), radius);
    
}

TEST_F(MembraneGeometryTest, Derivative) {
    m->updateGeometry(true);
    recordCoordinate(m);
    assignRandomForce(m, radius/1000); // Simple test shows that radius/500 induces a change not small enough

    size_t numEdges = m->getEdgeVector().size();
    size_t numTriangles = m->getTriangleVector().size();
    size_t numVertices = m->getVertexVector().size();

    // Then move every vertex a little bit
    moveAlongForce(m, 1.0);
    m->updateGeometry(false, 0.0); // use stretched calculation here to speed things up
    
    vector<double> edgeLength1(numEdges);
    for(size_t idx = 0; idx < numEdges; ++idx) {
        edgeLength1[idx] = m->getEdgeVector()[idx]->getGEdge()->getStretchedLength();
    }
    vector<double> triangleArea1(numTriangles);
    vector<array<double, 3>> triangleTheta1(numTriangles, {{}});
    vector<array<double, 3>> triangleCotTheta1(numTriangles, {{}});
    for(size_t idx = 0; idx < numTriangles; ++idx) {
        triangleArea1[idx] = m->getTriangleVector()[idx]->getGTriangle()->getStretchedArea();
        triangleTheta1[idx] = m->getTriangleVector()[idx]->getGTriangle()->getStretchedTheta();
        triangleCotTheta1[idx] = m->getTriangleVector()[idx]->getGTriangle()->getStretchedCotTheta();
    }
    vector<double> vCellArea1(numVertices);
    vector<double> vCellCurv1(numVertices);
    for(size_t idx = 0; idx < numVertices; ++idx) {
        vCellArea1[idx] = m->getVertexVector()[idx]->getGVoronoiCell()->getStretchedArea();
        vCellCurv1[idx] = m->getVertexVector()[idx]->getGVoronoiCell()->getStretchedCurv();
    }
    double volume1 = m->getGMembrane()->getStretchedVolume();

    // Now move every vertex in the opposite direction
    resetCoordinate(m);
    moveAlongForce(m, -1.0);
    m->updateGeometry(false, 0.0);

    vector<double> edgeLength2(numEdges);
    for(size_t idx = 0; idx < numEdges; ++idx) {
        edgeLength2[idx] = m->getEdgeVector()[idx]->getGEdge()->getStretchedLength();
    }
    vector<double> triangleArea2(numTriangles);
    vector<array<double, 3>> triangleTheta2(numTriangles, {{}});
    vector<array<double, 3>> triangleCotTheta2(numTriangles, {{}});
    for(size_t idx = 0; idx < numTriangles; ++idx) {
        triangleArea2[idx] = m->getTriangleVector()[idx]->getGTriangle()->getStretchedArea();
        triangleTheta2[idx] = m->getTriangleVector()[idx]->getGTriangle()->getStretchedTheta();
        triangleCotTheta2[idx] = m->getTriangleVector()[idx]->getGTriangle()->getStretchedCotTheta();
    }
    vector<double> vCellArea2(numVertices);
    vector<double> vCellCurv2(numVertices);
    for(size_t idx = 0; idx < numVertices; ++idx) {
        vCellArea2[idx] = m->getVertexVector()[idx]->getGVoronoiCell()->getStretchedArea();
        vCellCurv2[idx] = m->getVertexVector()[idx]->getGVoronoiCell()->getStretchedCurv();
    }
    double volume2 = m->getGMembrane()->getStretchedVolume();

    // Compare the results with derivative predictions
    // A(x+h) - A(x-h) = dotProduct(2h, dA/dx)
    for(size_t idx = 0; idx < numEdges; ++idx) {
        Edge* e = m->getEdgeVector()[idx];
        // Edge length
        double exDiff = 0.0;
        for(size_t vIdx = 0; vIdx < 2; ++vIdx) {
            exDiff += 2 * dotProduct(
                e->getVertices()[vIdx]->force,
                array2Vector<double, 3>(e->getGEdge()->getDLength()[vIdx])
            );
        }
        EXPECT_NEAR(edgeLength1[idx] - edgeLength2[idx], exDiff, abs(exDiff / 1000));
    }
	for(size_t idx = 0; idx < numTriangles; ++idx) {
        Triangle *t = m->getTriangleVector()[idx];
        // Triangle area
        double exDiff = 0.0;
        for(size_t vIdx = 0; vIdx < 3; ++vIdx) {
            exDiff += 2 * dotProduct(
                t->getVertices()[vIdx]->force,
                array2Vector<double, 3>(t->getGTriangle()->getDArea()[vIdx])
            );
        }
        EXPECT_NEAR(triangleArea1[idx] - triangleArea2[idx], exDiff, abs(exDiff / 1000));
        // Triangle angles
        for(size_t aIdx = 0; aIdx < 3; ++aIdx) {
            exDiff = 0.0;
            for(size_t vIdx = 0; vIdx < 3; ++vIdx) {
                exDiff += 2 * dotProduct(
                    t->getVertices()[vIdx]->force,
                    array2Vector<double, 3>(t->getGTriangle()->getDTheta()[aIdx][vIdx])
                );
            }
            EXPECT_NEAR(triangleTheta1[idx][aIdx] - triangleTheta2[idx][aIdx], exDiff, abs(exDiff / 1000));
        }
        // Triangle angles (cot)
        for(size_t aIdx = 0; aIdx < 3; ++aIdx) {
            exDiff = 0.0;
            for(size_t vIdx = 0; vIdx < 3; ++vIdx) {
                exDiff += 2 * dotProduct(
                    t->getVertices()[vIdx]->force,
                    array2Vector<double, 3>(t->getGTriangle()->getDCotTheta()[aIdx][vIdx])
                );
            }
            EXPECT_NEAR(triangleCotTheta1[idx][aIdx] - triangleCotTheta2[idx][aIdx], exDiff, abs(exDiff / 1000));
        }
    }
    for(size_t idx = 0; idx < numVertices; ++idx) {
        Vertex *v = m->getVertexVector()[idx];
        size_t numNeighbor = v->getNeighborNum();
        // Voronoi cell area
        double exDiff = 0.0;
        exDiff += 2 * dotProduct(
            v->force,
            array2Vector<double, 3>(v->getGVoronoiCell()->getDArea())
        );
        for(size_t vIdx = 0; vIdx < numNeighbor; ++vIdx) {
            exDiff += 2 * dotProduct(
                v->getNeighborVertices()[vIdx]->force,
                array2Vector<double, 3>(v->getGVoronoiCell()->getDNeighborArea()[vIdx])
            );
        }
        EXPECT_NEAR(vCellArea1[idx] - vCellArea2[idx], exDiff, abs(exDiff / 1000));
        // Voronoi cell curvature
        exDiff = 0.0;
        exDiff += 2 * dotProduct(
            v->force,
            array2Vector<double, 3>(v->getGVoronoiCell()->getDCurv())
        );
        for(size_t vIdx = 0; vIdx < numNeighbor; ++vIdx) {
            exDiff += 2 * dotProduct(
                v->getNeighborVertices()[vIdx]->force,
                array2Vector<double, 3>(v->getGVoronoiCell()->getDNeighborCurv()[vIdx])
            );
        }
        EXPECT_NEAR(vCellCurv1[idx] - vCellCurv2[idx], exDiff, abs(exDiff / 1000));
    }
    {
        // Total volume
        double exDiff = 0.0;
        for(size_t vIdx = 0; vIdx < numVertices; ++vIdx) {
            exDiff += 2 * dotProduct(
                m->getVertexVector()[vIdx]->force,
                array2Vector(m->getGMembrane()->getDVolume()[vIdx])
            );
        }
        EXPECT_NEAR(volume1 - volume2, exDiff, abs(exDiff / 1000));
    }

}

TEST_F(MembraneGeometryTest, MembraneHierarchy< Membrane >) {
    vector<Membrane*> ms(6);
    ms[0] = new Membrane(&s, 0, membraneDataOctahedron({2000, 2000, 2000}, 500));
    ms[1] = new Membrane(&s, 0, membraneDataOctahedron({1900, 2000, 2000}, 5));
    ms[2] = new Membrane(&s, 0, membraneDataOctahedron({2100, 2000, 2000}, 50));
    ms[3] = new Membrane(&s, 0, membraneDataOctahedron({1900, 2000, 2000}, 50));
    ms[4] = new Membrane(&s, 0, membraneDataOctahedron({4000, 2000, 2000}, 50));
    ms[5] = new Membrane(&s, 0, membraneDataOctahedron({4000, 2000, 2000}, 500));

    int idShift = ms[0]->getId();

    for(auto eachM: ms) eachM->updateGeometry(true);

    MembraneHierarchy< Membrane > root(nullptr);
    for(auto eachM: ms) {
        MembraneHierarchy< Membrane >::addMembrane(eachM, root);
    }

    // Check printSelf function manually
    cout << endl << "Testing MembraneHierarchy output function, it should look nice and correct." << endl
        << "Expected containing relationship: "
        << FlattenMembraneHierarchy(root)
        << endl;
    root.printSelf();

    // Check hierarchy
    {
        MembraneHierarchyFlattened exHie { -1, {
            { 0 + idShift, {
                { 2 + idShift },
                { 3 + idShift, {
                    { 1 + idShift }
                }}
            }},
            { 5 + idShift, {
                { 4 + idShift }
            }}
        }};
        EXPECT_EQ(FlattenMembraneHierarchy(root), exHie)
            << "Membrane hierarchy incorrect.";
    }

    MembraneHierarchy< Membrane >::removeMembrane(ms[0], root);
    {
        MembraneHierarchyFlattened exHie { -1, {
            { 5 + idShift, {
                { 4 + idShift }
            }},
            { 2 + idShift },
            { 3 + idShift, {
                { 1 + idShift }
            }}
        }};
        EXPECT_EQ(FlattenMembraneHierarchy(root), exHie)
            << "Membrane hierarchy incorrect after removing membrane with id " << 0 + idShift;
    }

    MembraneHierarchy< Membrane >::removeMembrane(ms[3], root);
    {
        MembraneHierarchyFlattened exHie { -1, {
            { 5 + idShift, {
                { 4 + idShift }
            }},
            { 2 + idShift },
            { 1 + idShift }
        }};
        EXPECT_EQ(FlattenMembraneHierarchy(root), exHie)
            << "Membrane hierarchy incorrect after removing membrane with id " << 3 + idShift;
    }

    MembraneHierarchy< Membrane >::addMembrane(ms[0], root);
    {
        MembraneHierarchyFlattened exHie { -1, {
            { 5 + idShift, {
                { 4 + idShift }
            }},
            { 0 + idShift, {
                { 2 + idShift },
                { 1 + idShift }
            }}
        }};
        EXPECT_EQ(FlattenMembraneHierarchy(root), exHie)
            << "Membrane hierarchy incorrect after adding back membrane with id " << 0 + idShift;
    }


    for(auto eachM: ms) delete eachM;
}

#  endif //DO_THIS_GEOMETRY_MEMBRANE_TEST
#endif //TESTING


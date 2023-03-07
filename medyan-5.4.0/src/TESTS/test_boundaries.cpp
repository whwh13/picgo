
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

//#define DO_THIS_BOUNDARY_TEST
#ifdef DO_THIS_BOUNDARY_TEST

#include "gtest/gtest.h"

#include "common.h"

#include "BoundaryImpl.h"
#include "BoundaryElementImpl.h"
#include "SysParams.h"

TEST(PlaneBoundaryElementTest, Distances) {
    
    SysParams::GParams.compartmentSizeX = 10.0;
    SysParams::GParams.compartmentSizeY = 10.0;
    SysParams::GParams.compartmentSizeZ = 10.0;
    
    SysParams::GParams.NX = 5;
    SysParams::GParams.NY = 5;
    SysParams::GParams.NZ = 5;
    
    SysParams::GParams.nDim = 3;
    
    GController g;
    g.initializeGrid();
    
    BoundaryElement* b = new PlaneBoundaryElement({10.0,10.0,10.0}, {1,0,0}, 1.0, 1.0);
    
    ///test distance calculations
    EXPECT_EQ(5.0 ,b->distance({15,10,10}));
    EXPECT_EQ(-5.0 ,b->distance({5,10,10}));
    
    EXPECT_EQ(15.0 ,b->stretchedDistance({15,15,10}, {1.0,0,0}, 10.0));
    EXPECT_EQ(10.0 ,b->stretchedDistance({10,5,10}, {1.0,0,0}, 10.0));
}

TEST(SphereBoundaryElementTest, Distances) {
    
    SysParams::GParams.compartmentSizeX = 10.0;
    SysParams::GParams.compartmentSizeY = 10.0;
    SysParams::GParams.compartmentSizeZ = 10.0;
    
    SysParams::GParams.NX = 5;
    SysParams::GParams.NY = 5;
    SysParams::GParams.NZ = 5;
    
    SysParams::GParams.nDim = 3;
    
    GController g;
    g.initializeGrid();
    
    BoundaryElement* b = new SphereBoundaryElement({25.0,25.0,25.0}, 10.0, 1.0, 1.0);
    
    ///test distance calculations
    EXPECT_EQ(10.0, b->distance({25.0,25.0,25.0}));
    EXPECT_EQ(5.0, b->distance({25.0,30.0,25.0}));
    EXPECT_EQ(-5.0, b->distance({25.0, 40.0, 25.0}));
    
    EXPECT_EQ(0 ,b->stretchedDistance({25,25,25}, {1.0,0,0}, 10.0));
    EXPECT_EQ(-10.0 ,b->stretchedDistance({25,35,25}, {0,1,0}, 10.0));

}

TEST(CylindricalZBoundaryElementTest, Distances){
    
    SysParams::GParams.compartmentSizeX = 10.0;
    SysParams::GParams.compartmentSizeY = 10.0;
    SysParams::GParams.compartmentSizeZ = 10.0;
    
    SysParams::GParams.NX = 5;
    SysParams::GParams.NY = 5;
    SysParams::GParams.NZ = 5;
    
    SysParams::GParams.nDim = 3;
    
    GController g;
    g.initializeGrid();
    
    BoundaryElement* b = new CylindricalZBoundaryElement({25.0,25.0,25.0}, 10.0, 20.0, 1.0, 1.0);
    
    ///test distance calculations
    EXPECT_EQ(10.0, b->distance({25.0, 25.0, 25.0}));
    EXPECT_EQ(10.0, b->distance({25.0, 25.0, 30.0}));
    EXPECT_EQ(5.0, b->distance({25.0, 30.0, 30.0}));
    
    EXPECT_EQ(-5.0, b->distance({25.0, 40.0, 30.0}));
    EXPECT_EQ(-7.0, b->distance({25.0, 42.0, 16.0}));
    
    EXPECT_EQ(numeric_limits<double>::infinity(), b->distance({25.0, 40.0, 45.0}));
    EXPECT_EQ(numeric_limits<double>::infinity(), b->distance({25.0, 40.0, 5.0}));
    
    EXPECT_EQ(0.0 ,b->stretchedDistance({25,25,25}, {1.0,0,0}, 10.0));
    EXPECT_EQ(-10.0 ,b->stretchedDistance({25,35,25}, {0,1,0}, 10.0)); 
}

TEST(HalfSphereZBoundaryElementTest, Distances){
    
    SysParams::GParams.compartmentSizeX = 10.0;
    SysParams::GParams.compartmentSizeY = 10.0;
    SysParams::GParams.compartmentSizeZ = 10.0;
    
    SysParams::GParams.NX = 5;
    SysParams::GParams.NY = 5;
    SysParams::GParams.NZ = 5;
    
    SysParams::GParams.nDim = 3;
    
    GController g;
    g.initializeGrid();
    
    BoundaryElement* b = new HalfSphereZBoundaryElement({25.0,25.0,25.0}, 10.0, false, 1.0, 1.0);
    
    ///test distance calculations
    EXPECT_EQ(10.0, b->distance({25.0, 25.0, 25.0}));
    EXPECT_EQ(5.0, b->distance({25.0, 25.0, 30.0}));
    EXPECT_EQ(5.0, b->distance({25.0, 30.0, 25.0}));
    
    EXPECT_EQ(-5.0, b->distance({25.0, 40.0, 25.0}));
    EXPECT_EQ(-7.0, b->distance({25.0, 42.0, 25.0}));
    
    EXPECT_EQ(-10.0, b->distance({25.0, 25.0, 45.0}));
    EXPECT_EQ(numeric_limits<double>::infinity(), b->distance({25.0, 40.0, 10.0}));
    
    EXPECT_EQ(0.0 ,b->stretchedDistance({25.0,25.0,25.0}, {1.0,0,0}, 10.0));
}


#endif //DO_THIS_BOUNDARY_TEST
#endif //TESTING


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

#  ifndef MEDYAN_TEST_PUBLIC_H
#  define MEDYAN_TEST_PUBLIC_H

/* Some public classes and functions are defined here for test use.
*/

#    include "common.h"

#    include "CController.h"
#    include "Composite.h"
#    include "Controller/GController.h"
#    include "SubSystem.h"

namespace test_public {

    // Data type definition
    using VertexData    = tuple<array<double, 3>, vector<size_t>>;
    using MembraneData  = vector<VertexData>;

    // Dummy solid class of Composite
    class CompositeDummy: public Composite {
    public:
        short type;

        CompositeDummy(short newType): type(newType) {}

        virtual void printSelf()const override {
            cout << "This is a dummy composite object with type "
                 << type << endl;
        }
    
        virtual int getType() override { return type; }

    };

    // Initialize global variables to get a playground for test.
    // There is no need to undo this function to clean the test.
    inline void quickSetupPlayground(SubSystem* s, double size=1e10, size_t nx=1, size_t ny=1, size_t nz=1) {
        SysParams::GParams.compartmentSizeX = size;
        SysParams::GParams.compartmentSizeY = size;
        SysParams::GParams.compartmentSizeZ = size;
        
        SysParams::GParams.NX = nx;
        SysParams::GParams.NY = ny;
        SysParams::GParams.NZ = nz;
        
        SysParams::GParams.nDim = 3;

        GController g(s); // Dummy variable to initialize the compartments
        g.initializeGrid();
    }

    inline void quickSetupChem(SubSystem* s, string chemAlgorithm="GILLESPIE") {
        CController c(s);
        ChemistryData cData; // Dummy data
        c.initialize(chemAlgorithm, cData);
    }

    // Test object and data generators
    inline MembraneData membraneDataOctahedron(const array<double, 3>& center, double radius) {

        return MembraneData{
            VertexData({center[0], center[1], center[2]+radius}, {1, 2, 3, 4}),
            VertexData({center[0]+radius, center[1], center[2]}, {0, 4, 5, 2}),
            VertexData({center[0], center[1]+radius, center[2]}, {0, 1, 5, 3}),
            VertexData({center[0]-radius, center[1], center[2]}, {0, 2, 5, 4}),
            VertexData({center[0], center[1]-radius, center[2]}, {0, 3, 5, 1}),
            VertexData({center[0], center[1], center[2]-radius}, {4, 3, 2, 1})
        };
    }

}

#  endif // Include guard
#endif // TESTING

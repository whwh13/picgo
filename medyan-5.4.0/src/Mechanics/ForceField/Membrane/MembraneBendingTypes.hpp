#ifndef MEDYAN_Mechanics_ForceField_Membrane_MembraneBendingTypes_hpp
#define MEDYAN_Mechanics_ForceField_Membrane_MembraneBendingTypes_hpp

#include <vector>

namespace medyan {

// This structure specifies the curvature mismatch data corresponding to each protein on each vertex.
struct VertexCurvatureMismatchParams {
    // Parameters set before minimization.
    //----------------------------------
    // Number of molecules. Using floats in case we want to use averaged number of molecules.
    double mol = 0;

    // Parameters set during energy computation as a side effect.
    double energy = 0;
};

struct AllVertexCurvatureMismatchParams {
    int numProteins = 0;
    int numVertices = 0;
    std::vector< VertexCurvatureMismatchParams > value;

    void resize( int numProteins, int numVertices ) {
        this->numProteins = numProteins;
        this->numVertices = numVertices;
        value.resize( numProteins * numVertices );
    }

    auto& operator()( int proteinIndex, int vertexIndex ) {
        return value[ vertexIndex * numProteins + proteinIndex ];
    }
    auto& operator()( int proteinIndex, int vertexIndex ) const {
        return value[ vertexIndex * numProteins + proteinIndex ];
    }
};

} // namespace medyan

#endif

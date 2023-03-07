#ifndef MEDYAN_MedyanMeta_hpp
#define MEDYAN_MedyanMeta_hpp

#include <iostream>

#ifndef GIT_COMMIT_HASH
    #define GIT_COMMIT_HASH "?"
#endif
#ifndef GIT_BRANCH
    #define GIT_BRANCH "?"
#endif
#ifndef GIT_LATEST_TAG
    #define GIT_LATEST_TAG "?"
#endif


namespace medyan {

inline void printHeader() {
    using namespace std;

    // Copyright and version.
    //--------------------------------------------------------------------------
    cout << endl;
    cout << "*********************** MEDYAN ************************" << endl;
    cout << "   Simulation package for the Mechanochemical Dynamics " << endl;
    cout << "         of Active Networks, Third Generation.         " << endl;
    cout << "         PAPOIAN LAB 2015, ALL RIGHTS RESERVED         " << endl;
    cout << "*******************************************************" << endl;
    cout << "Commit hash                          " << GIT_COMMIT_HASH << endl;
    cout << "Git branch                           " << GIT_BRANCH      << endl;
    cout << "Latest version:                      " << GIT_LATEST_TAG  << endl;

    // Compile-time options.
    //--------------------------------------------------------------------------
    cout << "Memory model:                        "<< static_cast<unsigned>(8 * sizeof(void*))<<" bit"<<endl;

    cout << "Coordinate/Force precision:          ";
    #ifdef FLOAT_PRECISION
        cout << "single" << endl;
    #else
        cout << "double" << endl;
    #endif

    cout << "Pair-wise list algorithm:            ";
    #if defined NLORIGINAL
        cout << "original (legacy)"<<endl;
    #elif defined NLSTENCILLIST
        cout << "stencil list based"<<endl;
    #elif defined HYBRID_NLSTENCILLIST
        cout << "hybrid**"<<endl;
        cout << "**single neighbor list for compatible distance ranges corresponding to the same filament type pairs"<<endl;
    #elif defined SIMDBINDINGSEARCH
        cout << "SIMD**"<<endl;
        cout << "SIMD instruction set:                ";
        #if defined __AVX512F__
            cout<<"AVX512 (MEDYAN will use AVX2 instructions)"<<endl;
        #elif defined __AVX2__
            cout<<"AVX2"<<endl;
        #elif defined __AVX__
            cout<<"AVX"<<endl;
        #else
            cout<<"none"<<endl;
        #endif
        cout<<"**SIMD accelerated hybrid search - single neighbor list for compatible distance ranges corresponding to the same filament type pairs"<<endl;
    #endif
    cout << "*******************************************************" << endl;


    // Important messages to user.
    //--------------------------------------------------------------------------
    cout << "Important messages:" << endl;
    cout << "⚠️ Starting MEDYAN v5.1.0, diffusion coefficient of diffusing species in chemistry input uses the 3D diffusion coefficient with dimension nm²/s, replacing the previous reaction rate (with dimension 1/s) computed from diffusion coefficient and the size of a cubic compartment. For example, if 500 nm compartments were used and the diffusion reaction rate was 1.0, it should be changed to 0.25E6." << endl;
    cout << "⚠️ Starting MEDYAN v5.1.0, brancher output is similar to other linkers, where TWO 3D coordinates (on mother/daughter filament, respectively) are recorded. Previously, only the starting coordinates (on the mother filament) were recorded." << endl;
    cout << "*******************************************************" << endl;
}

} // namespace medyan

#endif

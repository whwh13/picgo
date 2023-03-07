#ifndef MEDYAN_Parameter_Bubble_hpp
#define MEDYAN_Parameter_Bubble_hpp

#include <filesystem>

#include "common.h"
#include "Parameter/ReactionEmissionAbsorption.hpp"
#include "Util/Math/Vec.hpp"
#include "Util/Parser/StringParser.hpp"

namespace medyan {

/// Struct to hold Bubble setup information
struct BubbleSetup {
    
    std::filesystem::path inputFile;
    
    ///If want a random distribution, used if inputFile is left blank
    int numBubbles = 0;
    ///Bubble type to create
    short bubbleType = 0;
};

struct ReactionEmissionAbsorptionSetupInit {
    ReactionEmissionAbsorptionSetup setup;
    int initNumSpecies1 = 0;
};

// MTOC related.
struct MTOCInit {
    // Underlying bubble information.
    //----------------------------------
    int                   bubbleType = 0;
    Vec<3, floatingpoint> bubbleCoord {};
    bool                  bubbleFixed = false;

    // Attached filament parameters.
    //----------------------------------
    int                   filamentType = 0;
    int                   numFilaments = 0;
    int                   numCylindersPerFilament = 0;
    // These parameters specify the initial angle spread of filaments.
    // Read MTOC attached filament creation code for more information.
    floatingpoint         theta1 = 0;
    floatingpoint         theta2 = 1;
    floatingpoint         phi1 = 0;
    floatingpoint         phi2 = 1;

    // mechanical parameters.
    //---------------------------------
    // Filament attachment parameters.
    floatingpoint         attachmentStretchingK = 0;

    // Chemical parameters.
    //---------------------------------
    // Emission-absorption reactions.
    std::vector<ReactionEmissionAbsorptionSetupInit> vecRxnEmiAbs;
    // Dynamic emission rate parameters. Only works when dyRateType is set to "force".
    // Assumes a piecewise-linear activation function.
    floatingpoint         emiForce1 = 0; // Minimum force required to activate emission.
    floatingpoint         emiForce2 = 0; // Minimum force that maximizes emission rate.
};
struct MTOCSettings {
    std::vector< MTOCInit >  initVec;
};

// AFM related.
struct AFMInit {
    // Underlying bubble information.
    //----------------------------------
    int                   bubbleType = 0;
    Vec<3, floatingpoint> bubbleCoord {};
    bool                  bubbleFixed = false;

    // Attached filament parameters.
    //----------------------------------
    int                   filamentType = 0;
    int                   numFilaments = 0;
    int                   numCylindersPerFilament = 0;
    // These parameters specify the initial angle spread of filaments.
    // Read AFM attached filament creation code for more information.
    floatingpoint         theta1 = 0;
    floatingpoint         theta2 = 1;
    floatingpoint         phi1 = 0;
    floatingpoint         phi2 = 1;

    // mechanical parameters.
    //---------------------------------
    // Filament attachment parameters.
    floatingpoint         attachmentStretchingK = 0;

    // Chemical parameters.
    //---------------------------------
    // Emission-absorption reactions.
    std::vector<ReactionEmissionAbsorptionSetupInit> vecRxnEmiAbs;
    // Dynamic emission rate parameters. Only works when dyRateType is set to "force".
    // Assumes a piecewise-linear activation function.
    floatingpoint         emiForce1 = 0; // Minimum force required to activate emission.
    floatingpoint         emiForce2 = 0; // Minimum force that maximizes emission rate.
};
struct AFMSettings {
    std::vector< AFMInit >  initVec;
};

} // namespace medyan

#endif

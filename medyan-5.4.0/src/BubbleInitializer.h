
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

#ifndef MEDYAN_BubbleInitializer_h
#define MEDYAN_BubbleInitializer_h

#include "MathFunctions.h"
#include "SysParams.h"
#include "Structure/SubSystem.h"

namespace medyan {


/// An interface to initialize an initial configuration of [Bubbles](@ref Bubble)
/// in the SubSystem.
/*!
 *  Bubblenitiazer class should be inherited to provide an intial scheme for
 *  filling a SubSystem with [Bubbles](@ref Bubble). The bubbles could be
 *  completely random or distributed in other ways.
 */
inline BubbleData createBubblesRandomDist(
    SubSystem&        sys,
    int               numBubbles,
    int               bubbleType,
    const MechParams& mechParams
) {
    BubbleData ret;

    auto& b = *sys.getBoundary();

    int bubbleCounter = 0;
    while(bubbleCounter < numBubbles) {
        auto coord = sys.getCompartmentGrid()->getRandomCoordinates();
        if(
            b.within(mathfunc::vec2Vector(coord)) &&
            b.distance(mathfunc::vec2Vector(coord)) > mechParams.BubbleRadius[bubbleType]
        ) {
            ret.bubbles.push_back({ bubbleType, coord });
            ++bubbleCounter;
        }
    }
    return ret;
}

} // namespace medyan

#endif

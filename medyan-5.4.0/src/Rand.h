
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

#ifndef MEDYAN_Rand_h
#define MEDYAN_Rand_h

#include <cmath>
#include <limits>
#include <random>
#include <type_traits>

#include "common.h"
#include "Util/Math/Vec.hpp"

namespace medyan {
/// A random number generator class.
class Rand {
    
public:
    // Defined with default seed. Must be seeded at initialization.
    inline static std::mt19937 eng {};
    #ifdef DEBUGCONSTANTSEED
        inline static int intcounter       = 0;
        inline static int floatcounter     = 0;
        inline static int chemistrycounter = 0;
	#endif

    // Get a random floatingpoint between low and high.
    static inline floatingpoint randfloatingpoint(floatingpoint low, floatingpoint high) {
        #ifdef DEBUGCONSTANTSEED
        floatcounter++;
		#endif
        return std::uniform_real_distribution<floatingpoint>(low, high)(eng);
    }
    // Get a random integer between low and high (inclusive).
    static inline int randInteger(int low, int high) {
        #ifdef DEBUGCONSTANTSEED
        intcounter++;
        #endif
        return std::uniform_int_distribution<int>(low, high)(eng);
    }

    static inline std::uint64_t randUInt64bit() {
        return std::uniform_int_distribution<std::uint64_t>{}(eng);
    }

    // Generate a random float with mean 0 and standard deviation 1.
    template< typename Float = floatingpoint, std::enable_if_t< std::is_floating_point_v<Float> >* = nullptr >
    static auto randn() {
        static std::normal_distribution<Float> nd;
        return nd(eng);
    }

    // Generate a random 3D unit vector.
    template< typename Float = floatingpoint, std::enable_if_t< std::is_floating_point_v<Float> >* = nullptr >
    static auto randUnitVector3() {
        std::uniform_real_distribution<Float> distTheta(0, 2 * M_PI);
        std::uniform_real_distribution<Float> distZ(-1, 1);
        const auto theta = distTheta(eng);
        const auto z = distZ(eng);
        const auto r = std::sqrt(1 - z * z);
        return medyan::Vec<3, Float> { r * std::cos(theta), r * std::sin(theta), z };
    }
};

namespace rand {

// Safe exponential distribution
//
// Generates exponential distribution safely:
//   - Precondition: λ is not negative (not checked)
//   - if λ is zero: return infinity.
//   - if λ is infinity: return zero.
//   - if λ is positive: set the λ for the specified exponential distribution,
//     and generate a random number which is not infinity.
//
// Inputs:
//   - d:      the exponential distribution object
//   - lambda: the exponential distribution parameter
//   - g:      the random number generator
//
// Output:
//   - The generated random number with type ExpDist::result_type
//
// Note:
//   - When λ = 0 or inf, the input distribution and generator are not touched.
//   - When λ > 0, the generator can be used 1 or more times.
//
// Background:
//   - std::exponential_distribution requires that λ > 0, though most
//     implementations would return ∞ when λ = 0.
//   - Some implementations of std::exponential_distribution returns infinity
//     when λ is positive. See
//     http://open-std.org/JTC1/SC22/WG21/docs/lwg-active.html#2524
//-----------------------------------------------------------------------------
template< typename ExpDist, typename FloatLambda, typename Generator >
inline auto safeExpDist(
    ExpDist&    d,
    FloatLambda lambda,
    Generator&  g
) {
    using ResType = typename ExpDist::result_type;

    ResType res;
    if(lambda == 0) {
        res = std::numeric_limits< ResType >::infinity();
    } else if(lambda == std::numeric_limits< FloatLambda >::infinity()) {
        res = 0;
    } else {
        // Set d's parameter
        d.param(typename ExpDist::param_type( lambda ));

        // Generate non-inf random number
        do res = d(g); while(std::isinf(res));
    }

    return res;
}

} // namespace rand
} // namespace medyan

#endif

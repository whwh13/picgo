#ifndef MEDYAN_TESTS_Mechanics_ForceField_TestFFCommon_Hpp
#define MEDYAN_TESTS_Mechanics_ForceField_TestFFCommon_Hpp

#include <algorithm>
#include <cmath>

#include "Rand.h"
#include "Util/Math/Vec.hpp"

// Provides common methods for force field testing purposes
namespace medyan::test_ff_common {

// Fill vectors with random numbers
template< typename ForwardIt, typename Float >
inline void fillNormalRand(ForwardIt first, ForwardIt last, Float mean, Float stddev) {
    std::normal_distribution< Float > nd(mean, stddev);
    std::generate(first, last, [&nd] { return nd(Rand::eng); });
}

template< size_t dim, typename Float, typename Container >
inline void fillNormalRand(medyan::VecArray< dim, Float, Container >& v, Float mean, Float stddev) {
    fillNormalRand(v.value.begin(), v.value.end(), mean, stddev);
}
template< typename Float >
inline void fillNormalRand(std::vector< Float >& v, Float mean, Float stddev) {
    fillNormalRand(v.begin(), v.end(), mean, stddev);
}

// Vector increment
template< size_t dim, typename Float, typename Container >
inline void vecInc(medyan::VecArray< dim, Float, Container >& v1, const medyan::VecArray< dim, Float, Container >& v2) {
    v1 += v2;
}
template< typename Float >
inline void vecInc(std::vector< Float >& v1, const std::vector< Float >& v2) {
    const int sz = std::min(v1.size(), v2.size());
    for(std::size_t i = 0; i < sz; ++i) {
        v1[i] += v2[i];
    }
}

// Vector decrement
template< size_t dim, typename Float, typename Container >
inline void vecDec(medyan::VecArray< dim, Float, Container >& v1, const medyan::VecArray< dim, Float, Container >& v2) {
    v1 -= v2;
}
template< typename Float >
inline void vecDec(std::vector< Float >& v1, const std::vector< Float >& v2) {
    const int sz = std::min(v1.size(), v2.size());
    for(std::size_t i = 0; i < sz; ++i) {
        v1[i] -= v2[i];
    }
}

// Dot product
template< size_t dim, typename Float, typename Container >
inline Float vecDot(const medyan::VecArray< dim, Float, Container >& v1, const medyan::VecArray< dim, Float, Container >& v2) {
    return medyan::dot(v1, v2);
}
template< typename Float >
inline Float vecDot(const std::vector< Float >& v1, const std::vector< Float >& v2) {
    Float res = 0;
    for(std::size_t i = 0; i < v1.size(); ++i) {
        res += v1[i] * v2[i];
    }
    return res;
}
// Dot product but with first n terms
template< typename Float >
inline Float vecDotFirstNTerms(const std::vector< Float >& v1, const std::vector< Float >& v2, int n) {
    Float res = 0;
    for(int i = 0; i < n; ++i) {
        res += v1[i] * v2[i];
    }
    return res;
}


// Compare equal
template< typename Float, std::enable_if_t< std::is_floating_point_v<Float> >* = nullptr >
inline bool equalRelEps(Float a, Float b, Float relEps) {
    return std::abs(a - b) <= std::min(std::abs(a), std::abs(b)) * relEps;
}
template< typename Float, std::enable_if_t< std::is_floating_point_v<Float> >* = nullptr >
inline bool equalAbsRelEps(Float a, Float b, Float absEps, Float relEps) {
    return std::abs(a - b) <= std::min(std::abs(a), std::abs(b)) * relEps + absEps;
}

// Energy - force consistency test
template< typename CoordContainer, typename Float >
struct EnergyForceConsistencyReport {
    bool passed;
    CoordContainer dc;
    Float deActual;
    Float deExpected;
};


// Note:
//   - calcEnergy and calcForce should handle the dependency computation and force propagation.
template<
    typename CoordContainer,
    typename Float,
    typename FuncEnergy,
    typename FuncForce
> inline auto testEnergyForceConsistency(
    CoordContainer c0, int numDof,
    FuncEnergy&& calcEnergy, FuncForce&& calcForce,
    Float moveMag, Float energyRelEps
) {
    EnergyForceConsistencyReport< CoordContainer, Float > res;

    res.dc.resize(numDof);
    fillNormalRand(res.dc, (Float)0.0, moveMag);
    CoordContainer c1 = c0; vecDec(c1, res.dc);
    CoordContainer c2 = c0; vecInc(c2, res.dc);

    CoordContainer f;
    f.resize(c0.size());

    // Actual change in energy
    const Float e1 = calcEnergy(c1);
    const Float e2 = calcEnergy(c2);
    res.deActual = e2 - e1;

    // Expected change in energy
    calcForce(c0, f);
    res.deExpected = -2 * vecDotFirstNTerms(res.dc, f, numDof);

    // Judge equality
    res.passed = equalRelEps(res.deActual, res.deExpected, energyRelEps);

    return res;
}

// This version is used if the input coordinate container contain all independent variables.
template<
    typename CoordContainer,
    typename Float,
    typename FuncEnergy,
    typename FuncForce
> inline auto testEnergyForceConsistency(
    const CoordContainer& c0,
    FuncEnergy&& calcEnergy, FuncForce&& calcForce,
    Float moveMag, Float energyRelEps
) {
    EnergyForceConsistencyReport< CoordContainer, Float > res;

    res.dc.resize(c0.size());
    fillNormalRand(res.dc, (Float)0.0, moveMag);
    CoordContainer c1 = c0; vecDec(c1, res.dc);
    CoordContainer c2 = c0; vecInc(c2, res.dc);

    CoordContainer f;
    f.resize(c0.size());

    // Actual change in energy
    const Float e1 = calcEnergy(c1);
    const Float e2 = calcEnergy(c2);
    res.deActual = e2 - e1;

    // Expected change in energy
    calcForce(c0, f);
    res.deExpected = -2 * vecDot(res.dc, f);

    // Judge equality
    res.passed = equalRelEps(res.deActual, res.deExpected, energyRelEps);

    return res;
}


template< typename CoordContainer, typename Float >
struct DependentCoordinateConsistencyReport {
    bool passed;
    CoordContainer dcIndependent;
    CoordContainer forceOriginal; // not propagated
    Float deActual;
    Float deExpected;
};


// Test the dependent variables calculation and the force propagation.
template< typename CoordContainer, typename Float, typename FuncCalcDep, typename FuncBackProp >
inline auto testDependentCoordinateConsistency(
    CoordContainer c0, int numDof,
    FuncCalcDep&& calcDep, FuncBackProp&& backProp,
    Float moveMag, Float forceMag, Float energyRelEps
) {
    DependentCoordinateConsistencyReport< CoordContainer, Float > res;

    res.dcIndependent.resize(numDof);
    fillNormalRand(res.dcIndependent, (Float)0.0, moveMag);
    calcDep(c0); // needed because it may be used by backprop.
    CoordContainer c1 = c0; vecDec(c1, res.dcIndependent); calcDep(c1);
    CoordContainer c2 = c0; vecInc(c2, res.dcIndependent); calcDep(c2);
    // Find the change in coordinates, including dependent coordinates.
    auto diff12Full = c2; vecDec(diff12Full, c1);

    // Assign random forces.
    res.forceOriginal.resize(c0.size());
    fillNormalRand(res.forceOriginal, (Float)0.0, forceMag);

    // Expected change in energy, calculated before propagating forces.
    res.deExpected = - vecDot(diff12Full, res.forceOriginal);

    // Actual change in energy, after force propagation
    auto forcePropagated = res.forceOriginal;
    backProp(c0, forcePropagated);
    res.deActual = -2 * vecDotFirstNTerms(res.dcIndependent, forcePropagated, numDof);

    // Judge equality
    res.passed = equalRelEps(res.deActual, res.deExpected, energyRelEps);

    return res;
}


template< typename CoordContainer, typename Float >
struct PushForwardConsistencyReport {
    bool passed;
    CoordContainer dcIndependent; // The tangent vector X chosen.
    CoordContainer dcFull;        // The pushed forward vector Y.
    int numVar = 0;       // Total number of variables including dependent variables.
    int numVarPassed = 0; // Total number of variables passing the check.
};

// Test independent tangent vector pushing forward.
//
// Let X be an arbitrary independent tangent vector at x, and Y be the pushed-forward vector in all variable space at y(x).
// Then we need to make sure that for a small |X|, we have
//   y(x + X) - y(x - X) â‰ˆ 2Y.
//
// Parameters:
// - c0:              initial value of x, but with dimension of y.
// - numDof:          dimension of independent variables (x).
// - calcDep:         function y(x).
// - pushForward:     function Dy(X).
// - moveMag:         magnitude to determine X.
// - absEps/relEps:   Ïµ = min{|a|, |b|} * relEps + absEps.
template<
    typename CoordContainer,
    typename Float,
    typename FuncCalcDep,
    typename FuncPushForward
> inline auto testPushForwardConsistency(
    CoordContainer c0, int numDof,
    FuncCalcDep&& calcDep, FuncPushForward&& pushForward,
    Float moveMag, Float absEps, Float relEps
) {
    PushForwardConsistencyReport< CoordContainer, Float > res {};
    res.numVar = c0.size();

    res.dcIndependent.resize(numDof);
    fillNormalRand(res.dcIndependent, (Float)0.0, moveMag);
    calcDep(c0); // needed because it may be used by pushforward.
    CoordContainer c1 = c0; vecDec(c1, res.dcIndependent); calcDep(c1);
    CoordContainer c2 = c0; vecInc(c2, res.dcIndependent); calcDep(c2);

    res.dcFull = res.dcIndependent;
    res.dcFull.resize(c0.size());
    pushForward(c0, res.dcFull);

    for(int i = 0; i < c0.size(); ++i) {
        const auto diffActual   = c2[i] - c1[i];
        const auto diffExpected = 2 * res.dcFull[i];
        if(equalAbsRelEps(diffActual, diffExpected, absEps, relEps)) {
            ++res.numVarPassed;
        }
    }

    // Passed?
    res.passed = res.numVar == res.numVarPassed;

    return res;
}


} // namespace medyan::test_ff_common

#endif

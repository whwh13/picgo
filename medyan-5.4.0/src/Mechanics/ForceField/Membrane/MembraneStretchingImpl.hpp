#ifndef MEDYAN_Mechanics_ForceField_Membrane_MembraneStretchingImpl_Hpp
#define MEDYAN_Mechanics_ForceField_Membrane_MembraneStretchingImpl_Hpp

#include "MathFunctions.h"

namespace medyan {

// A harmonic potential used by the MembraneStretching
struct MembraneStretchingHarmonic {

    double energy(double area, double kElastic, double eqArea) const {
        // kElastic is the elastic modulus, which is independent of the actual eqArea

        double dist = area - eqArea;

        return 0.5 * kElastic * dist * dist / eqArea;

    }
    
    void forces(
        FP* force, double area, const Vec3d& dArea, double kElastic, double eqArea
    ) const {
        // F_i = -grad_i U = -k / A_0 * (A - A_0) * grad_i A
        // A(rea) and grad_i A(rea) are obtained as function parameters

        const auto deltaF = (- kElastic * (area - eqArea) / eqArea) * dArea;

        for(int i = 0; i < 3; ++i) force[i] += deltaF[i]; // v->force += deltaF;
    }

};

// A linear potential used by the MembraneStretching
// Note: In energy minimization, other forces are needed to balance this force,
// since the preferable area would be -inf.
struct MembraneStretchingLinear {

    double energy(double area, double tension) const {
        return tension * area;
    }

    void forces(FP* force, const Vec3& dArea, double tension) const {
        const auto deltaF = - tension * dArea;
        for(size_t i = 0; i < 3; ++i) force[i] += deltaF[i];
    }
};

} // namespace medyan

#endif

#ifndef MEDYAN_Mechanics_ForceField_Membrane_MembraneBendingHelfrich_hpp
#define MEDYAN_Mechanics_ForceField_Membrane_MembraneBendingHelfrich_hpp

#include "MathFunctions.h"

namespace medyan {
/// A harmonic potential used by the MembraneStretching
struct MembraneBendingHelfrich {

    double energy(double area, double curv, double kBending, double eqCurv) const {
        const auto dist = curv - eqCurv;
        return 2 * kBending * dist * dist * area;
    }

    void forces(
        FP* force,
        double area, const medyan::Vec3& dArea,
        double curv, const medyan::Vec3& dCurv,
        double kBending, double eqCurv
    ) const {
        // F_i = -grad_i U = -4k (H - c0) A (grad_i H) - 2k (H - c0)^2 (grad_i A)
        // A(area), grad_i A(area), H(curv) and grad_i H(curv) are obtained as function parameters

        const auto dist = curv - eqCurv;
        const auto coeff1 = -4 * kBending * dist * area;
        const auto coeff2 = -2 * kBending * dist * dist;

        for(size_t i = 0; i < 3; ++i) force[i] += coeff1 * dCurv[i] + coeff2 * dArea[i];
    }

    // Compute forces, but treating area as a constant.
    void forcesConstArea(
        FP* force,
        double area,
        double curv, const Vec3d& dCurv,
        double kBending, double eqCurv,
        double scaleFactor
    ) const {
        // F_i = -grad_i U = -4k (H - c0) A (grad_i H)

        const auto dist = curv - eqCurv;
        for(int i = 0; i < 3; ++i) force[i] += -4 * kBending * dist * area * dCurv[i] * scaleFactor;
    }

};

struct MembraneBendingHelfrichQuadratic {
    double energy(double area, double curv2, double kBending) const {
        return 2 * kBending * curv2 * area;
    }

    void forces(
        FP* force,
        double area, const medyan::Vec3& dArea,
        double curv2, const medyan::Vec3& dCurv2,
        double kBending
    ) const {
        // F_i = -grad_i U = -2k A (grad_i H^2) - 2k H^2 (grad_i A)

        const auto coeff1 = -2 * kBending * area;
        const auto coeff2 = -2 * kBending * curv2;

        for(size_t i = 0; i < 3; ++i) force[i] += coeff1 * dCurv2[i] + coeff2 * dArea[i];
    }

};

} // namespace medyan

#endif

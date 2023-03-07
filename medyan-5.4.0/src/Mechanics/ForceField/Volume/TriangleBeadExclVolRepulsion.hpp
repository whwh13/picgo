
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v3.1
//
//  Copyright (2017-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_Mechanics_ForceField_Volume_TriangleBeadExclVolRepulsion_Hpp
#define MEDYAN_Mechanics_ForceField_Volume_TriangleBeadExclVolRepulsion_Hpp

#include "MathFunctions.h"

namespace medyan {

/// Represents a repulsive excluded volume potential used by the
/// TriangleBeadExclVolume template.
struct TriangleBeadExclVolRepulsion {
    
    double energy(Vec3, Vec3, Vec3, Vec3 cb, double area, double kExVol) const;
    
    void forces(
        floatingpoint* f0, floatingpoint* f1, floatingpoint* f2, floatingpoint* fb,
        Vec3, Vec3, Vec3, Vec3 cb,
        double area, const Vec3&, const Vec3&, const Vec3&,
        double kExVol
    ) const;

    Vec3 loadForces(
        const Vec3& c0, const Vec3& c1, const Vec3& c2, const Vec< 3, floatingpoint >& coord,
        double area, double kExVol
    ) const;
};

} // namespace medyan

#endif


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

#include "Mechanics/ForceField/Volume/TriangleBeadExclVolRepulsion.hpp"

#include "MathFunctions.h"


namespace medyan {
using namespace mathfunc;

/**
 * Implements the volume exclusion force field of triangle and bead.
 * 
 * Please refer to the document files for the math derivation of this force
 * field. Developers of this force field are responsible of keeping the
 * related documents up to date.
 */

// Note: These functions require that the area of the triangle has already been calculated

double TriangleBeadExclVolRepulsion::energy(
    Vec3 c0, Vec3 c1, Vec3 c2, Vec3 cb, double area,
    double kExVol
) const {

    //check if in same plane
    auto cp = cross(c1 - c0, c2 - c0);
    if(areEqual(dot(cp, cb - c0), 0.0)) {
        
        // slightly move point. Using negative number here to move "into the cell".
        cb -= cp * (0.01 / magnitude(cp));
    }

    const double A = dot(c1 - c0, c1 - c0);
    const double B = dot(c2 - c1, c2 - c1);
    const double C = dot(c0 - cb, c0 - cb);
    const double D = 2 * dot(c1 - c0, c0 - cb);
    const double E = 2 * dot(c2 - c1, c0 - cb);
    const double F = 2 * dot(c1 - c0, c2 - c1);

    const double A1 = 2 * A*E - D*F;
    const double A2 = 2 * B*D - 2 * A*E + (D - E)*F;
    const double A3 = -4 * A*B - 2 * B*D + F*(E + F);

    const double B1 = 4 * A*C - D*D;
    const double B2 = 4 * A*C - D*D + 4 * B*C - E*E + 4 * C*F - 2 * D*E;
    const double B3 = 4 * B*A - F*F + 4 * B*C - E*E + 4 * B*D - 2 * E*F;
    const double BB1 = sqrt(B1);
    const double BB2 = sqrt(B2);
    const double BB3 = sqrt(B3);

    const double C1 = 2 * A + D;
    const double C2 = 2 * A + D + E + 2 * (B + F);
    const double C3 = 2 * B + E + F;
    const double D1 = D;
    const double D2 = D + E;
    const double D3 = E + F;

    const double E1 = atan(C1 / BB1);
    const double E2 = atan(C2 / BB2);
    const double E3 = atan(C3 / BB3);
    const double F1 = atan(D1 / BB1);
    const double F2 = atan(D2 / BB2);
    const double F3 = atan(D3 / BB3);

    const double G1 = A1 / BB1;
    const double G2 = A2 / BB2;
    const double G3 = A3 / BB3;

    const double numerator = G1*(E1 - F1) + G2*(E2 - F2) + G3*(E3 - F3);
    const double denominator = B*D*D + A*(-4 * B*C + E*E) + F*(-D*E + C*F);

    const double H = numerator / denominator;
    
    const double energy = kExVol * 2 * area * H;
    
    return energy;
}

void TriangleBeadExclVolRepulsion::forces(
    floatingpoint* f0, floatingpoint* f1, floatingpoint* f2, floatingpoint* fb,
    Vec3 c0, Vec3 c1, Vec3 c2, Vec3 cb,
    double area, const Vec3& dArea0, const Vec3& dArea1, const Vec3& dArea2,
    double kExVol
) const {
    
    //check if in same plane
    auto cp = cross(c1 - c0, c2 - c0);
    if(areEqual(dot(cp, cb - c0), 0.0)) {
        
        // slightly move point. Using negative number here to move "into the cell".
        cb -= cp * (0.01 / magnitude(cp));
    }

    const double A = dot(c1 - c0, c1 - c0);
    const double B = dot(c2 - c1, c2 - c1);
    const double C = dot(c0 - cb, c0 - cb);
    const double D = 2 * dot(c1 - c0, c0 - cb);
    const double E = 2 * dot(c2 - c1, c0 - cb);
    const double F = 2 * dot(c1 - c0, c2 - c1);

    // Derivative index is 0, 1, 2, b
    const array<Vec3d, 4> dA = {{
        2 * (c0 - c1),
        2 * (c1 - c0)
    }};
    const array<Vec3d, 4> dB = {{
        { 0, 0, 0 },
        2 * (c1 - c2),
        2 * (c2 - c1)
    }};
    const array<Vec3d, 4> dC = {{
        2 * (c0 - cb),
        { 0, 0, 0 },
        { 0, 0, 0 },
        2 * (cb - c0)
    }};
    const array<Vec3d, 4> dD = {{
        2 * (c1 + cb - 2*c0),
        2 * (c0 - cb),
        { 0, 0, 0 },
        2 * (c0 - c1)
    }};
    const array<Vec3d, 4> dE = {{
        2 * (c2 - c1),
        2 * (cb - c0),
        2 * (c0 - cb),
        2 * (c1 - c2)
    }};
    const array<Vec3d, 4> dF = {{
        2 * (c1 - c2),
        2 * (c0 + c2 - 2*c1),
        2 * (c1 - c0)
    }};

    const double A1 = 2 * A*E - D*F;
    const double A2 = 2 * B*D - 2 * A*E + (D - E)*F;
    const double A3 = -4 * A*B - 2 * B*D + F*(E + F);

    const double B1 = 4 * A*C - D*D;
    const double B2 = 4 * A*C - D*D + 4 * B*C - E*E + 4 * C*F - 2 * D*E;
    const double B3 = 4 * B*A - F*F + 4 * B*C - E*E + 4 * B*D - 2 * E*F;
    const double BB1 = sqrt(B1);
    const double BB2 = sqrt(B2);
    const double BB3 = sqrt(B3);

    const double C1 = 2 * A + D;
    const double C2 = 2 * A + D + E + 2 * (B + F);
    const double C3 = 2 * B + E + F;
    const double D1 = D;
    const double D2 = D + E;
    const double D3 = E + F;

    const double E1 = atan(C1 / BB1);
    const double E2 = atan(C2 / BB2);
    const double E3 = atan(C3 / BB3);
    const double F1 = atan(D1 / BB1);
    const double F2 = atan(D2 / BB2);
    const double F3 = atan(D3 / BB3);

    const double G1 = A1 / BB1;
    const double G2 = A2 / BB2;
    const double G3 = A3 / BB3;

    const double numerator = G1*(E1 - F1) + G2*(E2 - F2) + G3*(E3 - F3);
    const double denominator = B*D*D + A*(-4 * B*C + E*E) + F*(-D*E + C*F);

    const double H = numerator / denominator;
    
    array<Vec3d, 4> dH = {};
    for(size_t dIdx = 0; dIdx < 4; ++dIdx) {
        for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
            const double dAT = dA[dIdx][coordIdx];
            const double dBT = dB[dIdx][coordIdx];
            const double dCT = dC[dIdx][coordIdx];
            const double dDT = dD[dIdx][coordIdx];
            const double dET = dE[dIdx][coordIdx];
            const double dFT = dF[dIdx][coordIdx];

            const double dA1 = 2 * (dAT*E + A*dET) - (dDT*F + D*dFT);
            const double dA2 = 2 * (dBT*D + B*dDT) - 2 * (dAT*E + A*dET) + ((dDT - dET)*F + (D - E)*dFT);
            const double dA3 = -4 * (dAT*B + A*dBT) - 2 * (dBT*D + B*dDT) + (dFT*(E + F) + F*(dET + dFT));

            const double dB1 = 4 * (dAT*C + A*dCT) - 2 * D*dDT;
            const double dB2 = 4 * (dAT*C + A*dCT) - 2 * D*dDT + 4 * (dBT*C + B*dCT) - 2 * E*dET + 4 * (dCT*F + C*dFT) - 2 * (dDT*E + D*dET);
            const double dB3 = 4 * (dBT*A + B*dAT) - 2 * F*dFT + 4 * (dBT*C + B*dCT) - 2 * E*dET + 4 * (dBT*D + B*dDT) - 2 * (dET*F + E*dFT);
            const double dBB1 = dB1 / 2 / BB1;
            const double dBB2 = dB2 / 2 / BB2;
            const double dBB3 = dB3 / 2 / BB3;

            const double dC1 = 2 * dAT + dDT;
            const double dC2 = 2 * dAT + dDT + dET + 2 * (dBT + dFT);
            const double dC3 = 2 * dBT + dET + dFT;
            const double dD1 = dDT;
            const double dD2 = dDT + dET;
            const double dD3 = dET + dFT;

            const double dE1 = (BB1*dC1 - C1*dBB1) / (B1 + C1*C1);
            const double dE2 = (BB2*dC2 - C2*dBB2) / (B2 + C2*C2);
            const double dE3 = (BB3*dC3 - C3*dBB3) / (B3 + C3*C3);
            const double dF1 = (BB1*dD1 - D1*dBB1) / (B1 + D1*D1);
            const double dF2 = (BB2*dD2 - D2*dBB2) / (B2 + D2*D2);
            const double dF3 = (BB3*dD3 - D3*dBB3) / (B3 + D3*D3);

            const double dG1 = (BB1*dA1 - A1*dBB1) / B1;
            const double dG2 = (BB2*dA2 - A2*dBB2) / B2;
            const double dG3 = (BB3*dA3 - A3*dBB3) / B3;

            const double dNumerator = dG1*(E1 - F1) + G1*(dE1 - dF1) + dG2*(E2 - F2) + G2*(dE2 - dF2) + dG3*(E3 - F3) + G3*(dE3 - dF3);
            const double dDenominator = dBT * D*D + B * 2 * D*dDT + dAT*(-4 * B*C + E*E) + A*(-4 * (dBT*C + B*dCT) + 2 * E*dET) + dFT*(-D*E + C*F) + F*(-(dDT*E + D*dET) + dCT*F + C*dFT);

            dH[dIdx][coordIdx] = (denominator*dNumerator - numerator*dDenominator) / (denominator*denominator);
        }
    }

    for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
        fb[coordIdx] -= kExVol * dH[3][coordIdx] * 2 * area;
        f0[coordIdx] -= kExVol * 2 * (dArea0[coordIdx] * H + dH[0][coordIdx] * area);
        f1[coordIdx] -= kExVol * 2 * (dArea1[coordIdx] * H + dH[1][coordIdx] * area);
        f2[coordIdx] -= kExVol * 2 * (dArea2[coordIdx] * H + dH[2][coordIdx] * area);
    }
}

Vec3 TriangleBeadExclVolRepulsion::loadForces(
    const Vec3& c0, const Vec3& c1, const Vec3& c2, const Vec< 3, floatingpoint >& coord,
    double area, double kExVol
) const {
    Vec3 cb (coord);

    //check if in same plane
    auto cp = cross(c1 - c0, c2 - c0);
    if(areEqual(dot(cp, cb - c0), 0.0)) {
        
        // slightly move point. Using negative number here to move "into the cell".
        cb -= cp * (0.01 / magnitude(cp));
    }
    
    const double A = dot(c1 - c0, c1 - c0);
    const double B = dot(c2 - c1, c2 - c1);
    const double C = dot(c0 - cb, c0 - cb);
    const double D = 2 * dot(c1 - c0, c0 - cb);
    const double E = 2 * dot(c2 - c1, c0 - cb);
    const double F = 2 * dot(c1 - c0, c2 - c1);

    // Only consider derivative on the coordinate of the bead (cb)
    // Vec3 dA = {{}};
    // Vec3 dB = {{}};
    const Vec3d dC = 2 * (cb - c0);
    const Vec3d dD = 2 * (c0 - c1);
    const Vec3d dE = 2 * (c1 - c2);
    const Vec3d dF = 2 * (c1 - c0);

    const double A1 = 2 * A*E - D*F;
    const double A2 = 2 * B*D - 2 * A*E + (D - E)*F;
    const double A3 = -4 * A*B - 2 * B*D + F*(E + F);

    const double B1 = 4 * A*C - D*D;
    const double B2 = 4 * A*C - D*D + 4 * B*C - E*E + 4 * C*F - 2 * D*E;
    const double B3 = 4 * B*A - F*F + 4 * B*C - E*E + 4 * B*D - 2 * E*F;
    const double BB1 = sqrt(B1);
    const double BB2 = sqrt(B2);
    const double BB3 = sqrt(B3);

    const double C1 = 2 * A + D;
    const double C2 = 2 * A + D + E + 2 * (B + F);
    const double C3 = 2 * B + E + F;
    const double D1 = D;
    const double D2 = D + E;
    const double D3 = E + F;

    const double E1 = atan(C1 / BB1);
    const double E2 = atan(C2 / BB2);
    const double E3 = atan(C3 / BB3);
    const double F1 = atan(D1 / BB1);
    const double F2 = atan(D2 / BB2);
    const double F3 = atan(D3 / BB3);

    const double G1 = A1 / BB1;
    const double G2 = A2 / BB2;
    const double G3 = A3 / BB3;

    const double numerator = G1*(E1 - F1) + G2*(E2 - F2) + G3*(E3 - F3);
    const double denominator = B*D*D + A*(-4 * B*C + E*E) + F*(-D*E + C*F);

    // Only consider derivative on the coordinate of the bead (cb)
    Vec3 dH {};
    for(size_t coordIdx = 0; coordIdx < 3; ++coordIdx) {
        const double dCT = dC[coordIdx];
        const double dDT = dD[coordIdx];
        const double dET = dE[coordIdx];
        const double dFT = dF[coordIdx];

        const double dA1 = 2 * (A*dET) - (dDT*F + D*dFT);
        const double dA2 = 2 * (B*dDT) - 2 * (A*dET) + ((dDT - dET)*F + (D - E)*dFT);
        const double dA3 = - 2 * (B*dDT) + (dFT*(E + F) + F*(dET + dFT));

        const double dB1 = 4 * (A*dCT) - 2 * D*dDT;
        const double dB2 = 4 * (A*dCT) - 2 * D*dDT + 4 * (B*dCT) - 2 * E*dET + 4 * (dCT*F + C*dFT) - 2 * (dDT*E + D*dET);
        const double dB3 = - 2 * F*dFT + 4 * (B*dCT) - 2 * E*dET + 4 * (B*dDT) - 2 * (dET*F + E*dFT);
        const double dBB1 = dB1 / 2 / BB1;
        const double dBB2 = dB2 / 2 / BB2;
        const double dBB3 = dB3 / 2 / BB3;

        const double dC1 = dDT;
        const double dC2 = dDT + dET + 2 * dFT;
        const double dC3 = dET + dFT;
        const double dD1 = dDT;
        const double dD2 = dDT + dET;
        const double dD3 = dET + dFT;

        const double dE1 = (BB1*dC1 - C1*dBB1) / (B1 + C1*C1);
        const double dE2 = (BB2*dC2 - C2*dBB2) / (B2 + C2*C2);
        const double dE3 = (BB3*dC3 - C3*dBB3) / (B3 + C3*C3);
        const double dF1 = (BB1*dD1 - D1*dBB1) / (B1 + D1*D1);
        const double dF2 = (BB2*dD2 - D2*dBB2) / (B2 + D2*D2);
        const double dF3 = (BB3*dD3 - D3*dBB3) / (B3 + D3*D3);

        const double dG1 = (BB1*dA1 - A1*dBB1) / B1;
        const double dG2 = (BB2*dA2 - A2*dBB2) / B2;
        const double dG3 = (BB3*dA3 - A3*dBB3) / B3;

        const double dNumerator = dG1*(E1 - F1) + G1*(dE1 - dF1) + dG2*(E2 - F2) + G2*(dE2 - dF2) + dG3*(E3 - F3) + G3*(dE3 - dF3);
        const double dDenominator = B * 2 * D*dDT + A*(-4 * (B*dCT) + 2 * E*dET) + dFT*(-D*E + C*F) + F*(-(dDT*E + D*dET) + dCT*F + C*dFT);

        dH[coordIdx] = (denominator*dNumerator - numerator*dDenominator) / (denominator*denominator);
    }

    const Vec3d forceBead = (-kExVol * 2 * area) * dH;

    return forceBead;

}

} // namespace medyan

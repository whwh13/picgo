#include "Mechanics/ForceField/Branching/BranchingDihedralQuadratic.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint> // uint_fast8_t
#include <limits>

#include "Mechanics/ForceField/Branching/BranchingDihedral.h"
#include "Structure/BranchingPoint.h"
#include "Util/Io/Log.hpp"
#include "Util/Math/Vec.hpp"
//This version uses a the following vectors to determine dihedral angles.
//b1 = c2 - mp;
//b2 = c3 - mp;
//b3 = c4 - c3;
//n1 = b1 x b2;
//n2 = b3 x b2;
namespace medyan {
using namespace mathfunc;
using V3 = Vec< 3, floatingpoint >;

constexpr floatingpoint cosDihTol = 0.01;
constexpr floatingpoint angSinMin = 0.001;
constexpr floatingpoint dihSinMin = 0.00001;

floatingpoint BranchingDihedralQuadratic::energy(
    const floatingpoint *coord, size_t nint,
    const unsigned int *beadSet, const floatingpoint *kdih, const floatingpoint *pos
) const {

    // Beads per interaction
    constexpr std::uint_fast8_t bpi = 4;
    static_assert(
        bpi == BranchingDihedral< BranchingDihedralQuadratic >::n,
        "Number of beads per interaction in branching dihedral quadratic does not match"
    );
    floatingpoint *coord2temp = new floatingpoint[3];
    floatingpoint U = 0.0;

    for(size_t i = 0; i < nint; ++i) {
        const auto coord1 = makeRefVec< 3 >(coord + beadSet[bpi * i    ]);
        auto p = pos[i];
        Vec<3,floatingpoint> coord2;
        coord2 = makeRefVec< 3 >(coord + beadSet[bpi * i + 1]);
        if(areEqual(p, 1.0)) {
            p = (floatingpoint) 0.5;
            //cplus = cminus +p*(cplus_extended-cminus)
            coord2 = (1 / p) *(coord2 - (1 - p) * coord1);
        }
        const auto coord3 = makeRefVec< 3 >(coord + beadSet[bpi * i + 2]);
        const auto coord4 = makeRefVec< 3 >(coord + beadSet[bpi * i + 3]);

        // Brancher coordinate on the mother filament
        const auto mp = (1 - p) * coord1 + p * coord2;

        // Bonds
        const auto b1 = coord2 - mp;
        const auto b2 = coord3 - mp;
        const auto b3 = coord4 - coord3;

        // Triangle normals
        const auto n1    = cross(b1, b2);
        const auto n2    = cross(b3, b2);
        const auto n1mag = magnitude(n1);
        const auto n2mag = magnitude(n2);

        // cos(theta)
        auto ct = dot(n1, n2) / (n1mag * n2mag);

        if(ct > 1.0 + cosDihTol || ct < -1.0 - cosDihTol) {
            LOG(WARNING) << "Dihedral incorrect cos theta: " << ct;
            LOG(INFO) << "Interaction information:\n"
                << "Bead coords: " << coord1 << ' ' << coord2 << ' ' << coord3 << ' ' << coord4 << '\n'
                << "Position: " << p << ", brancher coord: " << mp;
        }

        if(ct < -1.0) ct = -1.0;
        if(ct >  1.0) ct =  1.0;

        const auto theta = std::acos(ct);

        const auto U_i = kdih[i] * theta * theta;

        U += U_i;
    } // End loop interactions
//	cout<<"Quadratic   "<<U<<endl;
    delete [] coord2temp;
    return U;

} // floatingpoint energy(...)

void BranchingDihedralQuadratic::forces(
    const floatingpoint *coord, floatingpoint *f, size_t nint,
    const unsigned int *beadSet, const floatingpoint *kdih, const floatingpoint *pos,
    floatingpoint *stretchforce
) const {

    // Beads per interaction
    constexpr std::uint_fast8_t bpi = 4;
    static_assert(
        bpi == BranchingDihedral< BranchingDihedralQuadratic >::n,
        "Number of beads per interaction in branching dihedral quadratic does not match"
    );

    for(size_t i = 0; i < nint; ++i) {
        const auto coord1 = makeRefVec< 3 >(coord + beadSet[bpi * i    ]);
        auto p = pos[i];
        Vec<3,floatingpoint> coord2;
        coord2 = makeRefVec< 3 >(coord + beadSet[bpi * i + 1]);
        if(areEqual(p, 1.0)) {
            p = (floatingpoint) 0.5;
            coord2 = (1 / p) *(coord2 - (1 - p) * coord1);
        }
        const auto coord3 = makeRefVec< 3 >(coord + beadSet[bpi * i + 2]);
        const auto coord4 = makeRefVec< 3 >(coord + beadSet[bpi * i + 3]);

        auto f1 = makeRefVec< 3, floatingpoint >(f + beadSet[bpi * i    ]);
        auto f2 = makeRefVec< 3, floatingpoint >(f + beadSet[bpi * i + 1]);
        auto f3 = makeRefVec< 3, floatingpoint >(f + beadSet[bpi * i + 2]);
        auto f4 = makeRefVec< 3, floatingpoint >(f + beadSet[bpi * i + 3]);

        // Brancher coordinate on the mother filament
        const auto mp = (1 - p) * coord1 + p * coord2;

        // Bonds
        const auto b1 = coord2 - mp;
        const auto b2 = coord3 - mp;
        const auto b3 = coord4 - coord3;

        //---------------------------------------------------------------------
        // E = E1(b1, b2, b3)
        //   = E2(c2, mp, c3, c4)
        //   = E3(c1, c2, c3, c4)
        //
        // If we have dE1, and let E1i = dE1 / dbi, then
        //
        //   dE2 / dc2 = E11
        //   dE2 / dmp = - E11 - E12
        //   dE2 / dc3 = E12 - E13
        //   dE2 / dc4 = E13
        //
        //   dE3 / dc1 = - (1-p)E11 - (1-p)E12
        //   dE3 / dc2 = (1-p)E11 - p E12
        //   dE3 / dc3 = E12 - E13
        //   dE3 / dc4 = E13
        //
        //---------------------------------------------------------------------
        //                 (b1 x b2) . (b3 x b2)
        // cos(theta) = ----------------------------
        //                 |b1 x b2|   |b3 x b2|
        //
        //                 (b1 . b3) (b2 . b2) - (b1 . b2) (b3 . b2)
        //            = -----------------------------------------------
        //                 |b1| |b3| |b2|^2  sin<b1, b2> sin<b3, b2>
        //
        //                 cos<b1, b3> - cos<b1, b2> cos<b3, b2>
        //            = -------------------------------------------
        //                       sin<b1, b2>  sin<b3, b2>
        //
        //---------------------------------------------------------------------
        // Useful formula
        //
        // d cos<x1, x2> / dx1
        //
        //           x2         (x1 . x2)  x1
        //    = ----------- - -----------------
        //       |x1| |x2|       |x1|^3 |x2|
        //
        //                  [    x2          x1    ]
        //    = cos<x1, x2> [ --------- - -------- ]
        //                  [  x1 . x2     |x1|^2  ]
        //
        // d sin<x1, x2> / dx1
        //
        //         cos<x1, x2>
        //    = - ------------- d cos<x1, x2> / dx1
        //         sin<x1, x2>
        //
        //---------------------------------------------------------------------
        const auto b1mag2 = magnitude2(b1);
        const auto b2mag2 = magnitude2(b2);
        const auto b3mag2 = magnitude2(b3);
        const auto b1mag  = std::sqrt(b1mag2);
        const auto b2mag  = std::sqrt(b2mag2);
        const auto b3mag  = std::sqrt(b3mag2);

        const auto b1mag2inv = (floatingpoint)1.0 / b1mag2;
        const auto b2mag2inv = (floatingpoint)1.0 / b2mag2;
        const auto b3mag2inv = (floatingpoint)1.0 / b3mag2;
        const auto mag13inv  = (floatingpoint)1.0 / (b1mag * b3mag);
        const auto mag12inv  = (floatingpoint)1.0 / (b1mag * b2mag);
        const auto mag32inv  = (floatingpoint)1.0 / (b3mag * b2mag);

        const auto dot13 = dot(b1, b3);
        const auto dot12 = dot(b1, b2);
        const auto dot32 = dot(b3, b2);

        // Bond angles
        const auto cos12     = dot12 * mag12inv;
        const auto sin12_2   = std::max< floatingpoint >(1 - cos12 * cos12, 0.0);
        const auto sin12     = std::max< floatingpoint >(std::sqrt(sin12_2), angSinMin);
        const auto sin12inv  = (floatingpoint)1.0 / sin12;
        const auto sin12inv2 = sin12inv * sin12inv;

        const auto cos32     = dot32 * mag32inv;
        const auto sin32_2   = std::max< floatingpoint >(1 - cos32 * cos32, 0.0);
        const auto sin32     = std::max< floatingpoint >(std::sqrt(sin32_2), angSinMin);
        const auto sin32inv  = (floatingpoint)1.0 / sin32;
        const auto sin32inv2 = sin32inv * sin32inv;

        const auto cos13 = dot13 * mag13inv;

        // ct -- cos(theta)
        const auto ctFac1 = cos13 - cos12 * cos32; // ct numerator
        const auto ctFac2 = sin12inv * sin32inv;   // ct inv denominator
        auto       ct     = ctFac1 * ctFac2;

        if(ct > 1.0 + cosDihTol || ct < -1.0 - cosDihTol) {
            LOG(WARNING) << "Dihedral incorrect cos theta: " << ct;
            LOG(INFO) << "Interaction information:\n"
                << "Bead coords: " << coord1 << ' ' << coord2 << ' ' << coord3 << ' ' << coord4 << '\n'
                << "Position: " << p << ", brancher coord: " << mp;
        }

        if(ct < -1.0) ct = -1.0;
        if(ct >  1.0) ct =  1.0;

        // derivative of the energy E = k theta^2
        // dE = eFac * d(ct), where eFac = -2 k theta / sin(theta)
        const auto theta = std::acos(ct);
        const auto st    = std::max< floatingpoint >(std::sin(theta), dihSinMin);
        const auto stInv = (floatingpoint)1.0 / st;
        const auto eFac  = -2 * kdih[i] * theta * stInv;

        // derivatives on ct - cosine theta; nu - numerator, de - denominator
        //---------------------------------------------------------------------
        // d(ct_nu) / db1
        //
        //           b3                      b1
        //    = ----------- - cos<b1, b3> --------
        //       |b1| |b3|                 |b1|^2
        //
        //                        b2                                   b1
        //      - cos<b3, b2> ----------- + cos<b3, b2> cos<b1, b2> --------
        //                     |b1| |b2|                             |b1|^2
        //
        // d(ct_nu) / db2
        //
        //                         b3                                  b2
        //    = - cos<b1, b2> ----------- + cos<b1, b2> cos<b3, b2> --------
        //                     |b2| |b3|                             |b2|^2
        //
        //                         b1                                  b2
        //      - cos<b3, b2> ----------- + cos<b1, b2> cos<b3, b2> --------
        //                     |b2| |b1|                             |b2|^2
        //
        // d(ln ct_de) / db1
        //
        //         cos^2 <b1, b2>  [     b2         b1    ]
        //    = - ---------------- [ --------- - -------- ]
        //         sin^2 <b1, b2>  [  b1 . b2     |b1|^2  ]
        //
        // d(ln ct_de) / db2
        //
        //         cos^2 <b3, b2>  [     b3         b2    ]
        //    = - ---------------- [ --------- - -------- ]
        //         sin^2 <b3, b2>  [  b3 . b2     |b2|^2  ]
        //
        //         cos^2 <b1, b2>  [     b1         b2    ]
        //      - ---------------- [ --------- - -------- ]
        //         sin^2 <b1, b2>  [  b1 . b2     |b2|^2  ]
        //
        //---------------------------------------------------------------------
        //          d(ct_nu)
        // d(ct) = ---------- - ct * d(ln ct_de)
        //            ct_de
        //
        // Let E1i = eFac * (cti1 * b1 + cti2 * b2 + cti3 * b3), then
        //
        //             ct       cos^2 <b1, b2>       ct                  ct
        // ct11 = - -------- - ---------------- * -------- = - -----------------------
        //           |b1|^2     sin^2 <b1, b2>     |b1|^2       sin^2 <b1, b2> |b1|^2
        //
        //             cos<b3, b2>       cos^2 <b1, b2>   ct
        // ct12 = - ----------------- + -------------------------
        //           ct_de |b1| |b2|     sin^2 <b1, b2>  b1 . b2
        //
        //             1      [    cos<b3, b2>     cos<b1, b2> ct  ]
        //      = ----------- [ - ------------- + ---------------- ]
        //         |b1| |b2|  [       ct_de        sin^2 <b1, b2>  ]
        //
        // ct13 = 1 / (|b1| |b3| ct_de)
        //
        //             cos<b3, b2>       cos^2 <b1, b2>   ct
        // ct21 = - ----------------- + ------------------------- = ct12
        //           ct_de |b1| |b2|     sin^2 <b1, b2>  b1 . b2
        //
        //         2 cos<b1, b2> cos<b3, b2>      cos^2 <b1, b2>  ct        cos^2 <b3, b2>  ct
        // ct22 = --------------------------- - ----------------------- - -----------------------
        //              ct_de   |b2|^2           sin^2 <b1, b2> |b2|^2     sin^2 <b3, b2> |b2|^2
        //
        //            1    [  2 cos<b1, b3>       [      cos^2 <b1, b2>     cos^2 <b3, b2>  ] ]
        //      = -------- [ --------------- - ct [ 2 + ---------------- + ---------------- ] ]
        //         |b2|^2  [      ct_de           [      sin^2 <b1, b2>     sin^2 <b3, b2>  ] ]
        //
        //            1    [  2 cos<b1, b3>       [         1                  1        ] ]
        //      = -------- [ --------------- - ct [ ---------------- + ---------------- ] ]
        //         |b2|^2  [      ct_de           [  sin^2 <b1, b2>     sin^2 <b3, b2>  ] ]
        //
        //             cos<b1, b2>       cos^2 <b3, b2>   ct
        // ct23 = - ----------------- + ------------------------- = ct32
        //           ct_de |b3| |b2|     sin^2 <b3, b2>  b3 . b2
        //
        //---------------------------------------------------------------------

        const auto e111 = -eFac * ct * sin12inv2 * b1mag2inv;
        const auto e112 = eFac * mag12inv * (-cos32 * ctFac2 + cos12 * ct * sin12inv2);
        const auto e113 = eFac * mag13inv * ctFac2;

        const auto e121 = e112;
        const auto e122 = eFac * b2mag2inv * (2 * cos13 * ctFac2 - ct * (sin12inv2 + sin32inv2));
        const auto e123 = eFac * mag32inv * (-cos12 * ctFac2 + cos32 * ct * sin32inv2);

        const auto e131 = e113;
        const auto e132 = e123;
        const auto e133 = -eFac * ct * sin32inv2 * b3mag2inv;

        // translate to actual force dependence on b1, b2 and b3
        // fij is factor for the jth component (b1, b2, b3) of force acting on bead i.
        const auto f11 = (1-p) * (e111 + e121);
        const auto f12 = (1-p) * (e112 + e122);
        const auto f13 = (1-p) * (e113 + e123);

        const auto f21 = (p-1) * e111 + p * e121;
        const auto f22 = (p-1) * e112 + p * e122;
        const auto f23 = (p-1) * e113 + p * e123;

        const auto f31 = -e121 + e131;
        const auto f32 = -e122 + e132;
        const auto f33 = -e123 + e133;

        const auto f41 = -e131;
        const auto f42 = -e132;
        const auto f43 = -e133;

        Vec<3,floatingpoint> force1, force2, force3, force4;
        //U2(c1,c2prime,c3,c4)<=>Utilda(c1,c2,c3,c4)
        if(areEqual(pos[i],(floatingpoint) 1.0)){
            //In this scenario, the energy is defined as U2(c1, c2prime, c3, c4). We are
            // trying to get U(c1, c2, c3, c4) from it.
            //c1-parent minusend | c2 - parent plusend/bindingsite | c2prime-parent
            // extendedplusend   | c3 - offspring minusend         | c4 - offspring plusend
            //[dU(c1,c2,c3,c4)]             [dU2(c1,c2prime,c3,c4)]   dc2prime
            //[---------------]        =    [---------------------] x --------
            //[    dc2        ]c1,c3,c4     [      dc2prime       ]     dc2
            //______________________________________________________________________________
            //[dU(c1,c2,c3,c4)]          [   [dU2(c1,c2prime,c3,c4)]   dc2prime
            //[---------------]        = [   [---------------------] x --------
            //[    dc1        ]c2,c3,c4  [   [       c2prime       ]     dc1
            //                           [
            //                           [        [dU2(c1,c2prime,c3,c4)]
            //                           [  +     [---------------------]
            //                           [        [         dc1         ]
            //We define c2 = c1 + s(c2prime - c1)
            // dc2prime    s - 1  | dc2prime    1
            // -------- = ------- | -------- = ---
            //   dc1         s    |   dc2       s
            force1 = f11 * b1 + f12 * b2 + f13 * b3;
            const auto force2prime = f21 * b1 + f22 * b2 + f23 * b3;
            force3 = f31 * b1 + f32 * b2 + f33 * b3;
            force4 = f41 * b1 + f42 * b2 + f43 * b3;
            const auto factor1 = (p-1)/p;
            const auto factor2 = (1/p);
            force1 += force2prime*factor1;
            force2 = force2prime*factor2;
        }
        // Default case. compute forces U(c1, c2, c3, c4)<=>Utilda(c2,mp,c3,c4)
        else {
            force1 = f11 * b1 + f12 * b2 + f13 * b3;
            force2 = f21 * b1 + f22 * b2 + f23 * b3;
            force3 = f31 * b1 + f32 * b2 + f33 * b3;
            force4 = f41 * b1 + f42 * b2 + f43 * b3;
        }

        // apply forces
        f1 += force1;
        f2 += force2;
        f3 += force3;
        f4 += force4;

        if(stretchforce) {
            for(short j = 0; j < 3; j++)
                stretchforce[3*i + j] = force3[j];
        }

    } // End loop interactions

} // void forces(...)

} // namespace medyan

#ifndef MEDYAN_Mechanics_ForceField_Membrane_MembraneMeshless_hpp
#define MEDYAN_Mechanics_ForceField_Membrane_MembraneMeshless_hpp

#include <algorithm> // fill
#include <vector>

#include "Structure/DofSerializer.hpp"
#include "Mechanics/ForceField/ForceField.h"

namespace medyan {

// Shiba, Hayato and Noguchi, Hiroshi. "Estimation of the bending rigidity and spontaneous curvature of fluid membranes in simulations" Phys. Rev. E, 2011.
struct MembraneMeshlessFF : ForceField {
    using Coord3Type    = Eigen::Matrix< floatingpoint, 3, 1 >;

    // Vectorized vertex info.
    struct VertexPairInfo {
        // coordIndex can be directly computed from loopIndex.
        Index coordIndex1 = 0;
        Index coordIndex2 = 0;
        Index loopIndex1 = 0;
        Index loopIndex2 = 0;

        // Temporary states, used by energy/force computation.
        floatingpoint s = 0;
    };

    struct VertexInfo {
        floatingpoint rho = 0;

        floatingpoint dRho_enAtt = 0;
    };

    // Temperature.
    double kT = ::kT;
    // Length scale of vertices, in medyan units.
    double sigma = 1;
    // Used in repulsion and attraction.
    double epsilon = 4;
    // Used in repulsion.
    double B = 0.126;
    // Used in attraction.
    double rhoStar = 6;
    // Used in repulsion.
    double scutRep = 1.2;
    // Used in attraction.
    double scutAtt = 2.1;
    // Used in orientation.
    double scc = 3;
    // Used in neighbor list cutoff during vectorization.
    double maxscut = std::max({ scutRep, scutAtt, scc });

    // Neighbor list and vectorization information.
    std::vector< VertexPairInfo > vertexPairsInfo;

    // Temporary storage for densities.
    std::vector< VertexInfo > verticesInfo;

    virtual std::string getName() override { return "MembraneMeshless"; }

    virtual void vectorize(const FFCoordinateStartingIndex& si) override {
        const auto& cl = si.ps->meshlessSpinVertexCellList;

        vertexPairsInfo.clear();
        const auto numCells = cl.numCells();
        for(int i = 0; i < numCells; ++i) {
            // vi, vj are the vertex system indices.
            cl.forAllPairs([&](int vi, int vj) {
                auto& v1 = si.ps->meshlessSpinVertices[vi];
                auto& v2 = si.ps->meshlessSpinVertices[vj];
                // Exclude out-of-range pairs.
                if((v1.coord - v2.coord).squaredNorm() < maxscut * maxscut * sigma * sigma) {
                    vertexPairsInfo.push_back({
                        findMeshlessSpinVertexCoordIndex(v1, si),
                        findMeshlessSpinVertexCoordIndex(v2, si),
                        v1.loopIndex,
                        v2.loopIndex,
                    });
                }
            });
        }

        verticesInfo.assign(si.ps->meshlessSpinVertices.size(), VertexInfo{});
    }

    virtual floatingpoint computeEnergy(floatingpoint *coord) override {
        floatingpoint cbd = 0;
        floatingpoint ktilt = 1;
        floatingpoint kbend = 1;

        floatingpoint enRep = 0;
        floatingpoint enAtt = 0;
        floatingpoint enTilt = 0;
        floatingpoint enBend = 0;

        // Clear temporary variables.
        std::fill(verticesInfo.begin(), verticesInfo.end(), VertexInfo{});

        // 1st pass, find all s and densities.
        for(auto& pairInfo : vertexPairsInfo) {
            auto c1 = Coord3Type::Map(coord + pairInfo.coordIndex1);
            auto c2 = Coord3Type::Map(coord + pairInfo.coordIndex2);
            auto r = (c1 - c2).norm();
            pairInfo.s = r / sigma;
            auto dens = densityPairAtt(pairInfo.s);
            verticesInfo[pairInfo.loopIndex1].rho += dens;
            verticesInfo[pairInfo.loopIndex2].rho += dens;
        }

        // 2nd pass, compute all energies except for attraction.
        for(const auto& pairInfo : vertexPairsInfo) {
            auto c1 = Coord3Type::Map(coord + pairInfo.coordIndex1);
            auto c2 = Coord3Type::Map(coord + pairInfo.coordIndex2);
            auto r12Normalized = (c2 - c1).normalized();

            const auto theta1 = coord[pairInfo.coordIndex1 + 3];
            const auto phi1   = coord[pairInfo.coordIndex1 + 4];
            const auto ct1    = std::cos(theta1);
            const auto st1    = std::sin(theta1);
            const auto cp1    = std::cos(phi1);
            const auto sp1    = std::sin(phi1);
            const auto theta2 = coord[pairInfo.coordIndex2 + 3];
            const auto phi2   = coord[pairInfo.coordIndex2 + 4];
            const auto ct2    = std::cos(theta2);
            const auto st2    = std::sin(theta2);
            const auto cp2    = std::cos(phi2);
            const auto sp2    = std::sin(phi2);
            Coord3Type u1{}; u1 << cp1 * st1, sp1 * st1, ct1;
            Coord3Type u2{}; u2 << cp2 * st2, sp2 * st2, ct2;

            enRep += energyPairRep(pairInfo.s);
            enTilt += energyPairTilt(r12Normalized, pairInfo.s, u1, u2, ktilt);
            enBend += energyPairBend(r12Normalized, pairInfo.s, u1, u2, kbend, cbd);
        }

        // Loop through vertices to find attraction energy.
        for(const auto& vInfo : verticesInfo) {
            enAtt += energySingleAtt(vInfo.rho);
        }

        return enRep + enAtt + enTilt + enBend;
    }

    virtual void computeForces(floatingpoint *coord, floatingpoint *force) override {
        floatingpoint cbd = 0;
        floatingpoint ktilt = 1;
        floatingpoint kbend = 1;

        // Clear temporary variables.
        std::fill(verticesInfo.begin(), verticesInfo.end(), VertexInfo{});

        // 1st pass, find all s and densities.
        for(auto& pairInfo : vertexPairsInfo) {
            auto c1 = Coord3Type::Map(coord + pairInfo.coordIndex1);
            auto c2 = Coord3Type::Map(coord + pairInfo.coordIndex2);
            auto r = (c1 - c2).norm();
            pairInfo.s = r / sigma;
            auto dens = densityPairAtt(pairInfo.s);
            verticesInfo[pairInfo.loopIndex1].rho += dens;
            verticesInfo[pairInfo.loopIndex2].rho += dens;
        }

        // Find derivatives on attraction.
        for(auto& vInfo : verticesInfo) {
            vInfo.dRho_enAtt = dRho_energySingleAtt(vInfo.rho);
        }

        // 2nd pass, compute all forces.
        for(const auto& pairInfo : vertexPairsInfo) {
            auto c1 = Coord3Type::Map(coord + pairInfo.coordIndex1);
            auto c2 = Coord3Type::Map(coord + pairInfo.coordIndex2);
            auto f1 = Coord3Type::Map(force + pairInfo.coordIndex1);
            auto f2 = Coord3Type::Map(force + pairInfo.coordIndex2);
            const auto r = (c2 - c1).norm();
            const Coord3Type r12Normalized = (c2 - c1) / r;
            // dc1_r = -r12Normalized.
            // dc2_r =  r12Normalized.

            const auto theta1  = coord[pairInfo.coordIndex1 + 3];
            const auto phi1    = coord[pairInfo.coordIndex1 + 4];
            auto&      ftheta1 = force[pairInfo.coordIndex1 + 3];
            auto&      fphi1   = force[pairInfo.coordIndex1 + 4];
            const auto ct1     = std::cos(theta1);
            const auto st1     = std::sin(theta1);
            const auto cp1     = std::cos(phi1);
            const auto sp1     = std::sin(phi1);
            const auto theta2  = coord[pairInfo.coordIndex2 + 3];
            const auto phi2    = coord[pairInfo.coordIndex2 + 4];
            auto&      ftheta2 = force[pairInfo.coordIndex2 + 3];
            auto&      fphi2   = force[pairInfo.coordIndex2 + 4];
            const auto ct2     = std::cos(theta2);
            const auto st2     = std::sin(theta2);
            const auto cp2     = std::cos(phi2);
            const auto sp2     = std::sin(phi2);
            Coord3Type u1{}; u1 << cp1 * st1, sp1 * st1, ct1;
            Coord3Type u2{}; u2 << cp2 * st2, sp2 * st2, ct2;

            const auto ds_enRep = ds_energyPairRep(pairInfo.s);
            const auto ds_dens  = ds_densityPairAtt(pairInfo.s);
            const auto ds_wg    = ds_weightGyration(pairInfo.s);

            // Repulsion.
            f1 -= (-ds_enRep / sigma) * r12Normalized;
            f2 -= ( ds_enRep / sigma) * r12Normalized;

            // Attraction.
            const auto sum_dRho_enAtt = verticesInfo[pairInfo.loopIndex1].dRho_enAtt +
                                        verticesInfo[pairInfo.loopIndex2].dRho_enAtt;
            f1 -= (- sum_dRho_enAtt * ds_dens / sigma) * r12Normalized;
            f2 -= (  sum_dRho_enAtt * ds_dens / sigma) * r12Normalized;

            // Tilt.
            {
                const auto dot1 = r12Normalized.dot(u1);
                const auto dot2 = r12Normalized.dot(u2);
                const auto enTiltPart1 = kT * ktilt / 2;
                const auto enTiltPart2 = dot1 * dot1 + dot2 * dot2;
                const auto enTiltPart3 = weightGyration(pairInfo.s);

                // Part 2.
                Coord3Type g2 = enTiltPart1 * enTiltPart3 * (
                    (2 * dot1 / r) * (- dot1 * r12Normalized + u1) +
                    (2 * dot2 / r) * (- dot2 * r12Normalized + u2)
                );
                // g1 = -g2.
                auto gtheta1 = enTiltPart1 * enTiltPart3 * 2 * dot1 * (
                    + r12Normalized[0] * cp1 * ct1
                    + r12Normalized[1] * sp1 * ct1
                    - r12Normalized[2] * st1
                );
                auto gphi1 = enTiltPart1 * enTiltPart3 * 2 * dot1 * (
                    - r12Normalized[0] * sp1 * st1
                    + r12Normalized[1] * cp1 * st1
                );
                auto gtheta2 = enTiltPart1 * enTiltPart3 * 2 * dot2 * (
                    + r12Normalized[0] * cp2 * ct2
                    + r12Normalized[1] * sp2 * ct2
                    - r12Normalized[2] * st2
                );
                auto gphi2 = enTiltPart1 * enTiltPart3 * 2 * dot2 * (
                    - r12Normalized[0] * sp2 * st2
                    + r12Normalized[1] * cp2 * st2
                );

                // Part 3.
                g2 += (enTiltPart1 * enTiltPart2 * ds_wg / sigma) * r12Normalized;

                // Conclude.
                f1 += g2;
                f2 -= g2;
                ftheta1 -= gtheta1;
                fphi1   -= gphi1;
                ftheta2 -= gtheta2;
                fphi2   -= gphi2;
            }

            // Bend.
            {
                const auto enBendPart1 = kT * kbend / 2;
                const Coord3Type diff = u1 - u2;
                const auto enBendPart2 = (diff - cbd * r12Normalized).squaredNorm();
                const auto enBendPart3 = weightGyration(pairInfo.s);

                // Part 2.
                const auto diffProj = diff.dot(r12Normalized);
                Coord3Type g2 = enBendPart1 * enBendPart3
                    * (-2 * cbd / r)
                    * (- diffProj * r12Normalized + diff);
                // g1 = -g2.
                const Coord3Type dotBy1 = 2 * (- u2 - cbd * r12Normalized);
                auto gtheta1 = enBendPart1 * enBendPart3 * (
                    + dotBy1[0] * cp1 * ct1
                    + dotBy1[1] * sp1 * ct1
                    - dotBy1[2] * st1
                );
                auto gphi1 = enBendPart1 * enBendPart3 * (
                    - dotBy1[0] * sp1 * st1
                    + dotBy1[1] * cp1 * st1
                );
                const Coord3Type dotBy2 = -2 * (u1 - cbd * r12Normalized);
                auto gtheta2 = enBendPart1 * enBendPart3 * (
                    + dotBy2[0] * cp2 * ct2
                    + dotBy2[1] * sp2 * ct2
                    - dotBy2[2] * st2
                );
                auto gphi2 = enBendPart1 * enBendPart3 * (
                    - dotBy2[0] * sp2 * st2
                    + dotBy2[1] * cp2 * st2
                );

                // Part 3.
                g2 += (enBendPart1 * enBendPart2 * ds_wg / sigma) * r12Normalized;

                // Conclude.
                f1 += g2;
                f2 -= g2;
                ftheta1 -= gtheta1;
                fphi1   -= gphi1;
                ftheta2 -= gtheta2;
                fphi2   -= gphi2;
            }
        }

    }



    // Functions to compute energies and forces.
    //----------------------------------

    // Precondition:
    // - s is positive.
    floatingpoint factorCut(floatingpoint s, floatingpoint a, floatingpoint scut, int n) const {
        return s < scut ? std::exp(a * (1 + 1 / (std::pow(s / scut, n) - 1))) : 0;
    }
    floatingpoint ds_factorCut(floatingpoint s, floatingpoint a, floatingpoint scut, int n) const {
        if(s < scut) {
            auto spow = std::pow(s / scut, n);
            auto invspowm1 = 1 / (spow - 1);
            return std::exp(a * (1 + invspowm1))
                * a
                * (- invspowm1 * invspowm1)
                * (n * spow / s);
        } else {
            return 0;
        }
    }

    floatingpoint weightGyration(floatingpoint s) const {
        const auto inv_scc = 1.0 / scc;
        const auto inv_sga = 1.0 / 1.5;
        const auto n = 12;
        return s < scc ? std::exp((s * inv_sga) * (s * inv_sga) / (std::pow(s * inv_scc, n) - 1)) : 0;
    }
    floatingpoint ds_weightGyration(floatingpoint s) const {
        if(s < scc) {
            const auto inv_scc = 1.0 / scc;
            const auto inv_sga = 1.0 / 1.5;
            const auto n = 12;
            auto spow = std::pow(s * inv_scc, n);
            auto invspowm1 = 1 / (spow - 1);
            return std::exp((s * inv_sga) * (s * inv_sga) * invspowm1)
                * (
                    invspowm1 * 2 * s * inv_sga * inv_sga
                    + (s * inv_sga) * (s * inv_sga)
                        * (- invspowm1 * invspowm1)
                        * (n * spow / s)
                );
        } else {
            return 0;
        }
    }

    // Precondition:
    // - s = rij / σ, is positive.
    floatingpoint energyPairRep(floatingpoint s) const {
        return kT * epsilon
            * std::exp(-20 * (s - 1) + B)
            * factorCut(s, 1, scutRep, 12);
    }
    floatingpoint ds_energyPairRep(floatingpoint s) const {
        auto expPart = std::exp(-20 * (s - 1) + B);
        return kT * epsilon * expPart * (
            (-20) * factorCut(s, 1, scutRep, 12)
            + ds_factorCut(s, 1, scutRep, 12)
        );
    }

    floatingpoint densityPairAtt(floatingpoint s) const {
        return factorCut(s, 3.715, scutAtt, 12);
    }
    floatingpoint ds_densityPairAtt(floatingpoint s) const {
        return ds_factorCut(s, 3.715, scutAtt, 12);
    }

    floatingpoint energySingleAtt(floatingpoint rho) const {
        // C = 0.25 * std::log(1 + std::exp(4 * rhoStar)) ≈ rhoStar
        return kT * epsilon * (0.25 * std::log(1 + std::exp(-4 * (rho - rhoStar))) - rhoStar);
    }
    floatingpoint dRho_energySingleAtt(floatingpoint rho) const {
        const auto expPart = std::exp(-4 * (rho - rhoStar));
        return kT * epsilon * (-expPart / (1 + expPart));
    }

    floatingpoint energyPairTilt(const Coord3Type& rijNormalized, floatingpoint s, const Coord3Type& ui, const Coord3Type& uj, floatingpoint ktilt) const {
        const auto dot1 = rijNormalized.dot(ui);
        const auto dot2 = rijNormalized.dot(uj);
        return kT * ktilt / 2
            * (dot1 * dot1 + dot2 * dot2)
            * weightGyration(s);
    }
    // cbd indicates spontaneous curvature.
    floatingpoint energyPairBend(const Coord3Type& rijNormalized, floatingpoint s, const Coord3Type& ui, const Coord3Type& uj, floatingpoint kbend, floatingpoint cbd) const {
        return kT * kbend / 2
            * (ui - uj - cbd * rijNormalized).squaredNorm()
            * weightGyration(s);
    }
};

} // namespace medyan

#endif

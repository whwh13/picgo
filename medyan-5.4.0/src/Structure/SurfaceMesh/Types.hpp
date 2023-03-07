#ifndef MEDYAN_Structure_SurfaceMesh_Types_hpp
#define MEDYAN_Structure_SurfaceMesh_Types_hpp

#include <stdexcept>
#include <string>

#include "Util/Parser/StringParser.hpp"

namespace medyan {

// If vertices represent fixed coordinates in a material coordinate system,
// a vertex represents the location of a specific lipid molecule.
// The local area elasticity is usually a necessity in this case.
// An exception is the vertices on the border connected to a lipid
// reservoir, where the vertices stand for the boundary locations which are
// usually fixed in the ambient space.
//
// If vertices represent fixed coordinates in a normal coordinate system,
// a vertex is only a representative point on the surface, where the motion
// of the vertex must be in the local normal direction.
// Local area elasticity cannot be directly defined on mesh elements.
//
// If vertices represent fixed coordinates in a normal coordinate system,
// then the system is similar to the "normal" coordinate case, but without
// the requirement for vertices to move only in the normal directions.
enum class MembraneMeshVertexSystem { material, normal, general };

template<>
struct StringSerializerTrait<MembraneMeshVertexSystem> {
    auto parse(std::string_view sv) const {
        if (sv == "material") {
            return MembraneMeshVertexSystem::material;
        } else if (sv == "normal") {
            return MembraneMeshVertexSystem::normal;
        } else if (sv == "general") {
            return MembraneMeshVertexSystem::general;
        } else {
            throw std::invalid_argument("Invalid membrane mesh vertex system");
        }
    }

    std::string toString(MembraneMeshVertexSystem val) const {
        switch(val) {
            case MembraneMeshVertexSystem::material: return "material";
            case MembraneMeshVertexSystem::normal:   return "normal";
            case MembraneMeshVertexSystem::general:  return "general";
            default:                                 return "";
        }
    }
};


// Specifies how membrane vertices are pinned.
// Pinned vertices should not move during the simulation.
enum class MembraneVertexPinning {
    none,
    // Pin border vertices.
    border1,
    // Pin border vertices and all their neighbors.
    border2,
};

template<>
struct StringSerializerTrait<MembraneVertexPinning> {
    auto parse(std::string_view sv) const {
        if (sv == "none") {
            return MembraneVertexPinning::none;
        } else if (sv == "border1") {
            return MembraneVertexPinning::border1;
        } else if (sv == "border2") {
            return MembraneVertexPinning::border2;
        } else {
            throw std::invalid_argument("Invalid membrane vertex pinning");
        }
    }

    std::string toString(MembraneVertexPinning val) const {
        switch(val) {
            case MembraneVertexPinning::none:    return "none";
            case MembraneVertexPinning::border1: return "border1";
            case MembraneVertexPinning::border2: return "border2";
            default:                             return "";
        }
    }
};


enum class SurfaceCurvaturePolicy {
    // Computes the signed curvature.
    //
    // This is generally needed when curvature is used in linear cases,
    // such as the bending energy when spontaneous curvature is non-zero.
    //
    // The default implementation is
    //         ∇ Area ⋅ ∇ Vol
    //   H = ------------------
    //        2 ∇ Vol ⋅ ∇ Vol
    //
    // where ∇ acts on the vertex of interest.
    //
    // The signed curvature will be stored in the curv variable, and their
    // derivatives will be stored in curv related places.
    withSign,

    // Computes the squared curvature.
    //
    // This is generally used when curvature is used in quadratic cases,
    // such as the bending energy when spontaneous curvature is zero.
    //
    // The default implementation is
    //           (∇ Area)^2
    //   H^2 = -------------
    //          4 (∇ Vol)^2
    //
    // The square of the curvature will be stored in the curv2 variable, and
    // their derivatives will be stored related to curv2.
    squared,

    // Described in the paper "Mem3DG: Modeling Membrane Mechanochemical Dynamics in 3D using Discrete Differential Geometry", by C. Zhu et al. (2021) https://doi.org/10.1101/2021.10.30.466618.
    //
    // In equation (10) of the preprint,
    //                  l_ij ϕ_ij
    //   H_i = sum_ij  -----------
    //                    4 A_i
    //
    // where l_ij is the length of the edge ij, ϕ_ij is the dihedral angle between the two faces of the edge ij, and A_i is the area around vertex i.
    //
    // Using this policy should reproduce the credit of the paper.
    mem3dg,
};

template<>
struct StringSerializerTrait<SurfaceCurvaturePolicy> {
    auto parse(std::string_view sv) const {
        if (sv == "with-sign") {
            return SurfaceCurvaturePolicy::withSign;
        } else if (sv == "squared") {
            return SurfaceCurvaturePolicy::squared;
        } else if (sv == "mem3dg") {
            return SurfaceCurvaturePolicy::mem3dg;
        } else {
            throw std::invalid_argument("Invalid surface curvature policy");
        }
    }

    std::string toString(SurfaceCurvaturePolicy val) const {
        switch(val) {
            case SurfaceCurvaturePolicy::withSign: return "with-sign";
            case SurfaceCurvaturePolicy::squared:  return "squared";
            case SurfaceCurvaturePolicy::mem3dg:   return "mem3dg";
            default:                               return "unknown";
        }
    }
};


struct SurfaceGeometrySettings {
    SurfaceCurvaturePolicy curvPol = SurfaceCurvaturePolicy::withSign;
};


// Parameters controlling remeshing.
struct MeshAdapterSettings {
    // Connectivity.
    Size   minDegree                         = 4;
    Size   maxDegree                         = 9;
    double edgeFlipMinDotNormal              = 0.9;
    double edgeCollapseMinQualityImprovement = 0.6;
    double edgeCollapseMinDotNormal          = 0.85;

    // Relaxation.
    Size   relaxationMaxIterRelocation       = 10;
    Size   relaxationMaxIterTotal            = 3; // (vertex relocation + edge flipping) as 1 iter

    // Size diffusion.
    double curvatureResolution               = 0.30; // cos of which should be slightly bigger than flip minDotNormal.
    double maxSize                           = 60; // Related to the resolution of the system.
    Size   diffuseIter                       = 4;

    // Main loop.
    Size   samplingAdjustmentMaxIter         = 10; // Max times of scanning all the edges for sampling adjustment.
    Size   mainLoopSoftMaxIter               = 8;
    Size   mainLoopHardMaxIter               = 16;
};

} // namespace medyan

#endif

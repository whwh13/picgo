#ifndef MEDYAN_Mechanics_ForceField_Membrane_MembraneTriangleProtect_hpp
#define MEDYAN_Mechanics_ForceField_Membrane_MembraneTriangleProtect_hpp

#include <array>
#include <type_traits>
#include <vector>

#include "Mechanics/ForceField/ForceField.h"
#include "Mechanics/ForceField/Types.hpp"
#include "Structure/SubSystem.h"
#include "Structure/SurfaceMesh/FuncMembraneGeo.hpp"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Util/Io/Log.hpp"
#include "Util/Math/Vec.hpp"

namespace medyan {

struct MembraneTriangleProtectInteraction {
    using VecDArea = decltype(GHalfEdge::dTriangleArea);

    const double*                  parea = nullptr;
    double                         initarea = 0;
    std::array<const VecDArea*, 3> pdarea {};
    std::array<Index, 3>           coordIndices {};
};

struct MembraneTriangleProtectFene {

    static FP energy(const MembraneTriangleProtectInteraction& in, double minRelArea, double k) {
        const auto relArea = *in.parea / (in.initarea * minRelArea);
        if(relArea >= 1.0) return 0.0;

        const auto relDiff = relArea - 1;
        return -0.5 * k * std::log(1 - relDiff * relDiff);
    }

    static void forces(FP* force, const MembraneTriangleProtectInteraction& in, double minRelArea, double k) {
        const auto& area = *in.parea;
        const auto minArea = in.initarea * minRelArea;
        if(area >= minArea) return;

        const auto relDiff = area / minArea - 1;
        const auto deda = k * relDiff / ((1 - relDiff * relDiff) * minArea);

        for(int i = 0; i < 3; ++i) {
            makeRefVec<3>(force + in.coordIndices[i]) -= deda * (*in.pdarea[i]);
        }
    }
};

template< typename Impl, bool warn = false >
struct MembraneTriangleProtect : ForceField {

    // Temporary topological information
    //---------------------------------
    std::vector<MembraneTriangleProtectInteraction>   ins;

    // Parameters
    //---------------------------------
    const double k           = 2000.0;
    const double minRelArea  = 0.35;
    const double warnRelArea = 0.10;

    virtual void vectorize(const FFCoordinateStartingIndex& si) override {
        using MT = Membrane::MeshType;

        ins.clear();

        for(const auto& m : si.ps->membranes) {
            const auto& mesh = m.getMesh();
            assertValidIndexCacheForFF(mesh);

            for(const auto& t : mesh.getTriangles()) {
                auto& in = ins.emplace_back();

                in.parea    = &t.attr.gTriangle.area;
                in.initarea = t.attr.gTriangle.area;

                Index bi = 0;
                mesh.forEachHalfEdgeInPolygon(t, [&](MT::HalfEdgeIndex hei) {
                    const auto vi = mesh.target(hei);
                    in.coordIndices[bi] = mesh.attribute(vi).cachedCoordIndex;
                    in.pdarea[bi]       = &mesh.attribute(hei).gHalfEdge.dTriangleArea;
                    ++bi;
                });
            }
        }
    }

    virtual FP computeEnergy(FP* coord) override {
        FP e = 0.0;

        for(Index ti = 0; ti < ins.size(); ++ti) {
            const auto& in = ins[ti];

            if(warn && *in.parea < in.initarea * warnRelArea) {
                log::warn("Triangle {} size becomes too small: area={} initarea={}", ti, *in.parea, in.initarea);
            }

            e += Impl::energy(in, minRelArea, k);
        }

        return e;
    }
    virtual void computeForces(FP* coord, FP* force) override {
        for(Index ti = 0; ti < ins.size(); ++ti) {
            const auto& in = ins[ti];

            if(warn && *in.parea < in.initarea * warnRelArea) {
                log::warn("Triangle {} size becomes too small: area={} initarea={}", ti, *in.parea, in.initarea);
            }

            Impl::forces(force, in, minRelArea, k);
        }
    }

    virtual std::string getName() override { return "MembraneTriangleProtect"; }

};

} // namespace medyan

#endif

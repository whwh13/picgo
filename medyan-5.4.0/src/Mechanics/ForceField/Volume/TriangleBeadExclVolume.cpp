#include "Mechanics/ForceField/Volume/TriangleBeadExclVolume.hpp"

#include <algorithm> // max

#include "MathFunctions.h"
#include "Mechanics/ForceField/Volume/TriangleBeadExclVolRepulsion.hpp"
#include "Structure/Bead.h"
#include "Structure/Cylinder.h"
#include "Structure/Filament.h"
#include "Structure/SubSystem.h"
#include "Structure/SurfaceMesh/FuncMembraneGeo.hpp"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Structure/SurfaceMesh/Triangle.hpp"
#include "Structure/SurfaceMesh/Vertex.hpp"
#include "Util/Math/RayTriangleIntersect.hpp"

namespace medyan {
using namespace mathfunc;

FP TriangleBeadExclVolume::computeEnergy(FP* coord) {
    using MT = Membrane::MeshType;

    double en = 0;

    for(auto t : ps->triangles) {

        const auto& mesh = t.getParent(*ps).getMesh();
        medyan::assertValidIndexCacheForFF(mesh);

        const MT::TriangleIndex ti { t.getTopoIndex() };
        const auto& ta = mesh.attribute(ti);
        const auto vi0 = ta.cachedCoordIndex[0];
        const auto vi1 = ta.cachedCoordIndex[1];
        const auto vi2 = ta.cachedCoordIndex[2];

        const auto area = ta.gTriangle.area;
        const double kExVol = SysParams::Mechanics().triangleBeadVolume.k;
        
        if(neighborList_->hasNeighborMech(t.sysIndex)) for(auto b : neighborList_->getNeighborsMech(t.sysIndex)) {

            const double enInter = interaction.energy(
                static_cast<Vec3d>(makeVec<3>(coord + vi0)),
                static_cast<Vec3d>(makeVec<3>(coord + vi1)),
                static_cast<Vec3d>(makeVec<3>(coord + vi2)),
                static_cast<Vec3d>(makeVec<3>(coord + 3 * b->getIndex() + beadStartingIndex)),
                area, kExVol);
            
            if(!std::isfinite(enInter) || enInter < 0.0) {
                
                return inffp;
            }
            else {
                en += enInter;
            }
        }
    }
    
    return en;
}

void TriangleBeadExclVolume::computeForces(FP* coord, FP* force) {
    using MT = Membrane::MeshType;

    for(auto t : ps->triangles) {

        const auto& mesh = t.getParent(*ps).getMesh();
        medyan::assertValidIndexCacheForFF(mesh);

        const MT::TriangleIndex ti { t.getTopoIndex() };
        const auto& ta = mesh.attribute(ti);
        const auto hei0 = ta.cachedHalfEdgeIndex[0];
        const auto hei1 = ta.cachedHalfEdgeIndex[1];
        const auto hei2 = ta.cachedHalfEdgeIndex[2];
        const auto vi0 = ta.cachedCoordIndex[0];
        const auto vi1 = ta.cachedCoordIndex[1];
        const auto vi2 = ta.cachedCoordIndex[2];

        const auto area = ta.gTriangle.area;
        const auto& dArea0 = mesh.attribute(hei0).gHalfEdge.dTriangleArea;
        const auto& dArea1 = mesh.attribute(hei1).gHalfEdge.dTriangleArea;
        const auto& dArea2 = mesh.attribute(hei2).gHalfEdge.dTriangleArea;
        const double kExVol = SysParams::Mechanics().triangleBeadVolume.k;

        if(neighborList_->hasNeighborMech(t.sysIndex)) for(auto b : neighborList_->getNeighborsMech(t.sysIndex)) {

            interaction.forces(
                force + vi0,
                force + vi1,
                force + vi2,
                force + 3 * b->getIndex() + beadStartingIndex,
                static_cast<Vec3>(makeVec<3>(coord + vi0)),
                static_cast<Vec3>(makeVec<3>(coord + vi1)),
                static_cast<Vec3>(makeVec<3>(coord + vi2)),
                static_cast<Vec3>(makeVec<3>(coord + 3 * b->getIndex() + beadStartingIndex)),
                area, dArea0, dArea1, dArea2, kExVol);
        }
    }

}

namespace {

template< typename InteractionType, typename VecType >
void exclVolLoadForce(
    const InteractionType& interaction, double area, double kExVol,
    const Bead& bo, Bead& bd, TriangleBeadExclVolume::LoadForceEnd end,
    const VecType& v0, const VecType& v1, const VecType& v2
) {
    using LoadForceEnd = TriangleBeadExclVolume::LoadForceEnd;

    auto& loadForces = (end == LoadForceEnd::Plus ? bd.loadForcesP : bd.loadForcesM);
    auto& lfi        = (end == LoadForceEnd::Plus ? bd.lfip        : bd.lfim       );

    // This is the direction of polymerization
    const auto dir = normalizedVector(bd.coordinate() - bo.coordinate());

    // Test intersection of the ray with the triangle
    const auto intersectRes = ray_tracing::MollerTrumboreIntersect<>()(
        bd.coordinate(), dir,
        v0, v1, v2
    );

    //array of coordinate values to update
    const auto monSize = SysParams::Geometry().monomerSize   [bd.getType()];
    auto       cylSize = SysParams::Geometry().cylinderNumMon[bd.getType()];

    // Set load force to inf right before intersection.
    if(intersectRes.intersect && intersectRes.t > 0) {
        const int maxCylSize = static_cast<int>(intersectRes.t / monSize);
        if(maxCylSize < cylSize) {
            loadForces[maxCylSize] = inffp;
            cylSize = maxCylSize;
        }
    }

    for (int i = 0; i < cylSize; i++) {

        const auto newCoord = bd.coordinate() + (i * monSize) * dir;

        const auto loadForce = interaction.loadForces(v0, v1, v2, newCoord, area, kExVol);
        const auto effLoadForce = std::max< double >(-dot(dir, loadForce), 0.0);

        loadForces[i] += effLoadForce;
    }
    //reset lfi
    lfi = 0;

} // void exclVolLoadForce(...)

} // namespace (anonymous)

void TriangleBeadExclVolume::computeLoadForces() {
    using MT = Membrane::MeshType;

    for(auto t : ps->triangles) {

        const auto& mesh = t.getParent(*ps).getMesh();
        const MT::TriangleIndex ti { t.getTopoIndex() };
        const auto hei0 = mesh.halfEdge(ti);
        const auto hei1 = mesh.next(hei0);
        const auto hei2 = mesh.next(hei1);
        const Vec3 v0 (mesh.attribute(mesh.target(hei0)).getCoordinate(*ps));
        const Vec3 v1 (mesh.attribute(mesh.target(hei1)).getCoordinate(*ps));
        const Vec3 v2 (mesh.attribute(mesh.target(hei2)).getCoordinate(*ps));

        const auto area = mesh.attribute(ti).gTriangle.area;
        double kExVol = SysParams::Mechanics().triangleBeadVolume.k;
        
        if(neighborList_->hasNeighbor(t.sysIndex)) for(auto b : neighborList_->getNeighbors(t.sysIndex)) {

            // potential acts on second cylinder bead if it is a plus  end
            // potential acts on first  cylinder bead if it is a minus end
            Cylinder& cPlus = *static_cast< Filament* >(b->getParent())->getPlusEndCylinder();
            if(b == cPlus.getSecondBead()) {
                exclVolLoadForce(
                    interaction, area, kExVol,
                    *cPlus.getFirstBead(), *b, LoadForceEnd::Plus,
                    v0, v1, v2
                );
            }

            Cylinder& cMinus = *static_cast< Filament* >(b->getParent())->getMinusEndCylinder();
            if(b == cMinus.getFirstBead()) {
                exclVolLoadForce(
                    interaction, area, kExVol,
                    *cMinus.getSecondBead(), *b, LoadForceEnd::Minus,
                    v0, v1, v2
                );
            }
            
        }
    }

} // void ...::computeLoadForces() const

void TriangleBeadExclVolume::computeLoadForce(SubSystem& sys, Cylinder* c, LoadForceEnd end) const {
    using MT = Membrane::MeshType;

    Bead* tip   = (end == LoadForceEnd::Plus ? c->getSecondBead() : c->getFirstBead());
    Bead* other = (end == LoadForceEnd::Plus ? c->getFirstBead() : c->getSecondBead());

    if(neighborList_->hasNeighbor(tip)) for(auto tiSys : neighborList_->getNeighbors(tip)) {

        const auto& t = sys.triangles[tiSys];
        const auto& mesh = t.getParent(sys).getMesh();
        const MT::TriangleIndex ti { t.getTopoIndex() };
        const auto hei0 = mesh.halfEdge(ti);
        const auto hei1 = mesh.next(hei0);
        const auto hei2 = mesh.next(hei1);
        const Vec3 v0 (mesh.attribute(mesh.target(hei0)).getCoordinate(sys));
        const Vec3 v1 (mesh.attribute(mesh.target(hei1)).getCoordinate(sys));
        const Vec3 v2 (mesh.attribute(mesh.target(hei2)).getCoordinate(sys));

        const auto area = mesh.attribute(ti).gTriangle.area;
        double kExVol = SysParams::Mechanics().triangleBeadVolume.k;

        exclVolLoadForce(
            interaction, area, kExVol,
            *other, *tip, end,
            v0, v1, v2
        );
    }

} // void ...::computeLoadForce(...) const

} // namespace medyan

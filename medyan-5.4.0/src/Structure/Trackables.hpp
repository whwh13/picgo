#ifndef MEDYAN_Structure_Trackables_hpp
#define MEDYAN_Structure_Trackables_hpp

#include "Structure/Bubble.h"
#include "Structure/Special/AFM.h"
#include "Structure/Special/MTOC.h"
#include "Structure/SurfaceMesh/Membrane.hpp"
#include "Structure/SurfaceMesh/MembraneHierarchy.hpp"
#include "Structure/SurfaceMesh/FixedVertexAttachment.hpp"
#include "Structure/SurfaceMesh/FuncMembraneGeo.hpp"

namespace medyan {

namespace trackable {

// The initialize function is used when invoking emplaceTrackable in SubSystemFunc.
// It specifies how a trackable should be initialized, and how other parts the system is affected by addition of this trackable.
//
// Notes:
// - Constructor is not called in this function. The trackable object has already been constructed and is passed to this function.
// - This function is overloaded mainly based on the type of the trackable.
//
// initialize will be called with the following arguments:
// - auto&&     sysFunc:       The subsystem functor object. It is templated to avoid circular dependencies.
// - SubSystem& sys:           The subsystem object.
// - T&         trackable:     The freshly constructed trackable object.
// - auto       index:         The index of the trackable in the associated container. An example is StableVector<T>::Index, when StableVector is used as the container.
// - auto&&...  args:          Other arguments, also used in the trackable's constructor.
constexpr Overload initialize {
    [](
        auto&& sysFunc, SubSystem& sys, Membrane& m, auto index,
        const MembraneSetup&                           memSetup,
        const std::vector< Membrane::CoordinateType >& vertexCoordinateList,
        const std::vector< std::array< int, 3 > >&     triangleVertexIndexList
    ) {
        using MT = Membrane::MeshType;

        auto& mesh = m.getMesh();

        // Build the meshwork topology using vertex and triangle information.
        mesh.init<typename MT::VertexTriangleInitializer>(
            ElementAttributeModifier<SubSystem, std::decay_t<decltype(sysFunc)>>(&sys),
            vertexCoordinateList.size(),
            triangleVertexIndexList,
            Membrane::MeshAttributeType::AttributeInitializerInfo{ vertexCoordinateList }
        );
        mesh.metaAttribute().membraneSysIndex = index;
        mesh.metaAttribute().vertexSystem = memSetup.vertexSystem;
        mesh.metaAttribute().hasLipidReservoir = memSetup.hasLipidReservoir;

        // Update geometry.
        updateGeometryValueForSystem(sys, mesh);

        // Add to membrane hierarchy (must have updated geometry).
        MembraneHierarchy<Membrane>::addMembrane(sys, index);

        // Set parameters already known.
        m.setup             = memSetup;
        m.mMembrane.kArea   = memSetup.areaElasticity;
        m.mMembrane.tension = memSetup.tension;
        m.mMembrane.kVolume = memSetup.volumeElasticity;
    },
    initializeVertex,

    [](
        auto&&, SubSystem&, Edge& e, auto index,
        StableVectorIndex<Membrane> parentSysIndex, Index topoIndex
    ) {
        e.setParentSysIndex(parentSysIndex);
        e.setTopoIndex(topoIndex);
        // Register the stable index inside edge.
        e.sysIndex = index;
    },

    [](
        auto&&, SubSystem&, Triangle& t, auto index,
        StableVectorIndex<Membrane> parentSysIndex, Index topoIndex
    ) {
        t.setParentSysIndex(parentSysIndex);
        t.setTopoIndex(topoIndex);
        // Register the stable index inside triangle.
        t.sysIndex = index;
    },
    initializeMeshlessSpinVertex,
    initializeFixedVertexAttachment,

    [](
        auto&&, SubSystem&, Bubble& b, auto index,
        const Bubble::CoordinateType& coord,
        int                           type,
        const MechParams&             mechParams
    ) {
        b.sysIndex = index;
        b.coord = coord;
        b.setType(type);
        b.setMechanicalProperties(mechParams);
    },
    [](auto&&, SubSystem&, MTOC& obj, auto index) {
        obj.sysIndex = index;
    },
    [](auto&&, SubSystem&, AFM& obj, auto index) {
        obj.sysIndex = index;
    },

    // Other trackables. Future: remove this.
    [](auto&&, SubSystem&, Trackable&, auto, auto&&...) {},
};

// The finalize function is used when invoking removeTrackable in SubSystemFunc.
// It specifies how a trackable should be finalized before removal, and how other parts the system is affected by removal of this trackable.
//
// Notes:
// - Destructor is not called in this function. The trackable object will be removed after calling this function.
// - This function is overloaded mainly based on the type of the trackable.
//
// finalize will be called with the following arguments:
// - auto&&     sysFunc:       The subsystem functor object. It is templated to avoid circular dependencies.
// - SubSystem& sys:           The subsystem object.
// - T&         trackable:     The freshly constructed trackable object.
// - auto       index:         The index of the trackable in the associated container. An example is StableVector<T>::Index, when StableVector is used as the container.
constexpr Overload finalize {
    [](auto&& sysFunc, SubSystem& sys, Membrane& m, auto index) {
        MembraneHierarchy<Membrane>::removeMembrane(index);

        m.getMesh().clear( ElementAttributeModifier<SubSystem, std::decay_t<decltype(sysFunc)>>(&sys) );
    },
    finalizeVertex,
    [](auto&&, SubSystem&, Edge&,     auto) {},
    [](auto&&, SubSystem&, Triangle&, auto) {},
    finalizeMeshlessSpinVertex,
    finalizeFixedVertexAttachment,

    [](auto&&, SubSystem&, Bubble&,   auto) {},
    [](auto&&, SubSystem& sys, MTOC&,     auto) {},
    [](auto&&, SubSystem& sys, AFM& obj,  auto) {
        obj.clearChemistry(sys);
    },

    // Other trackables. Future: remove this.
    [](auto&&, SubSystem&, Trackable&, auto) {},
};

} // namespace trackable

} // namespace medyan

#endif

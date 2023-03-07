
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#ifndef MEDYAN_Structure_SubSystemFunc_hpp
#define MEDYAN_Structure_SubSystemFunc_hpp

#include <type_traits>

#include "Structure/DynamicNeighbors.hpp"
#include "Structure/SubSystem.h"
#include "Structure/Movables.hpp"
#include "Structure/Reactables.hpp"
#include "Structure/Trackables.hpp"

namespace medyan {


struct SubSystemFunc {
    // Auxiliary function to get the corresponding element container.
    template< typename T, bool isConst >
    inline auto& getTrackableContainerAux(std::conditional_t<isConst, const SubSystem&, SubSystem&> sys) {
        if constexpr(std::is_same_v< T, Membrane >) return sys.membranes;
        if constexpr(std::is_same_v< T, Vertex >)   return sys.vertices;
        if constexpr(std::is_same_v< T, Edge >)     return sys.edges;
        if constexpr(std::is_same_v< T, Triangle >) return sys.triangles;
        if constexpr(std::is_same_v< T, MeshlessSpinVertex >) return sys.meshlessSpinVertices;
        if constexpr(std::is_same_v< T, FixedVertexAttachment >) return sys.fixedVertexAttachments;

        if constexpr(std::is_same_v< T, Bubble >)   return sys.bubbles;
        if constexpr(std::is_same_v< T, MTOC >)     return sys.mtocs;
        if constexpr(std::is_same_v< T, AFM >)      return sys.afms;
    }
    template< typename T >
    inline auto& getTrackableContainer(SubSystem&       sys) { return getTrackableContainerAux<T, false>(sys); }
    template< typename T >
    inline auto& getTrackableContainer(const SubSystem& sys) { return getTrackableContainerAux<T, true>(sys); }


    // Add a trackable to the system.
    // Returns the index of insertion.
    template< typename T, typename... Args >
    inline auto emplaceTrackable(SubSystem& sys, Args&&... args) {
        using namespace std;

        auto& container = getTrackableContainer<T>(sys);
        auto res = container.emplace();
        auto& obj = container[res];

        // Initialize trackable.
        trackable::initialize(*this, sys, obj, res, forward<Args>(args)...);

        // Initialize with traits.
        if constexpr(movable::isMovable<T>) { movable::updatePosition(sys, obj); }
        if constexpr(dynamic_neighbor::isDynamicNeighbor<T>) { dynamic_neighbor::addNeighbor(sys, obj); }

        return res;
    }


    // Remove a trackable from the system.
    template< typename T >
    inline void removeTrackable(SubSystem& sys, typename medyan::StableVector<T>::Index index) {
        using namespace std;

        auto& container = getTrackableContainer<T>(sys);
        auto& obj = container[index];

        // Finalize with traits.
        if constexpr(dynamic_neighbor::isDynamicNeighbor<T>) { dynamic_neighbor::removeNeighbor(sys, obj); }
        if constexpr(movable::isMovable<T>) { movable::unregister(sys, obj); }

        // Finalize trackable.
        trackable::finalize(*this, sys, obj, index);

        // Remove trackable.
        container.erase(index);
    }



    // Auxiliary function to decide whether to execute or ignore looping for a certain trackable.
    // Trait<T>::value is a bool indicating whether T satisfies the trait.
    template< typename T, template<typename> typename Trait, bool isConst, typename Func >
    inline void forEachTrackableTAux(std::conditional_t<isConst, const SubSystem&, SubSystem&> sys, Func&& func) {
        if constexpr(Trait<T>::value) {
            for(auto& x : getTrackableContainerAux<T, isConst>(sys)) {
                func(sys, x);
            }
        }
    }
    template< template<typename> typename Trait, bool isConst, typename Func >
    inline void forEachTrackableAux(std::conditional_t<isConst, const SubSystem&, SubSystem&> sys, Func&& func) {
        forEachTrackableTAux< Membrane, Trait, isConst >(sys, func);
        forEachTrackableTAux< Vertex,   Trait, isConst >(sys, func);
        forEachTrackableTAux< Edge,     Trait, isConst >(sys, func);
        forEachTrackableTAux< Triangle, Trait, isConst >(sys, func);
        forEachTrackableTAux< MeshlessSpinVertex, Trait, isConst >(sys, func);
        forEachTrackableTAux< FixedVertexAttachment, Trait, isConst >(sys, func);

        forEachTrackableTAux< Bubble,   Trait, isConst >(sys, func);
        forEachTrackableTAux< MTOC,     Trait, isConst >(sys, func);
        forEachTrackableTAux< AFM,      Trait, isConst >(sys, func);
    }
    template< template< typename > typename Trait, typename Func >
    inline void forEachTrackable(SubSystem& sys, Func&& func) {
        forEachTrackableAux< Trait, false >(sys, std::forward<Func>(func));
    }
    template< template< typename > typename Trait, typename Func >
    inline void forEachTrackable(const SubSystem& sys, Func&& func) {
        forEachTrackableAux< Trait, true >(sys, std::forward<Func>(func));
    }


    // Activate functions for certain types iteratively.
    template< typename Func >
    inline void forEachMovable(SubSystem& sys, Func&& func) {
        // Future: deprecate Movable class.
        for(auto m : Movable::getMovableList()) func(sys, *m);
        forEachTrackable< movable::IsMovable >(sys, func);
    }
    template< typename Func >
    inline void forEachReactable(SubSystem& sys, Func&& func) {
        // Future: deprecate Reactable class.
        for(auto r : Reactable::getReactableList()) func(sys, *r);
        forEachTrackable< reactable::IsReactable >(sys, func);
    }
};


} // namespace medyan

#endif

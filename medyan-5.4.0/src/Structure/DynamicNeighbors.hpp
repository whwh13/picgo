#ifndef MEDYAN_Structure_DynamicNeighbors_hpp
#define MEDYAN_Structure_DynamicNeighbors_hpp

#include <type_traits>

#include "Structure/Bubble.h"
#include "Structure/SurfaceMesh/Triangle.hpp"
#include "utility.h"

namespace medyan::dynamic_neighbor {

template< typename T >
constexpr bool isDynamicNeighbor =
    std::is_same_v< T, Triangle > ||
    std::is_same_v< T, Bubble >;
template< typename T >
struct IsDynamicNeighbor : std::integral_constant< bool, isDynamicNeighbor<T> > {};

constexpr Overload addNeighbor {
    [](auto&& sys, Triangle& t) {
        for(auto pnl : sys.getNeighborLists()) {
            pnl->addDynamicNeighbor(&t);
        }
    },
    [](auto&& sys, Bubble& b) {
        if(sys.opBoundaryBubbleNL.has_value()) {
            sys.opBoundaryBubbleNL->addDynamicNeighbor(sys, b.sysIndex);
        }
        if(sys.opBubbleBubbleNL.has_value()) {
            sys.opBubbleBubbleNL->addDynamicNeighbor(sys, b.sysIndex);
        }
        if(sys.opBubbleBeadNL.has_value()) {
            sys.opBubbleBeadNL->addDynamicNeighbor(sys, b.sysIndex);
        }
    },
};


constexpr Overload removeNeighbor {
    [](auto&& sys, Triangle& t) {
        for(auto pnl : sys.getNeighborLists()) {
            pnl->removeDynamicNeighbor(&t);
        }
    },
    [](auto&& sys, Bubble& b) {
        if(sys.opBoundaryBubbleNL.has_value()) {
            sys.opBoundaryBubbleNL->removeDynamicNeighbor(sys, b.sysIndex);
        }
        if(sys.opBubbleBubbleNL.has_value()) {
            sys.opBubbleBubbleNL->removeDynamicNeighbor(sys, b.sysIndex);
        }
        if(sys.opBubbleBeadNL.has_value()) {
            sys.opBubbleBeadNL->removeDynamicNeighbor(sys, b.sysIndex);
        }
    },
};

} // namespace medyan::dynamic_neighbor

#endif

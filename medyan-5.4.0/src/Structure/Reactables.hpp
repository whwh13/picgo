#ifndef MEDYAN_Structure_Reactables_hpp
#define MEDYAN_Structure_Reactables_hpp

#include "Structure/Reactable.h"
#include "Structure/Special/AFM.h"
#include "Structure/Special/MTOC.h"
#include "Structure/SurfaceMesh/FuncMembraneChem.hpp"

namespace medyan::reactable {

template< typename T >
constexpr bool isReactable =
    std::is_same_v< T, AFM > ||
    std::is_same_v< T, MTOC > ||
    std::is_same_v< T, Membrane > ||
    std::is_base_of_v< Reactable, T >;
template< typename T >
struct IsReactable : std::integral_constant< bool, isReactable<T> > {};

constexpr Overload updateReactionRates {
    // Future: replace with concrete types instead of "Reactable".
    [](auto&& sys, Reactable& r) {
        r.updateReactionRates();
    },
    [](auto&& sys, AFM& obj) {
        obj.updateReactionRates();
    },
    [](auto&& sys, MTOC& obj) {
        obj.updateReactionRates();
    },
    [](auto&& sys, Membrane& obj) {
        setReactionRates(sys, obj.getMesh());
    },
};

} // namespace medyan::reactable

#endif

#ifndef MEDYAN_Parameter_ReactionEmissionAbsorption_hpp
#define MEDYAN_Parameter_ReactionEmissionAbsorption_hpp

#include <string>

#include "common.h"
#include "Util/Parser/StringParser.hpp"

namespace medyan {

struct ReactionEmissionAbsorptionSetup {
    enum class DyRateType {
        none,
        force,
    };

    std::string speciesName1;
    std::string speciesName2;
    std::string speciesType2; // Must be "diffusing" or "bulk".
    floatingpoint emiRateConst = 0; // In unit of 1/s.
    floatingpoint absRateConst = 0; // In unit of nm^3/s.
    DyRateType  dyRateType = DyRateType::none;
};

template<>
struct StringSerializerTrait<ReactionEmissionAbsorptionSetup::DyRateType> {
    using Type = ReactionEmissionAbsorptionSetup::DyRateType;
    Type parse(std::string_view sv) const {
        if     (sv == "none")  return Type::none;
        else if(sv == "force") return Type::force;
        throw std::runtime_error("Invalid vertex column type");
    }
    std::string toString(Type type) const {
        switch(type) {
            case Type::none:      return "none";
            case Type::force:     return "force";
            default:              return "unknown";
        }
    }
};

} // namespace medyan

#endif

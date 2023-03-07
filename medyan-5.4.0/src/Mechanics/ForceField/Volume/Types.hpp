#ifndef MEDYAN_Mechanics_ForceField_Volume_Types_hpp
#define MEDYAN_Mechanics_ForceField_Volume_Types_hpp

#include "Util/Parser/StringParser.hpp"

namespace medyan {

enum class CylinderVolumeExclusionFFType {
    none,
    integral,
    monomer,
};

template<>
struct StringSerializerTrait<CylinderVolumeExclusionFFType> {
    auto parse(std::string_view sv) const {
        if(sv == "none") {
            return CylinderVolumeExclusionFFType::none;
        } else if (sv == "integral" || sv == "REPULSION") { // "REPULSION" is for backward compatibility.
            return CylinderVolumeExclusionFFType::integral;
        } else if (sv == "monomer") {
            return CylinderVolumeExclusionFFType::monomer;
        } else {
            throw std::invalid_argument("Invalid cylinder volume exclusion type.");
        }
    }

    std::string toString(CylinderVolumeExclusionFFType val) const {
        switch(val) {
            case CylinderVolumeExclusionFFType::none:     return "none";
            case CylinderVolumeExclusionFFType::integral: return "integral";
            case CylinderVolumeExclusionFFType::monomer:  return "monomer";
            default:                                      return "";
        }
    }
};

} // namespace medyan

#endif

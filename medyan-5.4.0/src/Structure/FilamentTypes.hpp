#ifndef MEDYAN_Structure_FilamentTypes_hpp
#define MEDYAN_Structure_FilamentTypes_hpp

#include <stdexcept>
#include <string>

#include "Util/Parser/StringParser.hpp"

namespace medyan {

enum class FilamentModel {
    beadCylinder, // Original MEDYAN model
    bezier,       // Bezier spline
    gek,          // Geodesic Extensible Kirchoff model.
    gc,           // Geodesic Cosserat model.
};

inline constexpr bool isGeodesic(FilamentModel model) {
    return model == FilamentModel::gek || model == FilamentModel::gc;
}

template<>
struct StringSerializerTrait<FilamentModel> {
    auto parse(std::string_view sv) const {
        if (sv == "bead-cylinder") {
            return FilamentModel::beadCylinder;
        } else if (sv == "bezier") {
            return FilamentModel::bezier;
        } else if (sv == "gek") {
            return FilamentModel::gek;
        } else if (sv == "gc") {
            return FilamentModel::gc;
        } else {
            throw std::invalid_argument("Invalid filament model.");
        }
    }

    std::string toString(FilamentModel val) const {
        switch(val) {
            case FilamentModel::beadCylinder: return "bead-cylinder";
            case FilamentModel::bezier:       return "bezier";
            case FilamentModel::gek:          return "gek";
            case FilamentModel::gc:           return "gc";
            default:                          return "error";
        }
    }
};

} // namespace medyan

#endif

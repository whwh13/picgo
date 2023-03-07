#ifndef MEDYAN_Util_Parser_StringParser_hpp
#define MEDYAN_Util_Parser_StringParser_hpp

#include <charconv>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>

#include "common.h"

namespace medyan {


// The trait for parsing text.
// The trait should contain
// - a function `parse` that takes a string_view and returns a parsed value.
// - a function `toString` that takes converts the value to a string.
template< typename T, typename = void > struct StringSerializerTrait;

// Integral or floating point types.
template< typename T >
struct StringSerializerTrait< T, std::enable_if_t< std::is_arithmetic_v<T> && !std::is_same_v<T, bool> > > {
    // Parse a string into a variable.
    T parse(std::string_view sv) const {
        T res;
#ifdef __cpp_lib_to_chars
        const auto [p, ec] = std::from_chars(sv.data(), sv.data() + sv.size(), res);
        if(ec != std::errc()) {
            throw std::runtime_error("Invalid argument.");
        }
#else
// GCC (as of version 8.4.0) does not support from_chars for floating point.
        if constexpr(std::is_floating_point_v<T>) {
            res = std::stod(std::string(sv));
        }
        else if constexpr(std::is_signed_v<T>) {
            res = std::stoll(std::string(sv));
        }
        else {
            res = std::stoull(std::string(sv));
        }
#endif

        return res;
    }

    // Convert a variable to a string.
    std::string toString(const T& val) const {
        return format("{}", val);
    }
};

// Bool.
template<>
struct StringSerializerTrait< bool > {
    bool parse(std::string_view sv) const {
        if(sv == "true" ) { return true; }
        if(sv == "false") { return false; }
        throw std::runtime_error("Invalid argument.");
    }
    std::string toString(bool val) const {
        return val ? "true" : "false";
    }
};

// String type.
template<>
struct StringSerializerTrait< std::string > {
    // Parse a string into a variable.
    std::string parse(std::string_view sv) const {
        return std::string(sv);
    }

    // Convert a variable to a string.
    std::string toString(const std::string& val) const {
        return val;
    }
};

// Path type.
template<>
struct StringSerializerTrait< std::filesystem::path > {
    // Parse a string into a variable.
    std::filesystem::path parse(std::string_view sv) const {
        return std::filesystem::path(sv);
    }

    // Convert a variable to a string.
    std::string toString(const std::filesystem::path& val) const {
        return val.string();
    }
};


// Actual serialization functions.
//---------------------------------

// Function to parse a string into a variable.
template< typename T, typename Trait = StringSerializerTrait< T > >
inline T parse(std::string_view sv, Trait&& trait = Trait{}) {
    return trait.parse(sv);
}

// Variable gets parsed in place.
template< typename T, typename Trait = StringSerializerTrait< T > >
inline void parse(T& var, std::string_view sv, Trait&& trait = Trait{}) {
    var = parse<T>(sv, std::forward<Trait>(trait));
}

// Function to convert a variable to a string.
template< typename T, typename Trait = StringSerializerTrait< T > >
inline auto toString(const T& val, Trait&& trait = Trait{}) {
    return trait.toString(val);
}

} // namespace medyan

#endif

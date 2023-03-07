#include <cmath>
#include <cstdint>
#include <numeric>

#include <catch2/catch.hpp>

#include "Util/Parser/StringParser.hpp"

namespace medyan {

TEST_CASE("Variable serialization from/to string", "[Parser]") {
    using namespace std;

    SECTION("Integers") {
        int var = 100;

        parse(var, "200");
        CHECK(var == 200);
        parse(var, "-300");
        CHECK(var == -300);

        CHECK(parse<int>("400") == 400);

        CHECK_THROWS(parse(var, "some-invalid-stuff"));

        const auto uint8zero = (std::uint8_t)0;
        const auto intmin    = numeric_limits<int>::min();
        const auto ulongmax  = numeric_limits<unsigned long>::max();
        CHECK(parse<std::uint8_t> (toString(uint8zero)) == uint8zero);
        CHECK(parse<int>          (toString(intmin))    == intmin);
        CHECK(parse<unsigned long>(toString(ulongmax))  == ulongmax);
    }

    SECTION("Boolean") {
        bool var = true;

        parse(var, "false");
        CHECK(!var);
        parse(var, "true");
        CHECK(var);

        CHECK(parse<bool>("true"));
        CHECK_THROWS(parse(var, "some-invalid-stuff"));

        CHECK(toString(true) == "true");
        CHECK(toString(false) == "false");
    }

    SECTION("Floating points") {
        float varf = 100;
        double vard = 200;

        parse(varf, "3e-6");
        CHECK(varf == Approx(3e-6));

        parse(vard, "-4e100");
        CHECK(vard == -4e100);

        CHECK(parse<double>("500.001") == Approx(500.001));

        CHECK_THROWS(parse(varf, "some-invalid-stuff"));

        const float floats[] {
            -numeric_limits<float>::max(),
            0.0f,
            numeric_limits<float>::min(),
            numeric_limits<float>::epsilon(),
            numeric_limits<float>::max(),
            numeric_limits<float>::infinity(),
            numeric_limits<float>::quiet_NaN(),
        };
        const double doubles[] {
            -numeric_limits<double>::max(),
            0.0,
            numeric_limits<double>::min(),
            numeric_limits<double>::epsilon(),
            numeric_limits<double>::max(),
            numeric_limits<double>::infinity(),
            numeric_limits<double>::quiet_NaN(),
        };
        for(auto x : floats) {
            if(isnan(x)) {
                CHECK(isnan(parse<float>(toString(x))));
            } else {
                CHECK(parse<float>(toString(x)) == x);
            }
        }
        for(auto x : doubles) {
            if(isnan(x)) {
                CHECK(isnan(parse<double>(toString(x))));
            } else {
                CHECK(parse<double>(toString(x)) == x);
            }
        }
    }

    SECTION("Strings") {
        std::string str;

        parse(str, "string-1");
        CHECK(str == "string-1");

        const char* strWithStuff = "   string with some random stuff !_\n\"/\\   ";
        CHECK(parse<std::string>(strWithStuff) == strWithStuff);

        CHECK(toString(std::string(strWithStuff)) == strWithStuff);
    }

    SECTION("File paths") {
        std::filesystem::path path;

        parse(path, "P:\\some\\path\\to\\file.txt");
        CHECK(path == "P:\\some\\path\\to\\file.txt");

        parse(path, "../path/to/file");
        CHECK(path == "../path/to/file");

        const char* pathStuff = "/a/b/c/d/e/f/g";
        CHECK(parse<std::filesystem::path>(pathStuff) == pathStuff);

        CHECK(toString(std::filesystem::path(pathStuff)) == pathStuff);
    }
}

} // namespace medyan

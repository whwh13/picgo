#include <catch2/catch.hpp>

#include "utility.h"

TEST_CASE("Scope guards", "[Utility]") {
    using namespace medyan;

    int x = 0;
    {
        ScopeGuard guard1{ [&]() { x = 1; } };
        ScopeGuard guard2{ [&]() { x = 2; } };
        REQUIRE(x == 0);

        {
            ScopeGuard guard3{ [&]() { x = 3; } };
            REQUIRE(x == 0);
        }
        REQUIRE(x == 3);

        // Anonymous (improper usage).
        ScopeGuard { [&]() { x = 4; } };
        REQUIRE(x == 4);
    }
    REQUIRE(x == 1);
}

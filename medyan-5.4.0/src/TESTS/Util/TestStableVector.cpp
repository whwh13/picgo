#include "catch2/catch.hpp"

#include "Util/StableVector.hpp"

TEST_CASE("Stable vector functions", "[Util][StableVector]") {
    using namespace std;
    using namespace medyan;

    struct Dummy {
        int x = 0;
        Dummy() = default;
        Dummy(int x) : x(x) {}
        Dummy(float halfX) : x(2 * halfX) {}
    };

    StableVector< Dummy > dvec;

    // Test adding and removing elements from the vectorized data.
    SECTION("insertion and removal") {
        REQUIRE(dvec.empty());

        {
            Dummy d3 {3};
            auto index = dvec.insert(d3);
            // d: [3]
            // deleted: []
            REQUIRE(index.value == 0);
            REQUIRE(dvec.at(index).x == 3);
            REQUIRE(dvec[index].x == 3);

            REQUIRE(!dvec.empty());
            REQUIRE(dvec.size() == 1);
            REQUIRE(dvec.value.size() == 1);
            REQUIRE(dvec.deletedIndices.size() == 0);
        }

        {
            auto index = dvec.insert(Dummy {5});
            // d: [3, 5]
            // deleted: []
            REQUIRE(index.value == 1);
            REQUIRE(dvec.at(index).x == 5);

            REQUIRE(dvec.size() == 2);
            REQUIRE(dvec.value.size() == 2);
            REQUIRE(dvec.deletedIndices.size() == 0);
        }

        {
            dvec.insert(Dummy {7});
            dvec.insert(Dummy {9});
            dvec.insert(Dummy {11});
            // d: [3, 5, 7, 9, 11]
            // deleted: []

            dvec.erase(3);
            dvec.emplace(6.6f); // int(2 * 6.6f) -> 13
            // d: [3, 5, 7, 13, 11]
            // deleted: []
            REQUIRE(dvec.at(3).x == 13);

            REQUIRE(dvec.size() == 5);
            REQUIRE(dvec.value.size() == 5);
            REQUIRE(dvec.deletedIndices.size() == 0);
        }

        {
            dvec.insert(Dummy {15});
            dvec.erase(3);
            dvec.erase(1);
            dvec.erase(2);
            // d: [3, -, -, -, 11, 15]
            // deleted: [3, 1, 2]
            REQUIRE(dvec.size() == 3);
            REQUIRE(dvec.value.size() == 6);
            REQUIRE(dvec.value[1].has_value() == false);
            REQUIRE(dvec.value[2].has_value() == false);
            REQUIRE(dvec.value[3].has_value() == false);
            REQUIRE_THROWS(dvec.at(2));
            REQUIRE_THROWS(dvec.at(6));
            REQUIRE(dvec.deletedIndices.size() == 3);
            REQUIRE(dvec.deletedIndices[0].value == 3);
            REQUIRE(dvec.deletedIndices[1].value == 1);
            REQUIRE(dvec.deletedIndices[2].value == 2);
        }

        {
            REQUIRE(dvec.begin()->x == 3);

            dvec.erase(dvec.begin());
            // d: [-, -, -, -, 11, 15]
            // deleted: [3, 1, 2, 0]
            REQUIRE(dvec.size() == 2);
            REQUIRE(dvec.value.size() == 6);
            REQUIRE(dvec.value[0].has_value() == false);
            REQUIRE(dvec.deletedIndices.size() == 4);
            REQUIRE(dvec.deletedIndices[3].value == 0);

            REQUIRE(dvec.begin()->x == 11);
        }
    }

    // Test iterators.
    SECTION("iterators") {
        dvec.insert(Dummy {3});
        dvec.insert(Dummy {5});
        dvec.insert(Dummy {7});
        dvec.insert(Dummy {9});
        dvec.insert(Dummy {11});
        dvec.insert(Dummy {13});
        dvec.erase(3);
        dvec.erase(0);
        dvec.erase(2);
        dvec.erase(5);
        // d: [-, 5, -, -, 11, -]
        // deleted: [3, 0, 2, 5]
        REQUIRE(dvec.size() == 2);

        REQUIRE(dvec.begin().index.value == 1);
        REQUIRE(dvec.begin()->x == 5);
        REQUIRE(dvec.end().index.value == 6);

        {
            auto it = dvec.begin();
            ++it;

            REQUIRE(it.index.value == 4);
            REQUIRE(dvec.indexat(it).value == 4);
            REQUIRE(it->x == 11);
            const auto it4 = it;

            ++it;
            REQUIRE(it.index.value == 6);
            REQUIRE(dvec.indexat(it).value == 6);
            REQUIRE(it == dvec.cend());

            --it;
            REQUIRE(it == it4);
            REQUIRE(it != dvec.begin());
            REQUIRE(it != dvec.end());
        }

        {
            auto it = dvec.begin();
            std::advance(it, 1);
            it->x = 17;
            // d: [-, 5, -, -, 17, -]
            // deleted: [3, 0, 2, 5]

            std::vector<Dummy> extracted(dvec.cbegin(), dvec.cend());
            REQUIRE(extracted.size() == 2);
            CHECK(extracted[0].x == 5);
            CHECK(extracted[1].x == 17);
        }

    }

    // Test that StableVector can be erased while iterating.
    SECTION("Erase while iterating") {
        // Let dummy.x = 10 * index.
        dvec.insert(Dummy {0});
        dvec.insert(Dummy {10});
        dvec.insert(Dummy {20});
        dvec.insert(Dummy {30});
        dvec.insert(Dummy {40});
        dvec.insert(Dummy {50});
        dvec.insert(Dummy {60});
        dvec.erase(5);
        dvec.erase(0);
        dvec.erase(1);

        REQUIRE(dvec.size() == 4);
        const int testIndexRemaining[] = {2, 3, 4, 6};
        int       testLoopIndex = 0;

        SECTION("Delete by index") {
            for(const auto& dummy : dvec) {
                REQUIRE(dummy.x == 10 * testIndexRemaining[testLoopIndex]);

                dvec.erase(testIndexRemaining[testLoopIndex]);
                CHECK(dvec.size() == 3 - testLoopIndex);

                ++testLoopIndex;
            }
        }
        SECTION("Delete by iterator") {
            for(auto it = dvec.begin(); it != dvec.end(); ++it) {
                REQUIRE(it->x == 10 * testIndexRemaining[testLoopIndex]);

                dvec.erase(it);
                CHECK(dvec.size() == 3 - testLoopIndex);

                ++testLoopIndex;
            }
        }

        REQUIRE(testLoopIndex == 4);
        REQUIRE(dvec.empty());
    }
}

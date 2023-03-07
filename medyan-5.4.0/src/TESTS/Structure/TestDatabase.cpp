#include <memory> // unique_ptr

#include "catch2/catch.hpp"

#include "Structure/Database.h"

namespace medyan {

TEST_CASE("Database tests", "[Database]") {
    SECTION("Database with dynamic indexing") {
        using DummyData = std::vector<int>;
        struct Dummy : Database< Dummy, false, DummyData > {
            Dummy(int x) : Database< Dummy, false, DummyData >(x) {}
        };

        // Test adding and removing elements from the vectorized data
        SECTION("vectorized data") {

            auto dummy0 = std::make_unique<Dummy>(0);
            // d: [0]
            REQUIRE(Dummy::getElements().size() == 1);
            REQUIRE(Dummy::getDbData().size() == 1);
            REQUIRE(dummy0->getIndex() == 0);

            auto dummy1 = std::make_unique<Dummy>(1);
            auto dummy2 = std::make_unique<Dummy>(2);
            // d: [0, 1, 2]
            REQUIRE(Dummy::getElements().size() == 3);
            REQUIRE(Dummy::getDbData().size() == 3);
            REQUIRE(dummy2->getIndex() == 2);
            REQUIRE(Dummy::getDbData()[2] == 2);

            dummy1.reset(nullptr);
            // d: [0, 2]
            REQUIRE(Dummy::getElements().size() == 2);
            REQUIRE(Dummy::getDbData().size() == 2);
            REQUIRE(dummy2->getIndex() == 1);
            REQUIRE(Dummy::getDbData()[1] == 2);

            auto dummy3 = std::make_unique<Dummy>(3);
            // d: [0, 2, 3]
            REQUIRE(Dummy::getElements().size() == 3);
            REQUIRE(Dummy::getDbData().size() == 3);
            REQUIRE(dummy3->getIndex() == 2);
            REQUIRE(Dummy::getDbData()[2] == 3);

            dummy3.reset(nullptr);
            // d: [0, 2]
            REQUIRE(Dummy::getElements().size() == 2);
            REQUIRE(Dummy::getDbData().size() == 2);
        }

        // Test copy/move constructors
        SECTION("copy/move constructors") {
            const auto dummy0PreviousVal = 100;
            auto dummy0 = std::make_unique<Dummy>(dummy0PreviousVal);
            const auto dummy0PreviousId = dummy0->getId();
            REQUIRE(dummy0->getIndex() == 0);

            auto dummy0Copy = std::make_unique<Dummy>(*dummy0);
            // id
            CHECK(dummy0->getId() == dummy0PreviousId);
            CHECK(dummy0->getId() != dummy0Copy->getId());
            // index
            CHECK(dummy0->getIndex() == 0);
            CHECK(dummy0Copy->getIndex() == 1);
            // reverse index
            REQUIRE(Dummy::numElements() == 2);
            CHECK(Dummy::getElements()[0] == dummy0.get());
            CHECK(Dummy::getElements()[1] == dummy0Copy.get());
            // value
            CHECK(Dummy::getDbData()[0] == dummy0PreviousVal);
            CHECK(Dummy::getDbData()[1] == dummy0PreviousVal);

            auto dummy0Move = std::make_unique<Dummy>(std::move(*dummy0));
            // id
            CHECK(dummy0Move->getId() == dummy0PreviousId);
            CHECK(dummy0->getId() != dummy0Move->getId());
            // index
            CHECK(dummy0->getIndex() == 2);
            CHECK(dummy0Move->getIndex() == 0);
            // reverse index
            REQUIRE(Dummy::numElements() == 3);
            CHECK(Dummy::getElements()[2] == dummy0.get());
            CHECK(Dummy::getElements()[0] == dummy0Move.get());
            // value
            CHECK(Dummy::getDbData()[0] == dummy0PreviousVal);
            CHECK(Dummy::getDbData()[1] == dummy0PreviousVal);
            CHECK(Dummy::getDbData()[2] == 0);
        }
    }

    SECTION("Database with stable indexing") {
        using DummyData = std::vector< int >;
        struct Dummy : Database< Dummy, true, DummyData > {
            Dummy(int x) : Database< Dummy, true, DummyData >(x) {}
        };

        // Test adding and removing elements from the vectorized data
        SECTION("vectorized data") {
            auto dummy0 = std::make_unique<Dummy>(0);
            // d: [0]
            // deleted: []
            REQUIRE(Dummy::getElements().size() == 1);
            REQUIRE(Dummy::getDbData().size() == 1);
            REQUIRE(dummy0->getIndex() == 0);
            REQUIRE(Dummy::getStableElement(0) == dummy0.get());
            REQUIRE(dummy0->getStableIndex() == 0);

            auto dummy1 = std::make_unique<Dummy>(1);
            auto dummy2 = std::make_unique<Dummy>(2);
            auto dummy3 = std::make_unique<Dummy>(3);
            auto dummy4 = std::make_unique<Dummy>(4);
            dummy3.reset(nullptr);
            auto dummy5 = std::make_unique<Dummy>(5);
            // dummies: 0, 1, 2, 4, 5
            // d: [0, 1, 2, 5, 4]
            // deleted: []
            REQUIRE(Dummy::getDbData().size() == 5);
            REQUIRE(dummy5->getIndex() == 4);
            REQUIRE(dummy5->getStableIndex() == 3);
            REQUIRE(Dummy::getStableElement(3) == dummy5.get());

            auto dummy6 = std::make_unique<Dummy>(6);
            dummy5.reset(nullptr);
            dummy1.reset(nullptr);
            dummy2.reset(nullptr);
            // dummies: 0, 6, 4
            // d: [0, 1, 2, 5, 4, 6]
            // deleted: [5, 1, 2]
            REQUIRE(Dummy::getElements().size() == 3);
            REQUIRE(Dummy::getDbData().size() == 6);

            REQUIRE(dummy4->getIndex() == 2);
            REQUIRE(dummy4->getStableIndex() == 4);
            REQUIRE(Dummy::getStableElement(4) == dummy4.get());

            REQUIRE(dummy6->getIndex() == 1);
            REQUIRE(dummy6->getStableIndex() == 5);
            REQUIRE(Dummy::getStableElement(5) == dummy6.get());

            Dummy::rearrange();
            // dummies: 0, 6, 4
            // d: [0, 4, 6]
            // deleted: []
            REQUIRE(Dummy::getElements().size() == 3);
            REQUIRE(Dummy::getDbData().size() == 3);

            REQUIRE(dummy4->getIndex() == 2);
            REQUIRE(Dummy::getStableElement(1) == dummy4.get());
            REQUIRE(dummy4->getStableIndex() == 1);

            REQUIRE(dummy6->getIndex() == 1);
            REQUIRE(Dummy::getStableElement(2) == dummy6.get());
            REQUIRE(dummy6->getStableIndex() == 2);
        }

        // Test copy/move constructors
        SECTION("copy/move constructors") {
            // Clear up the mess
            //-------------------------
            Dummy::rearrange();
            REQUIRE(Dummy::rawNumStableElements() == 0);

            // set up dummy 0
            //-------------------------
            const auto dummy0PreviousVal = 100;
            auto dummy0 = std::make_unique<Dummy>(dummy0PreviousVal);
            const auto dummy0PreviousId = dummy0->getId();
            REQUIRE(dummy0->getStableIndex() == 0);

            // copy/move with push_back
            //-------------------------
            auto dummy0Copy = std::make_unique<Dummy>(*dummy0);
            // stable index
            CHECK(dummy0->getStableIndex() == 0);
            CHECK(dummy0Copy->getStableIndex() == 1);
            // reverse stable index
            REQUIRE(Dummy::rawNumStableElements() == 2);
            CHECK(Dummy::getStableElement(0) == dummy0.get());
            CHECK(Dummy::getStableElement(1) == dummy0Copy.get());
            // value
            CHECK(Dummy::getDbData()[0] == dummy0PreviousVal);
            CHECK(Dummy::getDbData()[1] == dummy0PreviousVal);

            auto dummy0Move = std::make_unique<Dummy>(std::move(*dummy0));
            // stable index
            CHECK(dummy0->getStableIndex() == 2);
            CHECK(dummy0Move->getStableIndex() == 0);
            // reverse stable index
            REQUIRE(Dummy::rawNumStableElements() == 3);
            CHECK(Dummy::getStableElement(2) == dummy0.get());
            CHECK(Dummy::getStableElement(0) == dummy0Move.get());
            // value
            CHECK(Dummy::getDbData()[0] == dummy0PreviousVal);
            CHECK(Dummy::getDbData()[1] == dummy0PreviousVal);
            CHECK(Dummy::getDbData()[2] == 0);


            // intermediate steps
            //-------------------------
            // dummy 0 is now at stable index 2
            dummy0Copy.reset(nullptr);
            dummy0Move.reset(nullptr);
            REQUIRE(Dummy::rawNumStableElements() == 3);
            REQUIRE(Dummy::numElements() == 1);
            // Deleted indices: [1, 0]

            // set up dummy 1
            //-------------------------
            auto& dummy1 = dummy0;
            const auto dummy1PreviousVal = 200;
            Dummy::getDbData()[dummy1->getStableIndex()] = 200;

            // copy/move with hole filling
            //-------------------------
            auto dummy1Copy = std::make_unique<Dummy>(*dummy1);
            // Deleted indices: [1]
            // stable index
            CHECK(dummy1->getStableIndex() == 2);
            CHECK(dummy1Copy->getStableIndex() == 0);
            // reverse stable index
            REQUIRE(Dummy::rawNumStableElements() == 3);
            CHECK(Dummy::getStableElement(2) == dummy1.get());
            CHECK(Dummy::getStableElement(0) == dummy1Copy.get());
            // value
            CHECK(Dummy::getDbData()[0] == dummy1PreviousVal);
            CHECK(Dummy::getDbData()[2] == dummy1PreviousVal);

            auto dummy1Move = std::make_unique<Dummy>(std::move(*dummy1));
            // Deleted indices: []
            // stable index
            CHECK(dummy1->getStableIndex() == 1);
            CHECK(dummy1Move->getStableIndex() == 2);
            // reverse stable index
            REQUIRE(Dummy::rawNumStableElements() == 3);
            CHECK(Dummy::getStableElement(1) == dummy1.get());
            CHECK(Dummy::getStableElement(2) == dummy1Move.get());
            // value
            CHECK(Dummy::getDbData()[0] == dummy1PreviousVal);
            CHECK(Dummy::getDbData()[1] == 0);
            CHECK(Dummy::getDbData()[2] == dummy1PreviousVal);
        }
    }
}

} // namespace medyan

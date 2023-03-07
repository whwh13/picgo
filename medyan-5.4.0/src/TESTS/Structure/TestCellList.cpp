#include <memory> // unique_ptr
#include <utility> // move
#include <vector>

#include "catch2/catch.hpp"

#include "Rand.h"
#include "Structure/CellList.hpp"

TEST_CASE("Cell list tests", "[CellList]") {
    using namespace medyan;

    // Types
    //-------------------------------------------------------------------------
    struct DummyElement;
    struct DummyCell;
    using Manager = CellListManager< DummyElement*, DummyCell* >;

    struct DummyElement {
        CellListElementUser< DummyElement*, DummyCell* > ele;
    };
    struct DummyCell {
        CellListHeadUser< DummyElement*, DummyCell* > cell;
    };

    // Test objects
    //-------------------------------------------------------------------------
    Manager m;

    const auto registerDummyElement = [&](DummyElement& e) { e.ele .manager = &m; };
    const auto registerDummyCell    = [&](DummyCell   & c) { c.cell.manager = &m; };

    SECTION("Cell list operation and view") {
        DummyCell dc[3];
        for(auto& c : dc) { m.addHead(&c, c.cell); registerDummyCell(c); }
        REQUIRE(m.numHeads() == 3);

        DummyElement de[6];
        for(auto& e : de) registerDummyElement(e);

        const auto checkCellContent = [&](std::size_t dci, std::vector< std::size_t > content) {
            INFO("Checking cell " << dci << " content");

            const auto view = m.getElements(dc[dci].cell);
            std::vector< DummyElement* > deCell(view.begin(), view.end());

            REQUIRE(deCell.size() == content.size());
            const auto size = deCell.size();
            for(std::size_t i = 0; i < size; ++i)
                CHECK(deCell[i] == de + content[i]);
        };

        // Cell | 0   | 1   | 2   |
        //------|-----|-----|-----|
        // Ele  | 0 4 |     | 1   |
        m.addElement(de + 0, de[0].ele, dc[0].cell);
        m.addElement(de + 1, de[1].ele, dc[2].cell);
        m.addElement(de + 4, de[4].ele, dc[0].cell);

        checkCellContent(0, {0, 4});
        checkCellContent(1, {});
        checkCellContent(2, {1});

        // Cell | 0   | 1   | 2   |
        //------|-----|-----|-----|
        // Ele  | 4   |     | 1 2 |
        m.removeElement(de[0].ele);
        m.addElement(de + 5, de[5].ele, dc[2].cell);
        m.addElement(de + 2, de[2].ele, dc[2].cell);
        m.removeElement(de[5].ele);

        checkCellContent(0, {4});
        checkCellContent(1, {});
        checkCellContent(2, {1, 2});

        // Cell | 0   | 1   | 2   |
        //------|-----|-----|-----|
        // Ele  | 4   | 2 1 |     |
        m.updateElement(de[2].ele, dc[1].cell);
        m.updateElement(de[1].ele, dc[1].cell);

        checkCellContent(0, {4});
        checkCellContent(1, {2, 1});
        checkCellContent(2, {});

        // Check CellView
        {
            REQUIRE(m.getElements(dc[0].cell).size() == 1);
            REQUIRE(m.getElements(dc[2].cell).size() == 0);
            REQUIRE(m.getElements(dc[1].cell).empty() == false);
            REQUIRE(m.getElements(dc[2].cell).empty() == true);
        }

        // Check CellView iterators
        // Iterator dereference and comparison
        {
            auto it1 = m.getElements(dc[1].cell).begin();
            auto it1_cp = it1;
            ++it1;
            ++it1_cp;
            REQUIRE(*it1 == de + 1);
            REQUIRE(it1 == it1_cp);

            // End iterator
            ++it1;
            auto it1_end = m.getElements(dc[1].cell).end();
            REQUIRE(it1 == it1_end);
            --it1;
            --it1_end;
            REQUIRE(it1 == it1_end);

            // Empty cell
            REQUIRE(m.getElements(dc[2].cell).begin() == m.getElements(dc[2].cell).end());
        }

        // Check clearing.
        m.clearElements();
        CHECK(m.numElements() == 0);
        CHECK(m.numHeads() == 3);

        m.clearAll();
        CHECK(m.numElements() == 0);
        CHECK(m.numHeads() == 0);
    }

    SECTION("Random big number test") {
        std::vector< std::unique_ptr< DummyElement > > es;
        std::vector< std::unique_ptr< DummyCell    > > cs;

        const auto buildCells = [&](std::size_t num) {
            for(std::size_t i = 0; i < num; ++i) {
                cs.push_back(std::make_unique< DummyCell >());
                m.addHead(cs.back().get(), cs.back()->cell);
                registerDummyCell(*cs.back());
            }
        };

        const auto addElements = [&](std::size_t num) {
            std::uniform_int_distribution<> ci(0, cs.size() - 1); // cs.size() must > 0

            for(std::size_t i = 0; i < num; ++i) {
                es.push_back(std::make_unique< DummyElement >());
                m.addElement(es.back().get(), es.back()->ele, cs[ci(Rand::eng)]->cell);
                registerDummyElement(*es.back());
            }
        };
        const auto updateElements = [&](std::size_t num) {
            std::uniform_int_distribution<> ei(0, es.size() - 1); // es.size() must > 0
            std::uniform_int_distribution<> ci(0, cs.size() - 1); // cs.size() must > 0

            for(std::size_t i = 0; i < num; ++i) {
                m.updateElement(es[ei(Rand::eng)]->ele, cs[ci(Rand::eng)]->cell);
            }
        };
        const auto removeElements = [&](std::size_t num) {
            for(std::size_t i = 0; i < num; ++i) {
                std::uniform_int_distribution<> ei(0, es.size() - 1); // es.size() must > 0
                const auto idx = ei(Rand::eng);
                m.removeElement(es[idx]->ele);
                if(idx != es.size() - 1) {
                    es[idx] = std::move(es.back());
                }
                es.pop_back();
            }
        };

        const auto checkSizeConsistency = [&]() {
            std::size_t tot = 0;
            for(const auto& c : cs) tot += m.getElements(c->cell).size();
            REQUIRE(tot == es.size());
        };

        const auto numCells = 100;
        buildCells(numCells);

        addElements(5000);
        checkSizeConsistency();

        updateElements(2000);
        checkSizeConsistency();

        removeElements(4000);
        checkSizeConsistency();

        addElements(5000);
        checkSizeConsistency();

        updateElements(1000);
        checkSizeConsistency();

        removeElements(6000);
        checkSizeConsistency();
    }
}

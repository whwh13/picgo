
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "catch2/catch.hpp"

#include "Composite.h"

namespace medyan {
// Some concrete initializable types.
struct ConcreteComponent : Component {
    virtual void printSelf() const override {}
    virtual int getType() override { return 0; }
};
struct ConcreteComposite : Composite {
    virtual void printSelf() const override {}
    virtual int getType() override { return 0; }
};

TEST_CASE("Component test", "[Composite]") {
    ConcreteComponent cp;
    CHECK(nullptr == cp.getParent());
    CHECK(true == cp.isRoot());
    CHECK(0 == cp.numberOfChildren());
    CHECK(false == cp.isComposite());
    CHECK(0 == cp.countSpecies());
    CHECK(0 == cp.countReactions());
    CHECK("Component" == cp.getFullName());
}

TEST_CASE("Composite test", "[Composite]") {
    ConcreteComposite cp;
    CHECK(nullptr == cp.getParent());
    CHECK(true == cp.isRoot());
    CHECK(0 == cp.numberOfChildren());
    CHECK(true == cp.isComposite());
    CHECK(0 == cp.countSpecies());
    CHECK(0 == cp.countReactions());
    CHECK("Composite" == cp.getFullName());
    CHECK(0 == cp.children().size());
}

TEST_CASE("Composite children test", "[Composite]") {
    using namespace std;

    // Basic hierarchy. X is the root composite.
    auto X = make_unique<ConcreteComposite>();
    Component *Xptr = X.get();
    X->addChild(make_unique<ConcreteComponent>());
    Component *chA = X->children(0);
    X->addChild(make_unique<ConcreteComponent>());
    Component *chB = X->children(1);
    X->addChild(make_unique<ConcreteComponent>());
    Component *chC = X->children(2);
    CHECK(3 == X->numberOfChildren());
    X->removeChild(chA);
    CHECK(2 == X->numberOfChildren());
    CHECK(chB == X->children(0));
    CHECK(chC == X->children(1));
    
    CHECK(true == X->isRoot());
    CHECK(X.get() == chB->getParent());
    CHECK(X.get() == chC->getParent());

    // Multiple layers. X now has a parent Y.
    auto Y = make_unique<ConcreteComposite>();
    Y->addChild(move(X));
    Y->addChild(make_unique<ConcreteComponent>());
    Y->addChild(make_unique<ConcreteComponent>());
    CHECK(3 == Y->numberOfChildren());
    CHECK(5 == Y->countDescendents());
    
    CHECK(Y.get() == chC->getRoot());
    CHECK(Y.get() == Xptr->getRoot());
    CHECK(Xptr == chB->getParent());
    
    Y->removeChild(Xptr);
    CHECK(2 == Y->numberOfChildren());
    CHECK(2 == Y->countDescendents());

    // Child transfer.
    auto Z = make_unique<ConcreteComposite>();
    Z->addChild(make_unique<ConcreteComponent>());
    Z->addChild(make_unique<ConcreteComposite>());
    Z->addChild(make_unique<ConcreteComponent>());
    static_cast<ConcreteComposite*>(Z->children(1))->addChild(make_unique<ConcreteComponent>());
    REQUIRE(3 == Z->numberOfChildren());
    REQUIRE(4 == Z->countDescendents());

    Z->transferChild(Z->children(1), *Y);
    CHECK(2 == Z->numberOfChildren());
    CHECK(2 == Z->countDescendents());
    CHECK(3 == Y->numberOfChildren());
    CHECK(4 == Y->countDescendents());

    {
        ConcreteComponent irrelevantComponent;
        CHECK_THROWS_AS(Z->transferChild(&irrelevantComponent, *Y), std::out_of_range);
    }
}

} // namespace medyan


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

#include "Composite.h"

#include "Visitor.h"

namespace medyan {
bool Composite::apply (Visitor &v) {
    bool res_self = v.visit(this); //pre-order
    if(!res_self)
        return false;
    for (auto &c : children()) {
        bool res_child = c->apply(v);
        if(!res_child)
            return false;
    }
    return true;
}


bool Composite::apply (SpeciesVisitor &v) {
    bool res_self = apply_impl(v); //pre-order
    if(!res_self)
        return false;
    for (auto &c : children()) {
        bool res_child = c->apply(v);
        if(!res_child)
            return false;
    }
    return true;
}

bool Composite::apply (ReactionVisitor &v) {
    bool res_self = apply_impl(v); //pre-order
    if(!res_self)
        return false;
    for (auto &c : children()) {
        bool res_child = c->apply(v);
        if(!res_child)
            return false;
    }
    return true;
}

} // namespace medyan

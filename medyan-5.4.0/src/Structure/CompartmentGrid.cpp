
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

#include "CompartmentGrid.h"
#include "Chemistry/ChemSim.h"

#include "MathFunctions.h"
#include "SysParams.h"
#include "Controller/GController.h"

namespace medyan {
using namespace mathfunc;

void CompartmentGrid::addChemSimReactions(medyan::ChemSim* chem) {
    
    for(auto& c : compartmentList)
        c->addChemSimReactions(chem);
    
    for(auto &r : _bulkReactions.reactions())
        chem->addReaction(r.get());
    
}

species_copy_t CompartmentGrid::countDiffusingSpecies(const string& name) {
    
    species_copy_t copyNum = 0;

    for(auto &c : compartmentList) {
        
        auto s = c->findSpeciesByName(name);
        assert(s != nullptr && "Counting a diffusing species that does not exist.");
        
        copyNum += s->getN();
    }
    return copyNum;
}


species_copy_t CompartmentGrid::countBulkSpecies(const string& name) {
    
    auto s = findSpeciesBulkByName(name);
    assert(s != nullptr && "Counting a bulk species that does not exist.");
    
    return s->getN();
}

} // namespace medyan

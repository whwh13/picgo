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

#include "CMonomer.h"

#include "CBound.h"
#include "Structure/Compartment.h"
#include "SysParams.h"

namespace medyan {
CMonomer::CMonomer(const CMonomer& rhs, Compartment* c)

    : CMonomer(rhs._filamentType) {

    for(int i = 0; i < _numFSpecies[_filamentType]; i++) {
        
        //clone and add to array
        Species* s = rhs._speciesFilament[i];
        Species* sNew = s->clone();
        
        c->addSpeciesUnique(unique_ptr<Species>(sNew));
        _speciesFilament[i] = sNew;
    }
    
    //For bound species, transfer the CBound (if any)
    
    for(int i = 0; i < _numBSpecies[_filamentType]; i++) {
        
        //clone and add to array
        SpeciesBound* s = rhs._speciesBound[i];
        SpeciesBound* sNew = s->clone();
        
        c->addSpeciesUnique(unique_ptr<Species>(sNew));
        _speciesBound[i] = sNew;
        
        //update cbound
        CBound* cBound = s->getCBound();
        if(cBound != nullptr) {
            //set species
            if(cBound->getFirstSpecies() == s)
                cBound->setFirstSpecies(sNew);
            else
                cBound->setSecondSpecies(sNew);
        }
    }
}

//PRINT

void CMonomer::print()
{
    for(int i = 0; i < _numFSpecies[_filamentType]; i++) {
        Species* s = _speciesFilament[i];
        if(s != nullptr && areEqual(s->getN(), 1.0))
            cout << s->getName();
    }
    for(int i = 0; i < _numBSpecies[_filamentType]; i++) {
        SpeciesBound* s = _speciesBound[i];
        if(s != nullptr && areEqual(s->getN(), 1.0))
            cout << s->getName();
    }
}

//GETTERS

Species* CMonomer::speciesFilament(int index) {
    short offset = speciesFilamentIndex_[_filamentType][SPECIESFILAMENT].start;
    return _speciesFilament[index + offset];
}
Species* CMonomer::speciesPlusEnd (int index) {
    short offset = speciesFilamentIndex_[_filamentType][SPECIESPLUSEND].start;
    return _speciesFilament[index + offset];
}
Species* CMonomer::speciesMinusEnd(int index) {
    short offset = speciesFilamentIndex_[_filamentType][SPECIESMINUSEND].start;
    return _speciesFilament[index + offset];
}

SpeciesBound* CMonomer::speciesBound(int index) {
    short offset = speciesBoundIndex_[_filamentType][SPECIESBOUND].start;
    return _speciesBound[index + offset];
}
SpeciesBound* CMonomer::speciesLinker(int index) {
    short offset = speciesBoundIndex_[_filamentType][SPECIESLINKER].start;
    return _speciesBound[index + offset];
}
SpeciesBound* CMonomer::speciesMotor(int index) {
    short offset = speciesBoundIndex_[_filamentType][SPECIESMOTOR].start;
    return _speciesBound[index + offset];
}
SpeciesBound* CMonomer::speciesBrancher(int index) {
    short offset = speciesBoundIndex_[_filamentType][SPECIESBRANCHER].start;
    return _speciesBound[index + offset];
}


//GET ACTIVE

short CMonomer::activeSpeciesFilament() {
    short numFilamentSpecies = speciesFilamentIndex_[_filamentType][SPECIESFILAMENT].size;
    short offset             = speciesFilamentIndex_[_filamentType][SPECIESFILAMENT].start;
    
    for(int i = 0; i < numFilamentSpecies; i++) {
        Species* s = _speciesFilament[i + offset];
        if(s != nullptr && areEqual(s->getN(), 1.0)) return i;
    }
    return -1;
}
short CMonomer::activeSpeciesPlusEnd() {
    short numPlusEndSpecies = speciesFilamentIndex_[_filamentType][SPECIESPLUSEND].size;
    short offset            = speciesFilamentIndex_[_filamentType][SPECIESPLUSEND].start;
    
    for(int i = 0; i < numPlusEndSpecies; i++) {
        Species* s = _speciesFilament[i + offset];
        if(s != nullptr && areEqual(s->getN(), 1.0))
            return i;
    }
    return -1;
}
short CMonomer::activeSpeciesMinusEnd() {
    short numMinusEndSpecies = speciesFilamentIndex_[_filamentType][SPECIESMINUSEND].size;
    short offset             = speciesFilamentIndex_[_filamentType][SPECIESMINUSEND].start;
    
    for(int i = 0; i < numMinusEndSpecies; i++) {
        Species* s = _speciesFilament[i + offset];
        if(s != nullptr && areEqual(s->getN(), 1.0))
            return i;
    }
    return -1;
}

short CMonomer::activeSpeciesBrancher() {
    short numBrancherSpecies = speciesBoundIndex_[_filamentType][SPECIESBRANCHER].size;
    short offset             = speciesBoundIndex_[_filamentType][SPECIESBRANCHER].start;
    
    for(int i = 0; i < numBrancherSpecies; i++) {
        SpeciesBound* s = _speciesBound[i + offset];
        if(s != nullptr && areEqual(s->getN(), 1.0)) return i;
    }
    return -1;
}

bool CMonomer::isConsistent() {
    
    //check all species between 0 and 1 inclusive
    for(int i = 0; i < _numFSpecies[_filamentType]; i++) {
        
        if(!areEqual(_speciesFilament[i]->getN(), 1.0) &&
           !areEqual(_speciesFilament[i]->getN(), 0.0)) {
            
            cout << _speciesFilament[i]->getName() << " has an invalid copy number. It is = "
                 << _speciesFilament[i]->getN() << " and is at species index " << i << "." << endl;
            
            return false;
        }
    }
    //check filament species
    if(activeSpeciesFilament() != -1 &&
       (activeSpeciesPlusEnd() != -1 ||
        activeSpeciesMinusEnd() != -1)) {
           
        cout << "Has a simultaneously active filament and plus/minus end species." << endl;
           
        return false;
    }
    
    return true;
}

vector<short> CMonomer::_numFSpecies = vector<short>(MAX_FILAMENT_TYPES);
vector<short> CMonomer::_numBSpecies = vector<short>(MAX_FILAMENT_TYPES);

} // namespace medyan



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

#include "ChemManager.h"

#include <cmath>

#include "ChemCallbacks.h"
#include "CompartmentGrid.h"
#include "ReactionTemplate.h"
#include "BindingManager.h"

#include "SysParams.h"
#include "MathFunctions.h"
#include "Boundary.h"
#include "Structure/CompartmentChem.hpp"


namespace medyan {
using namespace mathfunc;

void ChemManager::setupBindingSites(ChemParams& chemParams, const medyan::SimulConfig& sc) {
    const auto& chemData = sc.chemistryData;

    //set binding indices
    //check if binding sites are valid and mark
    for(int filType = 0; filType < sc.chemParams.numFilaments; filType++) {


        if(chemData.B_BINDING_INDEX[filType] != "") {
            auto it = find(chemData.speciesBound[filType].begin(),
                           chemData.speciesBound[filType].end(),
                           chemData.B_BINDING_INDEX[filType]);

            if(it == chemData.speciesBound[filType].end()) {

                cout << "The brancher binding site listed is not a valid bound species. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else {
                chemParams.brancherBoundIndex[filType] = it - chemData.speciesBound[filType].begin();
            }
        }

        if(chemData.L_BINDING_INDEX[filType] != "") {

            auto it = find(chemData.speciesBound[filType].begin(),
                           chemData.speciesBound[filType].end(),
                           chemData.L_BINDING_INDEX[filType]);

            if(it == chemData.speciesBound[filType].end()) {

                cout << "The linker binding site listed is not a valid bound species. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else {
                chemParams.linkerBoundIndex[filType] = it - chemData.speciesBound[filType].begin();
            }
        }

        if(chemData.M_BINDING_INDEX[filType] != "") {

            auto it = find(chemData.speciesBound[filType].begin(),
                           chemData.speciesBound[filType].end(),
                           chemData.M_BINDING_INDEX[filType]);

            if(it == chemData.speciesBound[filType].end()) {

                cout << "The motor binding site listed is not a valid bound species. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
            else {
                chemParams.motorBoundIndex[filType] = it - chemData.speciesBound[filType].begin();
            }
        }

        //for initialization of cylinders
        chemParams.bindingIndices[filType].push_back(chemParams.brancherBoundIndex[filType]);

        if(chemParams.brancherBoundIndex[filType] != chemParams.linkerBoundIndex[filType])
            chemParams.bindingIndices[filType].push_back(chemParams.linkerBoundIndex[filType]);

        if(chemParams.brancherBoundIndex[filType] != chemParams.motorBoundIndex[filType] &&
           chemParams.linkerBoundIndex[filType]   != chemParams.motorBoundIndex[filType])
            chemParams.bindingIndices[filType].push_back(chemParams.motorBoundIndex[filType]);
    }
}

void ChemManager::configCMonomer(const medyan::SimulConfig& sc) {
    auto& chemData = sc.chemistryData;

    std::vector<std::vector< CMonomer::SpeciesSubarrayIndex >> speciesFilamentIndex(sc.chemParams.numFilaments);
    std::vector<std::vector< CMonomer::SpeciesSubarrayIndex >> speciesBoundIndex(sc.chemParams.numFilaments);
    for(int filType = 0; filType < sc.chemParams.numFilaments; filType++) {

        //set up static CMonomer things
        CMonomer::_numFSpecies[filType] = chemData.speciesFilament[filType].size() +
        chemData.speciesPlusEnd[filType].size()  +
        chemData.speciesMinusEnd[filType].size();

        CMonomer::_numBSpecies[filType] = chemData.speciesBound[filType].size()   +
        chemData.speciesLinker[filType].size()  +
        chemData.speciesMotor[filType].size()   +
        chemData.speciesBrancher[filType].size();

        //set up species offsets
        medyan::Index sf1 = chemData.speciesFilament[filType].size();
        medyan::Index sf2 = chemData.speciesPlusEnd[filType].size();
        medyan::Index sf3 = chemData.speciesMinusEnd[filType].size();

        medyan::Index sb1 = chemData.speciesBound[filType].size();
        medyan::Index sb2 = chemData.speciesLinker[filType].size();
        medyan::Index sb3 = chemData.speciesMotor[filType].size();
        medyan::Index sb4 = chemData.speciesBrancher[filType].size();

        //create offset vector for filament
        speciesFilamentIndex[filType] = {
            CMonomer::SpeciesSubarrayIndex { 0,         sf1 },
            CMonomer::SpeciesSubarrayIndex { sf1,       sf2 },
            CMonomer::SpeciesSubarrayIndex { sf1 + sf2, sf3 },
        };
        speciesBoundIndex[filType] = {
            CMonomer::SpeciesSubarrayIndex { 0,               sb1 },
            CMonomer::SpeciesSubarrayIndex { sb1,             sb2 },
            CMonomer::SpeciesSubarrayIndex { sb1 + sb2,       sb3 },
            CMonomer::SpeciesSubarrayIndex { sb1 + sb2 + sb3, sb4 },
        };
    }
    CMonomer::setSpeciesFilamentIndex(std::move(speciesFilamentIndex));
    CMonomer::setSpeciesBoundIndex(std::move(speciesBoundIndex));
}

void ChemManager::initCMonomer(CMonomer* m, short filamentType, Compartment* c, const ChemistryData& chemData) {

    // FILAMENT SPECIES

    int fIndex = 0;
    for(auto &f : chemData.speciesFilament[filamentType]) {
        Species* sf = c->addSpeciesUnique(
            std::make_unique< Species >(SpeciesNamesDB::genUniqueFilName(f), 0, 1, SpeciesType::FILAMENT, RSpeciesType::REG));
        m->_speciesFilament[fIndex] = sf;
        fIndex++;
    }
    for (auto &p : chemData.speciesPlusEnd[filamentType]) {
        Species* sp = c->addSpeciesUnique(
            std::make_unique< Species >(SpeciesNamesDB::genUniqueFilName(p), 0, 1, SpeciesType::PLUSEND, RSpeciesType::REG));
        m->_speciesFilament[fIndex] = sp;
        fIndex++;

    }
    for (auto &mi : chemData.speciesMinusEnd[filamentType]) {
        Species* smi = c->addSpeciesUnique(
            std::make_unique< Species >(SpeciesNamesDB::genUniqueFilName(mi), 0, 1, SpeciesType::MINUSEND, RSpeciesType::REG));
        m->_speciesFilament[fIndex] = smi;
        fIndex++;
    }

    // BOUND SPECIES

    int bIndex = 0;
    for (auto &b : chemData.speciesBound[filamentType]) {
        SpeciesBound* sb =
        c->addSpeciesBound(SpeciesNamesDB::genUniqueFilName(b));
        m->_speciesBound[bIndex] = sb;
        bIndex++;
    }
    for (auto &l : chemData.speciesLinker[filamentType]) {
        SpeciesBound* sl =
            c->addSpeciesBound(SpeciesNamesDB::genUniqueFilName(l), 0, 1, SpeciesType::LINKER);
        m->_speciesBound[bIndex] = sl;
        bIndex++;
    }
    for (auto &mo : chemData.speciesMotor[filamentType]) {
        SpeciesBound* sm =
            c->addSpeciesBound(SpeciesNamesDB::genUniqueFilName(mo), 0, 1, SpeciesType::MOTOR);
        m->_speciesBound[bIndex] = sm;
        bIndex++;
    }
    for (auto &br : chemData.speciesBrancher[filamentType]) {
        SpeciesBound* sbr =
            c->addSpeciesBound(SpeciesNamesDB::genUniqueFilName(br), 0, 1, SpeciesType::BRANCHER);
        m->_speciesBound[bIndex] = sbr;
        bIndex++;
    }
}

void ChemManager::genFilReactionTemplates(ChemParams& chemParams, const ChemistryData& chemData, DissipationTracker* pdt) {


    for(int filType = 0; filType < chemParams.numFilaments; filType++) {

        //set up reaction templates
        for(auto &r: chemData.polymerizationReactions[filType]) {

            vector<tuple<int, SpeciesType>> reactantTemplate;
            vector<tuple<int, SpeciesType>> productTemplate;
            FilamentReactionDirection d;

            //FIRST SPECIES MUST BE BULK OR DIFFUSING
            auto reactant = r.speciesReactantDiffusingBulk;
            if(reactant.find("BULK") != string::npos) {

                //Look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find_if(chemData.speciesBulk.begin(), chemData.speciesBulk.end(),
                                  [&name](auto&& element) {
                                      return element.name == name; });

                if(it == chemData.speciesBulk.end()) {
                    cout <<
                    "A bulk species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                reactantTemplate.push_back( tuple<int, SpeciesType>(SpeciesNamesDB::stringToInt(name), SpeciesType::BULK));
            }

            else if(reactant.find("DIFFUSING") != string::npos) {

                //Look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find_if(chemData.speciesDiffusing.begin(),chemData.speciesDiffusing.end(),
                                  [&name](auto&& element) {
                                      return element.name == name; });
                if(it == chemData.speciesDiffusing.end()) {
                    cout <<
                    "A diffusing species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                reactantTemplate.push_back(tuple<int, SpeciesType>(
                        SpeciesNamesDB::stringToInt(name), SpeciesType::DIFFUSING));
            }
            else {
                cout <<
                "First species listed in a polymerization reaction must be either bulk or diffusing. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }

            //SECOND SPECIES MUST BE PLUS OR MINUS END
            reactant = r.speciesReactantPlusEndMinusEnd;
            if(reactant.find("PLUSEND") != string::npos) {
                
                //look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find(chemData.speciesPlusEnd[filType].begin(), chemData.speciesPlusEnd[filType].end(), name);
                int position = 0;

                if(it != chemData.speciesPlusEnd[filType].end()) {

                    //get position of iterator
                    position = distance(chemData.speciesPlusEnd[filType].begin(), it);
                    reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::PLUSEND));

                    d = FilamentReactionDirection::FORWARD;
                }
                else {
                    cout <<
                    "A plus end species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }

            else if(reactant.find("MINUSEND") != string::npos) {

                //look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find(chemData.speciesMinusEnd[filType].begin(), chemData.speciesMinusEnd[filType].end(), name);
                int position = 0;

                if(it != chemData.speciesMinusEnd[filType].end()) {

                    //get position of iterator
                    position = distance(chemData.speciesMinusEnd[filType].begin(), it);
                    reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::MINUSEND));

                    d = FilamentReactionDirection::BACKWARD;
                }
                else {
                    cout <<
                    "A minus end species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }
            else {

                cout <<
                "Second species listed in a polymerization reaction must be either plusend or minusend. Exiting."
                << endl;
                exit(EXIT_FAILURE);

            }

            //FIRST PRODUCT SPECIES MUST BE FILAMENT SPECIES
            auto product = r.speciesProductFilament;
            if(product.find("FILAMENT") != string::npos) {

                //look up species, make sure in list
                string name = product.substr(0, product.find(":"));
                auto it = find(chemData.speciesFilament[filType].begin(), chemData.speciesFilament[filType].end(), name);
                int position = 0;

                if(it != chemData.speciesFilament[filType].end()) {

                    //get position of iterator
                    position = distance(chemData.speciesFilament[filType].begin(), it);
                    productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::FILAMENT));
                }
                else {
                    cout <<
                    "A filament species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }
            else {
                cout <<
                "Third species listed in a polymerization reaction must be filament. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }

            //SECOND PRODUCT SPECIES MUST BE PLUS OR MINUS END
            product = r.speciesProductPlusEndMinusEnd;
            //read strings, and look up type
            if(product.find("PLUSEND") != string::npos) {

                //look up species, make sure in list
                string name = product.substr(0, product.find(":"));
                auto it = find(chemData.speciesPlusEnd[filType].begin(), chemData.speciesPlusEnd[filType].end(), name);
                int position = 0;

                if(it != chemData.speciesPlusEnd[filType].end()) {

                    //get position of iterator
                    position = distance(chemData.speciesPlusEnd[filType].begin(), it);
                    productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::PLUSEND));
                }
                else {
                    cout <<
                    "A plusend species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }

            else if(product.find("MINUSEND") != string::npos) {

                //look up species, make sure in list
                string name = product.substr(0, product.find(":"));
                auto it = find(chemData.speciesMinusEnd[filType].begin(), chemData.speciesMinusEnd[filType].end(), name);
                int position = 0;

                if(it != chemData.speciesMinusEnd[filType].end()) {

                    //get position of iterator
                    position = distance(chemData.speciesMinusEnd[filType].begin(), it);
                    productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::MINUSEND));
                }
                else {
                    cout <<
                    "A plusend species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }
            else {
                cout <<
                "Fourth species listed in a polymerization reaction must be either plusend or minusend. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }

            chemParams.originalPolyPlusRate = r.rate;

            //Add polymerization managers
            if(d == FilamentReactionDirection::FORWARD)
                _filRxnTemplates[filType].emplace_back(
                new PolyPlusEndTemplate(filType, reactantTemplate, productTemplate, r.rate,
                        r.gnum, r.hrcdid, pdt));
            else
                _filRxnTemplates[filType].emplace_back(
                new PolyMinusEndTemplate(filType, reactantTemplate, productTemplate, r.rate,
                        r.gnum, r.hrcdid, pdt));
        }

        //set up reaction templates
        for(auto &r: chemData.depolymerizationReactions[filType]) {

            vector<tuple<int, SpeciesType>> reactantTemplate;
            vector<tuple<int, SpeciesType>> productTemplate;
            FilamentReactionDirection d;

            //FIRST REACTANT SPECIES MUST BE FILAMENT SPECIES
            auto reactant = r.speciesReactantFilament;
            if(reactant.find("FILAMENT") != string::npos) {

                //look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find(chemData.speciesFilament[filType].begin(), chemData.speciesFilament[filType].end(), name);
                int position = 0;

                if(it != chemData.speciesFilament[filType].end()) {

                    //get position of iterator
                    position = distance(chemData.speciesFilament[filType].begin(), it);
                    reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::FILAMENT));
                }
                else {
                    cout <<
                    "A filament species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }
            else {
                cout <<
                "First species listed in a depolymerization reaction must be filament. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }

            //SECOND REACTANT SPECIES MUST BE PLUS OR MINUS END
            reactant = r.speciesReactantPlusEndMinusEnd;
            //read strings, and look up type
            if(reactant.find("PLUSEND") != string::npos) {

                //look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find(chemData.speciesPlusEnd[filType].begin(), chemData.speciesPlusEnd[filType].end(), name);
                int position = 0;

                if(it != chemData.speciesPlusEnd[filType].end()) {

                    //get position of iterator
                    position = distance(chemData.speciesPlusEnd[filType].begin(), it);
                    reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::PLUSEND));

                    d = FilamentReactionDirection::BACKWARD;
                }
                else {
                    cout <<
                    "A plusend species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }

            else if(reactant.find("MINUSEND") != string::npos) {

                //look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find(chemData.speciesMinusEnd[filType].begin(), chemData.speciesMinusEnd[filType].end(), name);
                int position = 0;

                if(it != chemData.speciesMinusEnd[filType].end()) {

                    //get position of iterator
                    position = distance(chemData.speciesMinusEnd[filType].begin(), it);
                    reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::MINUSEND));
                    d = FilamentReactionDirection::FORWARD;
                }
                else {
                    cout <<
                    "A minusend species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }
            else {
                cout <<
                "Second species listed in a depolymerization reaction must be either plusend or minusend. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }


            //FIRST PRODUCT SPECIES MUST BE BULK OR DIFFUSING
            auto product = r.speciesProductDiffusingBulk;
            if(product.find("BULK") != string::npos) {

                //Look up species, make sure in list
                string name = product.substr(0, product.find(":"));
                auto it = find_if(chemData.speciesBulk.begin(), chemData.speciesBulk.end(),
                                  [&name](auto&& element) {
                                      return element.name == name; });

                if(it == chemData.speciesBulk.end()) {
                    cout <<
                    "A bulk species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                productTemplate.push_back(tuple<int, SpeciesType>(SpeciesNamesDB::stringToInt(name), SpeciesType::BULK));
            }

            else if(product.find("DIFFUSING") != string::npos) {

                //Look up species, make sure in list
                string name = product.substr(0, product.find(":"));
                auto it = find_if(chemData.speciesDiffusing.begin(),chemData.speciesDiffusing.end(),
                                  [&name](auto&& element) {
                                      return element.name == name; });
                if(it == chemData.speciesDiffusing.end()) {
                    cout <<
                    "A diffusing species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                productTemplate.push_back(tuple<int, SpeciesType>(
                SpeciesNamesDB::stringToInt(name), SpeciesType::DIFFUSING));
            }
            else {
                cout <<
                "Third species listed in a depolymerization reaction must be either bulk or diffusing. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }

            //SECOND PRODUCT SPECIES MUST BE PLUS OR MINUS END
            product = r.speciesProductPlusEndMinusEnd;
            if(product.find("PLUSEND") != string::npos) {

                //look up species, make sure in list
                string name = product.substr(0, product.find(":"));
                auto it = find(chemData.speciesPlusEnd[filType].begin(), chemData.speciesPlusEnd[filType].end(), name);
                int position = 0;

                if(it != chemData.speciesPlusEnd[filType].end()) {

                    //get position of iterator
                    position = distance(chemData.speciesPlusEnd[filType].begin(), it);
                    productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::PLUSEND));
                }
                else {
                    cout <<
                    "A plus end species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }

            else if(product.find("MINUSEND") != string::npos) {

                //look up species, make sure in list
                string name = product.substr(0, product.find(":"));
                auto it = find(chemData.speciesMinusEnd[filType].begin(), chemData.speciesMinusEnd[filType].end(), name);
                int position = 0;

                if(it != chemData.speciesMinusEnd[filType].end()) {

                    //get position of iterator
                    position = distance(chemData.speciesMinusEnd[filType].begin(), it);
                    productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::MINUSEND));
                }
                else {
                    cout <<
                    "A minus end species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }
            else {
                cout <<
                "Fourth species listed in a depolymerization reaction must be either plusend or minusend. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }

            //Add depolymerization managers
            if(d == FilamentReactionDirection::FORWARD)
                _filRxnTemplates[filType].emplace_back(
                new DepolyMinusEndTemplate(filType, reactantTemplate, productTemplate, r.rate, r.gnum, r.hrcdid, pdt));
            else
                _filRxnTemplates[filType].emplace_back(
                new DepolyPlusEndTemplate(filType, reactantTemplate, productTemplate, r.rate, r.gnum, r.hrcdid, pdt));
        }

        for(auto &r: chemData.motorWalkingReactions[filType]) {

            vector<tuple<int, SpeciesType>> reactantTemplate;
            vector<tuple<int, SpeciesType>> productTemplate;

            //read strings, and look up type
            ReactionType type;
            string species1;

            //FIRST REACTANT SPECIES MUST BE MOTOR
            auto reactant = r.speciesReactantMotor;
            if(reactant.find("MOTOR") != string::npos) {

                //look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find(chemData.speciesMotor[filType].begin(), chemData.speciesMotor[filType].end(), name);
                int position = 0;

                if(it != chemData.speciesMotor[filType].end()) {
                    species1 = name;

                    //check if forward or backward walking
                    if(reactant.find("N+1") != string::npos)
                        type = ReactionType::MOTORWALKINGBACKWARD;
                    else
                        type = ReactionType::MOTORWALKINGFORWARD;

                    //get position of iterator
                    position = distance(chemData.speciesMotor[filType].begin(), it);
                    reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::MOTOR));
                }
                else {
                    cout <<
                    "A motor species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }
            else{
                cout <<
                "First species listed in a motor walking reaction must be motor. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }


            //SECOND REACTANT SPECIES MUST BE EMPTY SITE
            reactant = r.speciesReactantEmptySite;
            //read strings, and look up type
            if(reactant.find("BOUND") != string::npos) {

                //look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find(chemData.speciesBound[filType].begin(), chemData.speciesBound[filType].end(), name);
                int position = 0;

                if(it != chemData.speciesBound[filType].end()) {

                    //get position of iterator
                    position = distance(chemData.speciesBound[filType].begin(), it);

                    if(position != chemParams.motorBoundIndex[filType]) {
                        cout <<
                        "Second species listed in a motor walking reaction must be the corresponding motor empty site. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }

                    reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::BOUND));
                }
                else {
                    cout <<
                    "A bound species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }
            else{
                cout <<
                "Second species listed in a motor walking reaction must be bound. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }

            //FIRST PRODUCT SPECIES MUST BE MOTOR
            auto product = r.speciesProductMotor;
            if(product.find("MOTOR") != string::npos) {

                //look up species, make sure in list
                string name = product.substr(0, product.find(":"));
                auto it = find(chemData.speciesMotor[filType].begin(), chemData.speciesMotor[filType].end(), name);
                int position = 0;

                if(name != species1) {
                    cout <<
                    "Motor species in reactants and products of motor walking reaction must be same. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }

                if((product.find("N+1") != string::npos && type != ReactionType::MOTORWALKINGFORWARD)) {

                    cout <<
                    "Motor walking reaction must have a direction (Check N and N+1 distinctions). Exiting."
                    <<endl;
                    exit(EXIT_FAILURE);

                }

                if(it != chemData.speciesMotor[filType].end()) {

                    //get position of iterator
                    position = distance(chemData.speciesMotor[filType].begin(), it);
                    productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::MOTOR));
                }
                else {
                    cout <<
                    "A motor species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }
            else {
                cout <<
                "Third species listed in a motor walking reaction must be motor. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }

            //SECOND PRODUCT SPECIES MUST BE EMPTY SITE
            product = r.speciesProductEmptySite;
            if(product.find("BOUND") != string::npos) {

                //look up species, make sure in list
                string name = product.substr(0, product.find(":"));
                auto it = find(chemData.speciesBound[filType].begin(), chemData.speciesBound[filType].end(), name);
                int position = 0;

                if(it != chemData.speciesBound[filType].end()) {

                    //get position of iterator
                    position = distance(chemData.speciesBound[filType].begin(), it);

                    if(position != chemParams.motorBoundIndex[filType]) {
                        cout <<
                        "Second species listed in a motor walking reaction must be the corresponding motor empty site. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }

                    productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::BOUND));
                }
                else {
                    cout <<
                    "A bound species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }
            else{
                cout <<
                "Fourth species listed in a motor walking reaction must be bound. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }

            //add reaction
            if(type == ReactionType::MOTORWALKINGFORWARD) {

                _filRxnTemplates[filType].emplace_back(
                new MotorWalkPTemplate(filType, reactantTemplate, productTemplate, r.rate, r.gnum, r.hrcdid, pdt));
            } else {
                _filRxnTemplates[filType].emplace_back(
                new MotorWalkMTemplate(filType, reactantTemplate, productTemplate, r.rate, -r.gnum, r.hrcdid+"m", pdt));
            }
        }

        //set up reaction templates
        for(auto &r: chemData.agingReactions[filType]) {

            vector<tuple<int, SpeciesType>> reactantTemplate;
            vector<tuple<int, SpeciesType>> productTemplate;

            //FIRST REACTANT SPECIES MUST BE FILAMENT, PLUS OR MINUS END
            auto reactant = r.speciesReactant;
            //read strings, and look up type
            if(reactant.find("FILAMENT") != string::npos) {

                //look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find(chemData.speciesFilament[filType].begin(), chemData.speciesFilament[filType].end(), name);
                int position = 0;

                if(it != chemData.speciesFilament[filType].end()) {

                    //get position of iterator
                    position = distance(chemData.speciesFilament[filType].begin(), it);
                    reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::FILAMENT));
                }
                else {
                    cout <<
                    "A filament species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }
            else if(reactant.find("PLUSEND") != string::npos) {

                //look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find(chemData.speciesPlusEnd[filType].begin(), chemData.speciesPlusEnd[filType].end(), name);
                int position = 0;

                if(it != chemData.speciesPlusEnd[filType].end()) {

                    //get position of iterator
                    position = distance(chemData.speciesPlusEnd[filType].begin(), it);
                    reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::PLUSEND));
                }
                else {
                    cout <<
                    "A plus end species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }
            else if(reactant.find("MINUSEND") != string::npos) {

                //look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find(chemData.speciesMinusEnd[filType].begin(), chemData.speciesMinusEnd[filType].end(), name);
                int position = 0;

                if(it != chemData.speciesMinusEnd[filType].end()) {

                    //get position of iterator
                    position = distance(chemData.speciesMinusEnd[filType].begin(), it);
                    reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::MINUSEND));
                }
                else {
                    cout <<
                    "A minus end species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }
            else{
                cout <<
                "First species listed in an aging reaction must be filament. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }

            //FIRST PRODUCT SPECIES MUST BE FILAMENT, PLUS, OR MINUS END
            auto product = r.speciesProduct;
            //read strings, and look up type
            if(product.find("FILAMENT") != string::npos) {

                //look up species, make sure in list
                string name = product.substr(0, product.find(":"));
                auto it = find(chemData.speciesFilament[filType].begin(), chemData.speciesFilament[filType].end(), name);
                int position = 0;

                if(it != chemData.speciesFilament[filType].end()) {

                    //get position of iterator
                    position = distance(chemData.speciesFilament[filType].begin(), it);
                    productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::FILAMENT));
                }
                else {
                    cout <<
                    "A filament species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }
            else if(product.find("PLUSEND") != string::npos) {

                //look up species, make sure in list
                string name = product.substr(0, product.find(":"));
                auto it = find(chemData.speciesPlusEnd[filType].begin(), chemData.speciesPlusEnd[filType].end(), name);
                int position = 0;

                if(it != chemData.speciesPlusEnd[filType].end()) {

                    //get position of iterator
                    position = distance(chemData.speciesPlusEnd[filType].begin(), it);
                    productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::PLUSEND));
                }
                else {
                    cout <<
                    "A plus end species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }
            else if(product.find("MINUSEND") != string::npos) {

                //look up species, make sure in list
                string name = product.substr(0, product.find(":"));
                auto it = find(chemData.speciesMinusEnd[filType].begin(), chemData.speciesMinusEnd[filType].end(), name);
                int position = 0;

                if(it != chemData.speciesMinusEnd[filType].end()) {

                    //get position of iterator
                    position = distance(chemData.speciesMinusEnd[filType].begin(), it);
                    productTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::MINUSEND));
                }
                else {
                    cout <<
                    "A minus end species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }
            else{
                cout <<
                "Second species listed in an aging reaction must be bound. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }

            //add reaction
            _filRxnTemplates[filType].emplace_back(new AgingTemplate(filType, reactantTemplate, productTemplate, r.rate, r.gnum, r.hrcdid, pdt));
        }


        //set up reaction templates
        for(auto &r: chemData.destructionReactions[filType]) {

            vector<tuple<int, SpeciesType>> reactantTemplate;
            vector<tuple<int, SpeciesType>> productTemplate;

            //FIRST SPECIES MUST BE PLUS END
            auto reactant = r.speciesReactantPlusEnd;
            if(reactant.find("PLUSEND") != string::npos) {

                //look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find(chemData.speciesPlusEnd[filType].begin(), chemData.speciesPlusEnd[filType].end(), name);

                if(it != chemData.speciesPlusEnd[filType].end()) {
                    //get position of iterator
                    int position = distance(chemData.speciesPlusEnd[filType].begin(), it);
                    reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::PLUSEND));
                }
                else {
                    cout <<
                    "A plus end species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }
            else {
                cout <<
                "First species listed in a destruction reaction must be plus end. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }

            //SECOND SPECIES MUST BE MINUS END
            reactant = r.speciesReactantMinusEnd;
            if(reactant.find("MINUSEND") != string::npos) {

                //look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find(chemData.speciesMinusEnd[filType].begin(), chemData.speciesMinusEnd[filType].end(), name);

                if(it != chemData.speciesMinusEnd[filType].end()) {
                    //get position of iterator
                    int position = distance(chemData.speciesMinusEnd[filType].begin(), it);
                    reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::MINUSEND));
                }
                else {
                    cout <<
                    "A minus end species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }
            else {
                cout <<
                "Second species listed in a destruction reaction must be minus end. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }

            ///ALL PRODUCTS MUST BE BULK OR DIFFUSING
            const auto addProduct = [&](const std::string& product) {

                if(product.find("BULK") != string::npos) {

                    //Look up species, make sure in list
                    string name = product.substr(0, product.find(":"));
                    auto it = find_if(chemData.speciesBulk.begin(), chemData.speciesBulk.end(),
                                      [&name](auto&& element) {
                                          return element.name == name; });

                    if(it == chemData.speciesBulk.end()) {
                        cout <<
                        "A bulk species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                    productTemplate.push_back(tuple<int, SpeciesType>(SpeciesNamesDB::stringToInt(name), SpeciesType::BULK));
                }

                else if(product.find("DIFFUSING") != string::npos) {

                    //Look up species, make sure in list
                    string name = product.substr(0, product.find(":"));
                    auto it = find_if(chemData.speciesDiffusing.begin(),chemData.speciesDiffusing.end(),
                                      [&name](auto&& element) {
                                          return element.name == name; });
                    if(it == chemData.speciesDiffusing.end()) {
                        cout <<
                        "A diffusing species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                    productTemplate.push_back(tuple<int, SpeciesType>(
                    SpeciesNamesDB::stringToInt(name), SpeciesType::DIFFUSING));
                }
                else {
                    cout <<
                    "Third species listed in a destruction reaction must be either bulk or diffusing. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            };
            addProduct(r.speciesProduct1);
            addProduct(r.speciesProduct2);

            //add reaction
            _filRxnTemplates[filType].emplace_back(new DestructionTemplate(filType, reactantTemplate, productTemplate, r.rate));
        }

        //set up reaction templates
        for(auto &r: chemData.severingReactions[filType]) {

            vector<tuple<int, SpeciesType>> reactantTemplate;
            vector<tuple<int, SpeciesType>> productTemplate;

            string reactant = get<0>(r);
            //read strings, and look up type


            // SPECIES MUST BE FILAMENT
            if(reactant.find("FILAMENT") != string::npos) {

                //look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find(chemData.speciesFilament[filType].begin(), chemData.speciesFilament[filType].end(), name);
                int position = 0;

                if(it != chemData.speciesFilament[filType].end()) {
                    //get position of iterator
                    position = distance(chemData.speciesFilament[filType].begin(), it);
                    reactantTemplate.push_back(tuple<int, SpeciesType>(position, SpeciesType::FILAMENT));
                }
                else {
                    cout <<
                    "A filament species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }
            else {
                cout <<
                "Reactant species listed in a severing reaction must be filament. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }

            //add reaction
            _filRxnTemplates[filType].emplace_back(new SeveringTemplate(filType, reactantTemplate, productTemplate, get<1>(r)));
        }
    }
}

void ChemManager::genFilBindingReactions(SubSystem& sys, const medyan::SimulConfig& sc, DissipationTracker* pdt) {
    using namespace std;

    auto grid = sys.getCompartmentGrid();
    auto& chemData = sc.chemistryData;

    //init subsystem ptr
    FilamentBindingManager::_subSystem = &sys;
    bool status = false;
	#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
	//If linker and motor reactions exist, create HybridBindingSearchManager
    short totalreactions = chemData.linkerReactions.size() + chemData.motorReactions.size();
	if(totalreactions)
		status = true;
	for (auto& C : grid->getCompartments()) {
		auto *Hbsn = new medyan::HybridBindingSearchManager(C.get());
		C->addHybridBindingSearchManager(Hbsn);
	}
	#endif


    for(int filType = 0; filType < sc.chemParams.numFilaments; filType++) {

        //loop through all compartments
        for(auto& C : grid->getCompartments()) {

            int managerIndex = 0;

            for(auto &r: chemData.branchingReactions[filType]) {

/*                cout<<"Considering compartment "<<C->getId()<<" coords "<<C->coordinates
                ()[0]<<" "<<C->coordinates()[1]<<" "<<C->coordinates()[2]<<" volFrac "<<
                C->getVolumeFrac()<<endl;*/
                
                //filament creation is not allowed in partially activated compartments
                //that volume fraction < threshold, be careful with teh threshold
//                if(C->getVolumeFrac() < 0.5) continue;

                vector<Species*> reactantSpecies;
                vector<Species*> productSpecies;

                vector<string> reactants = get<0>(r);
                vector<string> products = get<1>(r);

                //Checks on number of reactants, products
                if(reactants.size() != BRANCHINGREACTANTS ||
                   products.size() != BRANCHINGPRODUCTS) {
                    cout << "Invalid branching reaction. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                string brancherName;

                //FIRST PRODUCT SPECIES MUST BE BRANCHER
                short brancherInt;

                auto product = products[0];
                if(product.find("BRANCHER") != string::npos) {

                    //look up species, make sure in list
                    string name = product.substr(0, product.find(":"));
                    auto it = find(chemData.speciesBrancher[filType].begin(), chemData.speciesBrancher[filType].end(), name);

                    if(it != chemData.speciesBrancher[filType].end()) {

                        brancherName = name;

                        //get position of iterator
                        brancherInt = distance(chemData.speciesBrancher[filType].begin(), it);
                    }
                    else {
                        cout <<
                        "A brancher species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                }
                else{
                    cout <<
                    "Fourth species listed in a branching reaction must be brancher. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }

                //SECOND PRODUCT SPECIES MUST BE PLUS END
                short plusEnd;

                product = products[1];
                if(product.find("PLUSEND") != string::npos) {

                    //look up species, make sure in list
                    string name = product.substr(0, product.find(":"));
                    auto it = find(chemData.speciesPlusEnd[filType].begin(), chemData.speciesPlusEnd[filType].end(), name);

                    if(it != chemData.speciesPlusEnd[filType].end()) {

                        //get position of iterator
                        plusEnd = distance(chemData.speciesPlusEnd[filType].begin(), it);
                    }
                    else {
                        cout <<
                        "A plus end species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                }
                else {
                    cout <<
                    "Second product species listed in a branching reaction must be plus end. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }

                int numDiffusingReactant = 0; // Used in determining volume dependence

                //FIRST AND SECOND REACTANTS MUST BE BULK OR DIFFUSING
                auto reactant = reactants[0];

                if(reactant.find("BULK") != string::npos) {

                    //Look up species, make sure in list
                    string name = reactant.substr(0, reactant.find(":"));

                    auto it = find_if(chemData.speciesBulk.begin(), chemData.speciesBulk.end(),
                                      [&name](auto&& element) {
                                          return element.name == name; });

                    if(it == chemData.speciesBulk.end()) {
                        cout <<
                        "A bulk species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                    reactantSpecies.push_back(grid->findSpeciesBulkByName(name));
                }

                else if(reactant.find("DIFFUSING") != string::npos) {

                    //Look up species, make sure in list
                    string name = reactant.substr(0, reactant.find(":"));

                    auto it = find_if(chemData.speciesDiffusing.begin(),chemData.speciesDiffusing.end(),
                                      [&name](auto&& element) {
                                          return element.name == name; });
                    if(it == chemData.speciesDiffusing.end()) {
                        cout <<
                        "A diffusing species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                    reactantSpecies.push_back(C->findSpeciesByName(name));

                    ++numDiffusingReactant;
                }
                else {
                    cout <<
                    "First species listed in a branching reaction must be either bulk or diffusing. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }

                reactant = reactants[1];

                if(reactant.find("BULK") != string::npos) {

                    //Look up species, make sure in list
                    string name = reactant.substr(0, reactant.find(":"));
                    auto it = find_if(chemData.speciesBulk.begin(), chemData.speciesBulk.end(),
                                      [&name](auto&& element) {
                                          return element.name == name; });

                    if(it == chemData.speciesBulk.end()) {
                        cout <<
                        "A bulk species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                    reactantSpecies.push_back(grid->findSpeciesBulkByName(name));
                }

                else if(reactant.find("DIFFUSING") != string::npos) {

                    //Look up species, make sure in list
                    string name = reactant.substr(0, reactant.find(":"));
                    auto it = find_if(chemData.speciesDiffusing.begin(),chemData.speciesDiffusing.end(),
                                      [&name](auto&& element) {
                                          return element.name == name; });
                    if(it == chemData.speciesDiffusing.end()) {
                        cout <<
                        "A diffusing species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                    reactantSpecies.push_back(C->findSpeciesByName(name));

                    ++numDiffusingReactant;
                }
                else {
                    cout <<
                    "Second species listed in a branching reaction must be either bulk or diffusing. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }

                //THIRD REACTANT SPECIES MUST BE BOUND
                reactant = reactants[2];
                if(reactant.find("BOUND") != string::npos) {

                    //look up species, make sure in list
                    string name = reactant.substr(0, reactant.find(":"));
                    auto it = find(chemData.speciesBound[filType].begin(), chemData.speciesBound[filType].end(), name);
                    int position = 0;

                    if(it != chemData.speciesBound[filType].end()) {

                        //get position of iterator
                        position = distance(chemData.speciesBound[filType].begin(), it);

                        if(position != sc.chemParams.brancherBoundIndex[filType]) {
                            cout <<
                            "Third species listed in a branching reaction must be the corresponding brancher empty site. Exiting."
                            << endl;
                            exit(EXIT_FAILURE);
                        }

                        //find the species single binding, push
                        string bename = SpeciesNamesDB::genBindingName(brancherName, name);

                        reactantSpecies.push_back(C->findSpeciesByName(bename));
                    }
                    else {
                        cout <<
                        "A bound species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                }
                else{
                    cout <<
                    "Third species listed in a branching reaction must be bound. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }

                //Create reaction
                float onRate = get<2>(r);
                float offRate = get<3>(r);
                auto temp=SysParams::BUBBareRate;
                if(temp.size()>0)
                    temp[brancherInt]=offRate;
                else
                    temp.push_back(offRate);
                SysParams::BUBBareRate=temp;
                //get nucleation zone
                string nzstr = get<4>(r);
                NucleationZoneType nucleationZone;
                if(nzstr == "ALL")
                    nucleationZone = NucleationZoneType::ALL;
                else if(nzstr == "BOUNDARY")
                    nucleationZone = NucleationZoneType::BOUNDARY;
                else if(nzstr == "TOPBOUNDARY")
                    nucleationZone = NucleationZoneType::TOPBOUNDARY;
                else if(nzstr == "SIDEBOUNDARY")
                    nucleationZone = NucleationZoneType::SIDEBOUNDARY;
                else if (nzstr == "RIGHTBOUNDARY")
                    nucleationZone = NucleationZoneType::RIGHTBOUNDARY;
                else if(nzstr == "MEMBRANE")
                    nucleationZone = NucleationZoneType::MEMBRANE;
                else {
                    cout << "Nucleation zone type specified in a branching reaction not valid. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                floatingpoint nucleationDist = get<5>(r);

                ReactionBase* rxn = new Reaction<3,0>(reactantSpecies, onRate, false, C->getVolumeFrac(), -numDiffusingReactant);
                rxn->setReactionType(ReactionType::BRANCHING);

                C->addInternalReaction(rxn);

                vector<short> filTypevec = {short(filType), short(filType)};
                //create manager
                BranchingManager* bManager = new BranchingManager(rxn, C.get(), brancherInt, filTypevec, nucleationZone, nucleationDist);
                C->addFilamentBindingManager(bManager);

                #if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
                C->addBranchingBindingManager(bManager);
				#endif

                bManager->setMIndex(managerIndex++);


                //attach callback
                BranchingCallback bcallback(bManager, plusEnd, onRate, offRate, &sys);
                rxn->connect(bcallback);
            }
        }
    }

    for(auto& C : grid->getCompartments()) {

        int linkerReactionIndex = 0;
        int motorReactionIndex = 0;

        // Linker reactions.
        for(int ri = 0; ri < chemData.linkerReactions.size(); ++ri) {
            // The reaction index ri also defines the type of an actual linker.
            auto& r = chemData.linkerReactions[ri];

            vector<Species*> reactantSpecies;

            string linkerName1;
            string boundName1;
            int filType1 = 0;
            int linIndex1 = 0;
            string linkerName2;
            string boundName2;
            int filType2 = 0;
            int linIndex2 = 0;

            //FIRST TWO SPECIES IN PRODUCTS MUST BE LINKER
            {
                auto& product = r.productInfo.linkerBound1;

                if(product.find("LINKER") != string::npos) {

                    //look up species, make sure in list
                    linkerName1 = product.substr(0, product.find(":"));
                    tie(filType1, linIndex1) = locateSpecies(chemData.speciesLinker, linkerName1, "linker");
                }
                else {
                    cout <<
                    "Fourth species listed in a linker reaction must be linker. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }
            {
                auto& product = r.productInfo.linkerBound2;

                if(product.find("LINKER") != string::npos) {

                    //look up species, make sure in list
                    linkerName2 = product.substr(0, product.find(":"));
                    tie(filType2, linIndex2) = locateSpecies(chemData.speciesLinker, linkerName2, "linker");
                }
                else {
                    cout <<
                    "Fifth species listed in a linker reaction must be linker. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }

            //FIRST TWO REACTANTS SHOULD BE BOUND
            {
                auto& reactant = r.reactantInfo.speciesBound1;
                if(reactant.find("BOUND") != string::npos) {

                    //look up species, make sure in list
                    boundName1 = reactant.substr(0, reactant.find(":"));
                    auto [filType, boundIndex] = locateSpecies(chemData.speciesBound, boundName1, "bound");

                    if(filType != filType1) {
                        LOG(ERROR) << "Filament type of bound species " << boundName1
                            << " (" << filType << ") "
                            << " does not match linker species " << linkerName1
                            << " (" << filType1 << ").";
                        throw runtime_error("Inconsistent linker species.");
                    }
                    if(boundIndex != sc.chemParams.linkerBoundIndex[filType]) {
                        cout <<
                            "First species listed in a linker reaction must be the corresponding linker empty site. Exiting."
                            << endl;
                        throw logic_error("Incorrect bound index for linker.");
                    }
                }
                else {
                    cout <<
                    "First species listed in a linker reaction must be bound. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }
            {
                auto& reactant = r.reactantInfo.speciesBound2;
                if(reactant.find("BOUND") != string::npos) {

                    //look up species, make sure in list
                    boundName2 = reactant.substr(0, reactant.find(":"));
                    auto [filType, boundIndex] = locateSpecies(chemData.speciesBound, boundName2, "bound");

                    if(filType != filType2) {
                        LOG(ERROR) << "Filament type of bound species " << boundName2
                            << " (" << filType << ") "
                            << " does not match linker species " << linkerName2
                            << " (" << filType2 << ").";
                        throw runtime_error("Inconsistent linker species.");
                    }
                    if(boundIndex != sc.chemParams.linkerBoundIndex[filType]) {
                        cout <<
                            "First species listed in a linker reaction must be the corresponding linker empty site. Exiting."
                            << endl;
                        throw logic_error("Incorrect bound index for linker.");
                    }
                }
                else {
                    cout <<
                    "Second species listed in a linker reaction must be bound. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }


            // Get the name of the species pair binding.
            auto linkerBindingName = SpeciesNamesDB::genBindingName(linkerName1, boundName1, linkerName2, boundName2);
            reactantSpecies.push_back(C->findSpeciesByName(linkerBindingName));


            int numDiffusingReactant = 0; // Used in determining volume dependence

            //THIRD REACTANT SPECIES SHOULD BE BULK OR DIFFUSING
            {
                auto& reactant = r.reactantInfo.linkerDiffusing;
                if(reactant.find("BULK") != string::npos) {

                    //Look up species, make sure in list
                    string name = reactant.substr(0, reactant.find(":"));
                    auto it = find_if(chemData.speciesBulk.begin(), chemData.speciesBulk.end(),
                                      [&name](auto&& element) {
                                          return element.name == name; });

                    if(it == chemData.speciesBulk.end()) {
                        cout <<
                        "A bulk species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                    reactantSpecies.push_back(grid->findSpeciesBulkByName(name));
                }

                else if(reactant.find("DIFFUSING") != string::npos) {

                    //Look up species, make sure in list
                    string name = reactant.substr(0, reactant.find(":"));
                    auto it = find_if(chemData.speciesDiffusing.begin(),chemData.speciesDiffusing.end(),
                                      [&name](auto&& element) {
                                          return element.name == name; });
                    if(it == chemData.speciesDiffusing.end()) {
                        cout <<
                        "A diffusing species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                    reactantSpecies.push_back(C->findSpeciesByName(name));

                    ++numDiffusingReactant;
                }
                else {
                    cout << "Third species listed in a linker reaction must be bulk or diffusing. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }

            // Create the binding reaction.
            auto rxnUnique = make_unique<Reaction<2,0>>(reactantSpecies, r.onRate, false, C->getVolumeFrac(), -numDiffusingReactant);
            auto rxn = rxnUnique.get();
            C->addInternalReaction(move(rxnUnique));
            rxn->setReactionType(ReactionType::LINKERBINDING);

            // Dissipation
            if(sc.chemParams.dissTracking){
                rxn->setGNumber(r.gnum);
                rxn->setHRCDID(r.hrcdid);
            }


            // Create manager, using the reaction's index as the linker type.
            LinkerBindingManager* lManager = new LinkerBindingManager(rxn, C.get(), ri, {(short)filType1, (short)filType2}, linIndex1, linIndex2, r.rMax, r.rMin);
            C->addFilamentBindingManager(lManager);

            lManager->setNLIndex(linkerReactionIndex++);
            lManager->setMIndex(C->getFilamentBindingManagers().size() - 1);

            //attach callback
            LinkerBindingCallback lcallback(lManager, r.onRate, r.offRate, &sys, pdt);
            rxn->connect(lcallback);
#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
            auto Hbsm = C->getHybridBindingSearchManager();
            //1 refers to linker binding manager
            Hbsm->setbindingsearchparameter(lManager, 1, filType1, filType2, r.rMax, r.rMin);
#endif
        }

        // Motor reactions.
        for(int ri = 0; ri < chemData.motorReactions.size(); ++ri) {
            // The reaction index ri also defines the type of an actual motor.
            auto& r = chemData.motorReactions[ri];

            vector<Species*> reactantSpecies;

            string linkerName1;
            string boundName1;
            int filType1 = 0;
            int linIndex1 = 0;
            string linkerName2;
            string boundName2;
            int filType2 = 0;
            int linIndex2 = 0;

            //FIRST TWO SPECIES IN PRODUCTS MUST BE LINKER
            {
                auto& product = r.productInfo.linkerBound1;

                if(product.find("MOTOR") != string::npos) {

                    //look up species, make sure in list
                    linkerName1 = product.substr(0, product.find(":"));
                    tie(filType1, linIndex1) = locateSpecies(chemData.speciesMotor, linkerName1, "motor");
                }
                else {
                    cout <<
                    "Fourth species listed in a motor reaction must be motor. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }
            {
                auto& product = r.productInfo.linkerBound2;
                if(product.find("MOTOR") != string::npos) {

                    //look up species, make sure in list
                    linkerName2 = product.substr(0, product.find(":"));
                    tie(filType2, linIndex2) = locateSpecies(chemData.speciesMotor, linkerName2, "motor");
                }
                else {
                    cout <<
                    "Fifth species listed in a motor reaction must be motor. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }

            //FIRST TWO REACTANTS SHOULD BE BOUND
            {
                auto& reactant = r.reactantInfo.speciesBound1;
                if(reactant.find("BOUND") != string::npos) {

                    //look up species, make sure in list
                    boundName1 = reactant.substr(0, reactant.find(":"));
                    auto [filType, boundIndex] = locateSpecies(chemData.speciesBound, boundName1, "bound");

                    if(filType != filType1) {
                        LOG(ERROR) << "Filament type of bound species " << boundName1
                            << " (" << filType << ") "
                            << " does not match motor species " << linkerName1
                            << " (" << filType1 << ").";
                        throw runtime_error("Inconsistent motor species.");
                    }
                    if(boundIndex != sc.chemParams.motorBoundIndex[filType]) {
                        cout <<
                            "First species listed in a motor reaction must be the corresponding motor empty site. Exiting."
                            << endl;
                        throw logic_error("Incorrect bound index for motor.");
                    }
                }
                else {
                    cout <<
                    "First species listed in a motor reaction must be bound. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }
            {
                auto& reactant = r.reactantInfo.speciesBound2;
                if(reactant.find("BOUND") != string::npos) {

                    //look up species, make sure in list
                    boundName2 = reactant.substr(0, reactant.find(":"));
                    auto [filType, boundIndex] = locateSpecies(chemData.speciesBound, boundName2, "bound");

                    if(filType != filType2) {
                        LOG(ERROR) << "Filament type of bound species " << boundName2
                            << " (" << filType << ") "
                            << " does not match motor species " << linkerName2
                            << " (" << filType2 << ").";
                        throw runtime_error("Inconsistent motor species.");
                    }
                    if(boundIndex != sc.chemParams.motorBoundIndex[filType]) {
                        cout <<
                            "First species listed in a motor reaction must be the corresponding motor empty site. Exiting."
                            << endl;
                        throw logic_error("Incorrect bound index for motor.");
                    }
                }
                else {
                    cout <<
                    "Second species listed in a motor reaction must be bound. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }

            // Get the name of the species pair binding.
            auto linkerBindingName = SpeciesNamesDB::genBindingName(linkerName1, boundName1, linkerName2, boundName2);
            reactantSpecies.push_back(C->findSpeciesByName(linkerBindingName));


            int numDiffusingReactant = 0; // Used in determining volume dependence

            //THIRD REACTANT SPECIES SHOULD BE BULK OR DIFFUSING
            {
                auto& reactant = r.reactantInfo.linkerDiffusing;
                if(reactant.find("BULK") != string::npos) {

                    //Look up species, make sure in list
                    string name = reactant.substr(0, reactant.find(":"));
                    auto it = find_if(chemData.speciesBulk.begin(), chemData.speciesBulk.end(),
                                      [&name](auto&& element) {
                                          return element.name == name; });

                    if(it == chemData.speciesBulk.end()) {
                        cout <<
                        "A bulk species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                    reactantSpecies.push_back(grid->findSpeciesBulkByName(name));
                }

                else if(reactant.find("DIFFUSING") != string::npos) {

                    //Look up species, make sure in list
                    string name = reactant.substr(0, reactant.find(":"));
                    auto it = find_if(chemData.speciesDiffusing.begin(),chemData.speciesDiffusing.end(),
                                      [&name](auto&& element) {
                                          return element.name == name; });
                    if(it == chemData.speciesDiffusing.end()) {
                        cout <<
                        "A diffusing species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                    reactantSpecies.push_back(C->findSpeciesByName(name));

                    ++numDiffusingReactant;
                }
                else {
                    cout << "Third species listed in a motor reaction must be bulk or diffusing. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
            }

            //multiply by num heads to get rate
            ///CHANGED
            floatingpoint nh1 = sc.chemParams.motorNumHeadsMin[ri];
            floatingpoint nh2 = sc.chemParams.motorNumHeadsMax[ri];


            auto rxnUnique = make_unique<Reaction<2,0>>(reactantSpecies, r.onRate * (nh1 + nh2) / 2.0, false, C->getVolumeFrac(), -numDiffusingReactant);
            auto rxn = rxnUnique.get();
            C->addInternalReaction(move(rxnUnique));
            rxn->setReactionType(ReactionType::MOTORBINDING);

            // Dissipation
            if(sc.chemParams.dissTracking){
                rxn->setGNumber(r.gnum);
                rxn->setHRCDID(r.hrcdid);
            }


            // Create manager, using the reaction's index as the motor type.
            MotorBindingManager* mManager = new MotorBindingManager(rxn, C.get(), ri, {(short)filType1, (short)filType2}, linIndex1, linIndex2, r.rMax, r.rMin);
            C->addFilamentBindingManager(mManager);

            mManager->setNLIndex(motorReactionIndex++);
            mManager->setMIndex(C->getFilamentBindingManagers().size() - 1);

            //attach callback
            MotorBindingCallback mcallback(mManager, r.onRate, r.offRate, &sys);
            rxn->connect(mcallback);
#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
            auto Hbsm = C->getHybridBindingSearchManager();
            //2 let's it identify with a motor binding manager
            Hbsm->setbindingsearchparameter(mManager, 2, filType1, filType2, r.rMax, r.rMin);
#endif
        }

    } // Loop through Compartment

    for(int filType = 0; filType < sc.chemParams.numFilaments; filType++) {
        //init neighbor lists
        //if NOT DEFINED
#if !defined(HYBRID_NLSTENCILLIST) || !defined(SIMDBINDINGSEARCH)
        //get a compartment
        Compartment* C0 = grid->getCompartments()[0].get();
        for(auto &manager : C0->getFilamentBindingManagers()) {

            LinkerBindingManager* lManager;
            MotorBindingManager* mManager;

            if((lManager = dynamic_cast<LinkerBindingManager*>(manager.get()))) {

                auto nl =
                        new CylinderCylinderNL(lManager->getRMax() + sc.geoParams.cylinderSize[filType],
                                               0.0, true);

                //add to subsystem and manager
                LinkerBindingManager::_neighborLists.push_back(nl);
                sys.addNeighborList(nl);
            }

            else if((mManager = dynamic_cast<MotorBindingManager*>(manager.get()))) {

                auto nl =
                        new CylinderCylinderNL(mManager->getRMax() + sc.geoParams.cylinderSize[filType],
                                               0.0, true);

                //add to subsystem and manager
                MotorBindingManager::_neighborLists.push_back(nl);
                sys.addNeighborList(nl);
//#ifdef CUDAACCL_NL
//                mManager->assigncudavars();
//#endif
            }
        }
#endif
    } //Loop through Filament types
    //@@@ Defn:
    //1. BindingManager ->defined for each rxn
    //2. HybridBindingSearchManager -> defined for each compartment
    //3. HybridNeighborList -> defined for each unique distance pair.

    //@@@Current information
    //HybridNeighborList<=======HybridBindingSearchManager=====>BindingManager

    //@@@Missing connections
    //1. BindingManager=====X_WHICH_ID_TO_USE?_X=>>>>HybridNeighborList
    //2. BindingManager=====X_WHICH_ID_TO_USE?_X=>>>>HybridBindingSearchManager

    //At this point, HybridBindingSearchManager has access to individual BindingManagers
    // and an instance of HybridNeighborList.
    //A. HybridNeighborList has no information on distances for which
    // Neighborlists should be computed. Once that is done, HybridBindingSearchManager needs
    // to know which neighborlist ID within HybridNeighborlist to access.
    //B. Each BindingManager needs to have access to both HybridNeighborList ID and also
    // needs to know it's position in the BindingManager matrix.

    #if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
    Compartment *C0 = grid->getCompartments()[0].get();
    //status checks if there are linker and motor binding reactions for all filament Types
    if (status) {
        medyan::HybridBindingSearchManager::_HneighborList = sys.getHNeighborList();
        auto Hmanager = C0->getHybridBindingSearchManager();
        //A is accomplished with a single call to addtoHNeighborList
        Hmanager->addtoHNeighborList();
        //B is accomplished by looping through all HybridBindingSeachManagers and copying
        // back necessary information.
        for (auto& C : grid->getCompartments()) {
            C->getHybridBindingSearchManager()->copyInfotoBindingManagers();
        }
    }
    sys.initializeHNeighborList();
    #endif
}

void ChemManager::genSpecies(CompartmentGrid& grid, Compartment& protoCompartment, const medyan::SimulConfig& sc) {
    using namespace std;

    auto& chemData = sc.chemistryData;
    auto& chemParams = sc.chemParams;

    // add diffusing species (zero copy number for now)
    for(auto &sd : chemData.speciesDiffusing) {

        auto name = sd.name;
        auto diffCoeff = sd.diffusionCoefficient;
        auto rtypeStr = sd.rspeciesType;
        auto numEvents = sd.numEvents;

        RSpeciesType type;
        string rsptype(rtypeStr);

        if(rsptype == "AVG")
            type = RSpeciesType::AVG;
        else if (rsptype == "REG")
            type = RSpeciesType::REG;

        Species* s = protoCompartment.addSpeciesUnique(
            std::make_unique< Species >(name, 0, max_ulim, SpeciesType::DIFFUSING, type));

        //set num events if averaging
        if(rsptype == "AVG")
            ((RSpeciesAvg*)&s->getRSpecies())->setNumEvents(numEvents);

        protoCompartment.setDiffusionCoefficient(name, diffCoeff);
    }

    // add bulk species (zero copy number for now)
    for(auto &sb : chemData.speciesBulk) {

        auto& name = sb.name;
        auto& rtypeStr = sb.rspeciesType;

        RSpeciesType type;

        if(rtypeStr == "CONST")
            type = RSpeciesType::CONST;
        else if (rtypeStr == "REG")
            type = RSpeciesType::REG;

        grid.addSpeciesBulk(name, 0, max_ulim, SpeciesType::BULK, type);
    }

    for(int filType = 0; filType < sc.chemParams.numFilaments; filType++) {

        // create single binding and pair binding species
        for(auto &sb : chemData.speciesBrancher[filType]) {

            //look at brancher reaction that is associated with this species
            for(auto &rb : chemData.branchingReactions[filType]) {

                auto reactants = get<0>(rb);
                auto products = get<1>(rb);

                auto sb_bound = products[0].substr(0, products[0].find(":"));

                //                std::cout << reactants.size() << " " << products.size() << endl;

                //basic check because we have not yet checked reactions
                if(reactants.size() != BRANCHINGREACTANTS ||
                   products.size() != BRANCHINGPRODUCTS) {
                    cout << "Invalid branching reaction. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }

                if(sb_bound == sb) {
                    //look at bound species associated
                    string bound = reactants[2].substr(0, reactants[2].find(":"));

                    auto it = find(chemData.speciesBound[filType].begin(), chemData.speciesBound[filType].end(), bound);

                    //quick check for validity
                    if(it == chemData.speciesBound[filType].end()) {
                        cout <<
                        "A bound species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }

                    //add a single binding species with name sb + bound
                    protoCompartment.addSpeciesUnique(std::make_unique< Species >(
                        SpeciesNamesDB::genBindingName(sb, bound),
                        0, max_ulim,
                        SpeciesType::singleBinding, RSpeciesType::REG
                    ));
                }
            }
        }
    } // End loop filament type.


    // Add all pair bindings for linker reactions.
    for(auto &r : chemData.linkerReactions) {

        auto slinker1 = r.productInfo.linkerBound1.substr(0, r.productInfo.linkerBound1.find(':'));
        auto slinker2 = r.productInfo.linkerBound2.substr(0, r.productInfo.linkerBound2.find(':'));
        auto sbound1 = r.reactantInfo.speciesBound1.substr(0, r.reactantInfo.speciesBound1.find(':'));
        auto sbound2 = r.reactantInfo.speciesBound2.substr(0, r.reactantInfo.speciesBound2.find(':'));

        protoCompartment.addSpeciesUnique(std::make_unique< Species >(
            SpeciesNamesDB::genBindingName(slinker1, sbound1, slinker2, sbound2),
            0, max_ulim,
            SpeciesType::pairBinding, RSpeciesType::REG
        ));
    }

    // Add all pair bindings for motor reactions.
    for(auto &r : chemData.motorReactions) {

        auto slinker1 = r.productInfo.linkerBound1.substr(0, r.productInfo.linkerBound1.find(':'));
        auto slinker2 = r.productInfo.linkerBound2.substr(0, r.productInfo.linkerBound2.find(':'));
        auto sbound1 = r.reactantInfo.speciesBound1.substr(0, r.reactantInfo.speciesBound1.find(':'));
        auto sbound2 = r.reactantInfo.speciesBound2.substr(0, r.reactantInfo.speciesBound2.find(':'));

        protoCompartment.addSpeciesUnique(std::make_unique< Species >(
            SpeciesNamesDB::genBindingName(slinker1, sbound1, slinker2, sbound2),
            0, max_ulim,
            SpeciesType::pairBinding, RSpeciesType::REG
        ));
    }

}

void ChemManager::restartreleaseandremovaltime(floatingpoint _minimizationTime, ChemistryData& chemData){

    //look at copy number for each species
    for(auto &s : chemData.speciesDiffusing) {

        auto name = s.name;
        auto copyNumber = s.initialCopyNumber;
        auto releaseTime = s.releaseTime;
        auto removalTime = s.removalTime;
        if(tau()-releaseTime >= _minimizationTime) {
            //set zero copy number
            s.initialCopyNumber = 0;
        }
        if(tau()-removalTime >= _minimizationTime && !areEqual(removalTime,0.0) &&
           s.initialCopyNumber != -1) {
            ///set as removed by marking copy number to -1
            s.initialCopyNumber = -1;
        }
    }

}

void ChemManager::updateCopyNumbers(SubSystem& sys, ChemistryData& chemData, const medyan::SimulConfig& sc) {
    //Special protocol if move boundary protocol exists
    int tsaxis = sc.boundParams.transfershareaxis;
    floatingpoint cmpsize = 0.0;
    //X axis
    if(tsaxis == 0)
        cmpsize = sc.geoParams.compartmentSizeX;
        //Yaxis
    else if(tsaxis == 1)
        cmpsize = sc.geoParams.compartmentSizeY;
        //Z axis
    else if(tsaxis == 2)
        cmpsize = sc.geoParams.compartmentSizeZ;
    auto grid = sys.getCompartmentGrid();

    //look at copy number for each species
    for(auto &s : chemData.speciesDiffusing) {

        auto name = s.name;
        auto copyNumber = s.initialCopyNumber;
        auto releaseTime = s.releaseTime;
        auto removalTime = s.removalTime;
        auto cpynummanipulationType = s.copyNumberManipulationType;
        auto holdmolarity = s.holdMolarity;
        floatingpoint factor = sc.geoParams.compartmentSizeX * sc.geoParams
                        .compartmentSizeY * sc.geoParams.compartmentSizeZ * 6.023*1e-7;
        int updatedbasecopynumber = (int)(holdmolarity * factor);

        if(tau() >= releaseTime) {

            //add randomly in compartment
            while (copyNumber > 0) {

                //find a random compartment within the boundary
                Compartment* randomCompartment;
                while(true) {
                    randomCompartment = &GController::getRandomCompartment();
                    if(randomCompartment->isActivated()) break;
                }
                //find the species, increase copy number
                Species* species = randomCompartment->findSpeciesByName(name);

                //If user requests to use different copy numbers mentioned in chem file
                // during restart in place of those in restart file, do not do anything
                // in restart phase.
                if(sc.filamentSetup.USECHEMCOPYNUM && SysParams::RUNSTATE == false){
                    species->updateReactantPropensities();
                    copyNumber = 0;
                }
                else {
                    species->up();
                    species->updateReactantPropensities();
                    copyNumber--;

                    //set zero copy number
                    if (copyNumber == 0) s.initialCopyNumber = 0;
                }
            }

            //Change copy number if moveboundary is defined and if species is NOT past
            // removal.
            if(tsaxis >=0 && SysParams::RUNSTATE) {
                if (tsaxis < 3 && s.initialCopyNumber != -1 &&
                    cpynummanipulationType == "BASECONC") {
                    //set the coordinate that will help you find the necessary Base compartment
                    floatingpoint distancetocompare = 0.0;
                    if (sc.boundParams.planestomove == 2 &&
                        cpynummanipulationType != "NONE") {
                        cout << "Cannot set base concentration if both end planes are mobile as"
                                   " specified in BOUNDARYMOVE. Exiting." << endl;
                        exit(EXIT_FAILURE);
                    }
                        //if you are moving right, top or back boundaries, use left, bottom or
                        // front boundaries as the base.
                    else if (sc.boundParams.planestomove == 0)
                        distancetocompare = cmpsize / 2;
                        //if you are moving  left, bottom or front boundaries, use right, top or
                        // back boundaries as the base.
                    else if (sc.boundParams.planestomove == 1) {
                        floatingpoint systemspan = 0.0;
                        if (tsaxis == 0)
                            systemspan = sc.geoParams.NX * sc.geoParams
                                    .compartmentSizeX;
                        else if (tsaxis == 1)
                            systemspan = sc.geoParams.NY * sc.geoParams
                                    .compartmentSizeY;
                        else if (tsaxis == 2)
                            systemspan = sc.geoParams.NZ *
                                         sc.geoParams.compartmentSizeZ;
                        distancetocompare = -cmpsize / 2 + systemspan;
                    }
                    //Find the base compartment and set copy number.
                    for (auto& c:sys.getCompartmentGrid()->getCompartments()) {
                        if (c->coordinates()[tsaxis] == distancetocompare) {
                            //find the species, increase copy number
                            Species *species = c->findSpeciesByName(name);
                            int copynum = updatedbasecopynumber - species->getN();
                            while (copynum < 0) {
                                species->down();
                                species->updateReactantPropensities();
                                copynum++;
                            }
                            while (copynum > 0) {
                                species->up();
                                species->updateReactantPropensities();
                                copynum--;
                            }
                            std::cout << c->coordinates()[tsaxis] << " "
                                      << species->getFullName() << " " << species->getN()
                                      << endl;
                        }
                    }
                }
            }
        }

        if(SysParams::RUNSTATE && tau() >= removalTime && !areEqual(removalTime,0.0) &&
        s.initialCopyNumber != -1) {

            ///remove species from all compartments
            for(auto& C : grid->getCompartments()) {

                Species* species = C->findSpeciesByName(name);

                while(species->getN() > 0) {

                    species->down();
                    species->updateReactantPropensities();
                }
            }
            ///set as removed by marking copy number to -1
            s.initialCopyNumber = -1;
        }
    }

    for(auto &s : chemData.speciesBulk) {

        auto& name = s.name;
        auto& copyNumber = s.initialCopyNumber;
        auto& releaseTime = s.releaseTime;
        auto& removalTime = s.removalTime;
        auto& cpynummanipulationType = s.copyNumberManipulationType;
        auto& holdmolarity = s.holdMolarity;
        floatingpoint factor = Boundary::systemvolume * 6.023*1e-7;
        //If system is being restarted, do not update Copynumbers
        if(SysParams::RUNSTATE == false){
            //activate reactions
            //find the species, set copy number
            Species* species = grid->findSpeciesBulkByName(name);
            species->setN(0);
            species->activateReactantReactions();

        }
        if(SysParams::RUNSTATE && tau() >= releaseTime && copyNumber != 0) {

            //find the species, set copy number
            Species* species = grid->findSpeciesBulkByName(name);
            species->setN(copyNumber);

            //activate reactions
            species->activateReactantReactions();

            //set zero copy number
            s.initialCopyNumber = 0;
        }
        //if copy number changes with concentration
        if(SysParams::RUNSTATE && tau() >= releaseTime && cpynummanipulationType != "NONE"){
            //find the species, set copy number
            if(cpynummanipulationType == "BULKCONC") {
                int updatedcpynumber = (int)(holdmolarity * factor);
                Species *species = grid->findSpeciesBulkByName(name);
                species->setN(updatedcpynumber);

                //activate reactions
                species->activateReactantReactions();
            }

        }
        if(SysParams::RUNSTATE && tau() >= removalTime && !areEqual(removalTime,0.0) && s.initialCopyNumber != -1) {

            Species* species = grid->findSpeciesBulkByName(name);

            species->setN(0);

            //passivate reactions
            species->passivateReactantReactions();

            ///set as removed by marking copy number to -1
            s.initialCopyNumber = -1;
        }
    }
}

void ChemManager::genGeneralReactions(CompartmentGrid& grid, Compartment& protoCompartment, const ChemParams& chemParams, const ChemistryData& chemData) {

    //go through reactions, add each
    for(auto &r: chemData.genReactions) {

        vector<Species*> reactantSpecies;
        vector<Species*> productSpecies;

        vector<string> reactants = get<0>(r);
        vector<string> products = get<1>(r);

        int numDiffusingReactant = 0; // Used in determining volume dependence

        for(auto &reactant : reactants) {

            if(reactant.find("BULK") != string::npos) {

                //Look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find_if(chemData.speciesBulk.begin(), chemData.speciesBulk.end(),
                                  [&name](auto&& element) {
                                   return element.name == name; });

                if(it == chemData.speciesBulk.end()) {
                    cout <<
                    "A bulk species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                reactantSpecies.push_back(grid.findSpeciesBulkByName(name));
            }

            else if(reactant.find("DIFFUSING") != string::npos) {

                //Look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it =
                find_if(chemData.speciesDiffusing.begin(), chemData.speciesDiffusing.end(),
                        [&name](auto&& element) {
                            return element.name == name; });
                if(it == chemData.speciesDiffusing.end()) {
                    cout <<
                    "A diffusing species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                reactantSpecies.push_back(protoCompartment.findSpeciesByName(name));

                ++numDiffusingReactant;
            }
            else {
                cout <<
                "All reactants and products in a general reaction must be bulk or diffusing. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }

        for(auto &product : products) {

            if(product.find("BULK") != string::npos) {

                //Look up species, make sure in list
                string name = product.substr(0, product.find(":"));
                auto it = find_if(chemData.speciesBulk.begin(), chemData.speciesBulk.end(),
                                  [&name](auto&& element) {
                                      return element.name == name; });

                if(it == chemData.speciesBulk.end()) {
                    cout <<
                    "A bulk species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                productSpecies.push_back(grid.findSpeciesBulkByName(name));
            }

            else if(product.find("DIFFUSING") != string::npos) {

                //Look up species, make sure in list
                string name = product.substr(0, product.find(":"));
                auto it =
                find_if(chemData.speciesDiffusing.begin(), chemData.speciesDiffusing.end(),
                        [&name](auto&& element) {
                            return element.name == name; });
                if(it == chemData.speciesDiffusing.end()) {
                    cout <<
                    "A diffusing species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                productSpecies.push_back(protoCompartment.findSpeciesByName(name));
            }
            else {
                cout <<
                "All reactants and products in a general reaction must be bulk or diffusing. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        //add the reaction
        vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());

        ReactionBase* rxn;
        //create the reaction

        //<1,1>
        if(reactantSpecies.size() == 1 && productSpecies.size() == 1)
            rxn = new Reaction<1,1>(species, get<2>(r), true, 1.0, -numDiffusingReactant);
        //<2,1>
        else if(reactantSpecies.size() == 2 && productSpecies.size() == 1)
            rxn = new Reaction<2,1>(species, get<2>(r), true, 1.0, -numDiffusingReactant);
        //<1,2>
        else if(reactantSpecies.size() == 1 && productSpecies.size() == 2)
            rxn = new Reaction<1,2>(species, get<2>(r), true, 1.0, -numDiffusingReactant);
        //<2,0>
        else if(reactantSpecies.size() == 2 && productSpecies.size() == 0)
            rxn = new Reaction<2,0>(species, get<2>(r), true, 1.0, -numDiffusingReactant);
        //<2,2>
        else if(reactantSpecies.size() == 2 && productSpecies.size() == 2)
            rxn = new Reaction<2,2>(species, get<2>(r), true, 1.0, -numDiffusingReactant);
        //<1,3>
        else if(reactantSpecies.size() == 1 && productSpecies.size() == 3)
            rxn = new Reaction<1,3>(species, get<2>(r), true, 1.0, -numDiffusingReactant);
        //<2,2>
        else if(reactantSpecies.size() == 2 && productSpecies.size() == 3)
            rxn = new Reaction<2,3>(species, get<2>(r), true, 1.0, -numDiffusingReactant);
        //<3,2>
        else if(reactantSpecies.size() == 3 && productSpecies.size() == 2)
            rxn = new Reaction<3,2>(species, get<2>(r), true, 1.0, -numDiffusingReactant);
        else {
            cout <<
            "General reaction specified does not match any existing templates. Exiting."
            <<endl;
            exit(EXIT_FAILURE);
        }


        //add to compartment
        protoCompartment.addInternalReaction(rxn);
        rxn->setReactionType(ReactionType::REGULAR);

        // Dissipation
        if(chemParams.dissTracking){
            rxn->setGNumber(get<3>(r));
            rxn->setHRCDID(get<4>(r));
        }


    }
}

void ChemManager::genBulkReactions(CompartmentGrid& grid, const ChemistryData& chemData) {

    //go through reactions, add each
    for(auto &r: chemData.bulkReactions) {

        vector<Species*> reactantSpecies;
        vector<Species*> productSpecies;

        vector<string> reactants = get<0>(r);
        vector<string> products = get<1>(r);

        for(auto &reactant : reactants) {
            if(reactant.find("BULK") != string::npos) {

                //Look up species, make sure in list
                string name = reactant.substr(0, reactant.find(":"));
                auto it = find_if(chemData.speciesBulk.begin(), chemData.speciesBulk.end(),
                                  [&name](auto&& element) {
                                      return element.name == name; });

                if(it == chemData.speciesBulk.end()) {
                    cout <<
                    "A bulk species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                reactantSpecies.push_back(grid.findSpeciesBulkByName(name));
            }
            else {
                cout <<
                "All reactants and products in a bulk reaction must be bulk. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }

        for(auto &product : products) {
            if(product.find("BULK") != string::npos) {

                //Look up species, make sure in list
                string name = product.substr(0, product.find(":"));
                auto it = find_if(chemData.speciesBulk.begin(), chemData.speciesBulk.end(),
                                  [&name](auto&& element) {
                                      return element.name == name; });

                if(it == chemData.speciesBulk.end()) {
                    cout <<
                    "A bulk species that was included in a reaction was not initialized. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }
                productSpecies.push_back(grid.findSpeciesBulkByName(name));
            }
            else {
                cout <<
                "All reactants and products in a bulk reaction must be bulk. Exiting."
                << endl;
                exit(EXIT_FAILURE);
            }
        }
        //add the reaction
        vector<Species*> species = reactantSpecies;
        species.insert(species.end(), productSpecies.begin(), productSpecies.end());

        ReactionBase* rxn;
        //create the reaction

        //<1,1>
        if(reactantSpecies.size() == 1 && productSpecies.size() == 1)
            rxn = new Reaction<1,1>(species, get<2>(r));
        //<2,1>
        else if(reactantSpecies.size() == 2 && productSpecies.size() == 1)
            rxn = new Reaction<2,1>(species, get<2>(r));
        //<1,2>
        else if(reactantSpecies.size() == 1 && productSpecies.size() == 2)
            rxn = new Reaction<1,2>(species, get<2>(r));
        //<2,0>
        else if(reactantSpecies.size() == 2 && productSpecies.size() == 0)
            rxn = new Reaction<2,0>(species, get<2>(r));
        //<2,2>
        else if(reactantSpecies.size() == 2 && productSpecies.size() == 2)
            rxn = new Reaction<2,2>(species, get<2>(r));
        //<1,3>
        else if(reactantSpecies.size() == 1 && productSpecies.size() == 3)
            rxn = new Reaction<1,3>(species, get<2>(r));
        //<2,2>
        else if(reactantSpecies.size() == 2 && productSpecies.size() == 3)
            rxn = new Reaction<2,3>(species, get<2>(r));
        //<3,2>
        else if(reactantSpecies.size() == 3 && productSpecies.size() == 2)
            rxn = new Reaction<3,2>(species, get<2>(r));
        else {
            cout << "Bulk reaction specified does not match any existing templates. Exiting."
            <<endl;
            exit(EXIT_FAILURE);
        }

        //add to grid
        grid.addBulkReactionUnique(unique_ptr<ReactionBase>(rxn));
        rxn->setReactionType(ReactionType::REGULAR);
    }
}

void ChemManager::genNucleationReactions(SubSystem& sys, CompartmentGrid& grid, const medyan::SimulConfig& sc) {
    const auto& chemData = sc.chemistryData;

    for(int filType = 0; filType < sc.chemParams.numFilaments; filType++) {

        //loop through all compartments
        for(auto& C : grid.getCompartments()) {

            //go through reactions, add each
            for(auto &r: chemData.nucleationReactions[filType]) {
                
                //filament creation is not allowed in partially activated compartments
                //that volume fraction < threshold, be careful with the threshold
                if(C->getVolumeFrac() < 0.5) continue;

                vector<Species*> reactantSpecies;

                vector<string> reactants = get<0>(r);
                vector<string> products = get<1>(r);

                if(reactants.size() != NUCLEATIONREACTANTS ||
                   products.size() != NUCLEATIONPRODUCTS) {
                    cout << "Invalid nucleation reaction. Exiting." << endl;
                    exit(EXIT_FAILURE);
                }
                bool diffusing = false;

                int numDiffusingReactant = 0; // Used in determining volume dependence

                for(auto &reactant : reactants) {
                    if(reactant.find("BULK") != string::npos) {

                        //Look up species, make sure in list
                        string name = reactant.substr(0, reactant.find(":"));
                        auto it = find_if(chemData.speciesBulk.begin(), chemData.speciesBulk.end(),
                                          [&name](auto&& element) {
                                              return element.name == name; });

                        if(it == chemData.speciesBulk.end()) {
                            cout <<
                            "A bulk species that was included in a reaction was not initialized. Exiting."
                            << endl;
                            exit(EXIT_FAILURE);
                        }
                        reactantSpecies.push_back(grid.findSpeciesBulkByName(name));
                    }

                    else if(reactant.find("DIFFUSING") != string::npos) {

                        //Look up species, make sure in list
                        string name = reactant.substr(0, reactant.find(":"));
                        auto it =
                                find_if(chemData.speciesDiffusing.begin(), chemData.speciesDiffusing.end(),
                                        [&name](auto&& element) {
                                            return element.name == name; });
                        if(it == chemData.speciesDiffusing.end()) {
                            cout <<
                            "A diffusing species that was included in a reaction was not initialized. Exiting."
                            << endl;
                            exit(EXIT_FAILURE);
                        }
                        reactantSpecies.push_back(C->findSpeciesByName(name));
                        diffusing = true;

                        ++numDiffusingReactant;
                    }
                    else {
                        cout <<
                        "All reactants and products in a general reaction must be bulk or diffusing. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                }

                //add the reaction. The products will only be involved in creating the
                //callback needed to create a new filament
                ReactionBase* rxn = new Reaction<2,0>(reactantSpecies, get<2>(r), false, C->getVolumeFrac(), -numDiffusingReactant);
                rxn->setReactionType(ReactionType::FILAMENTCREATION);

                C->addInternalReaction(rxn);

                reactantSpecies.clear();

                //now, loop through products, add callback
                short plusEnd;
                short minusEnd;
                short filament;

                //FIRST SPECIES MUST BE PLUS END
                auto product = products[0];
                if(product.find("PLUSEND") != string::npos) {

                    //look up species, make sure in list
                    string name = product.substr(0, product.find(":"));
                    auto it = find(chemData.speciesPlusEnd[filType].begin(), chemData.speciesPlusEnd[filType].end(), name);

                    if(it != chemData.speciesPlusEnd[filType].end()) {
                        //get position of iterator
                        plusEnd = distance(chemData.speciesPlusEnd[filType].begin(), it);
                    }
                    else {
                        cout <<
                        "A plus end species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                }
                else {
                    cout <<
                    "First product species listed in a nucleation reaction must be plus end. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }

                //SECOND SPECIES MUST BE FILAMENT
                product = products[1];
                if(product.find("FILAMENT") != string::npos) {

                    //look up species, make sure in list
                    string name = product.substr(0, product.find(":"));
                    auto it = find(chemData.speciesFilament[filType].begin(), chemData.speciesFilament[filType].end(), name);

                    if(it != chemData.speciesFilament[filType].end()) {
                        //get position of iterator
                        filament = distance(chemData.speciesFilament[filType].begin(), it);
                    }
                    else {
                        cout <<
                        "A filament species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                }
                else {
                    cout <<
                    "Second product species listed in a nucleation reaction must be filament. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }

                //THIRD SPECIES MUST BE MINUSEND
                product = products[2];
                if(product.find("MINUSEND") != string::npos) {

                    //look up species, make sure in list
                    string name = product.substr(0, product.find(":"));
                    auto it = find(chemData.speciesMinusEnd[filType].begin(), chemData.speciesMinusEnd[filType].end(), name);

                    if(it != chemData.speciesMinusEnd[filType].end()) {
                        //get position of iterator
                        minusEnd = distance(chemData.speciesMinusEnd[filType].begin(), it);
                    }
                    else {
                        cout <<
                        "A minus end species that was included in a reaction was not initialized. Exiting."
                        << endl;
                        exit(EXIT_FAILURE);
                    }
                }
                else {
                    cout <<
                    "Third product species listed in a nucleation reaction must be minus end. Exiting."
                    << endl;
                    exit(EXIT_FAILURE);
                }

                //if the reaction had any diffusing species, create the filament
                //in a random position within that compartment
                Compartment* creationCompartment;
                if(diffusing) creationCompartment = C.get();
                else creationCompartment = &GController::getRandomCompartment();

                //now, add the callback
                FilamentCreationCallback
                        fcallback(plusEnd, minusEnd, filament, filType, &sys, creationCompartment);
                rxn->connect(fcallback);
            }
        }
    }
}

void ChemManager::initializeSystem(medyan::ChemSim* chemSim, SubSystem& sys, medyan::SimulConfig& sc) {

    auto& grid = *sys.getCompartmentGrid();
    auto& dt = *chemSim->getDT();

    setupBindingSites(sc.chemParams, sc);
    configCMonomer(sc);

    //Setup all species diffusing and bulk
    Compartment& cProto = grid.getProtoCompartment();

    genSpecies(grid, cProto, sc);

    //will print reactions as well
    genGeneralReactions(grid, cProto, sc.chemParams, sc.chemistryData);
    genBulkReactions(grid, sc.chemistryData);

    //initialize all compartments equivalent to cproto
    //will copy all general and bulk reactions
    for(auto& C : grid.getCompartments())
        *C = cProto;

    for(int cindex = 0; cindex < grid.getCompartments().size(); ++cindex) {
        medyan::generateAllDiffusionReactions(grid, cindex, true);
    }

    //try initial copy number setting
    updateCopyNumbers(sys, sc.chemistryData, sc);

    genNucleationReactions(sys, grid, sc);
    genFilBindingReactions(sys, sc, &dt);

    //add reactions in compartment grid to chemsim
    grid.addChemSimReactions(chemSim);

    genFilReactionTemplates(sc.chemParams, sc.chemistryData, &dt);

    // Backward compatibility section.
    //---------------------------------
    SysParams::CParams = sc.chemParams;
    chemDataBackup_ = sc.chemistryData;
}


void ChemManager::initializeCCylinder(
    CCylinder* cc,
    bool extensionFront,
    bool extensionBack,
    bool initialization,
    int nummonomers,
    int firstmonomer,
    int lastmonomer,
    bool minusendstatus,
    bool plusendstatus,
    short minusendtype,
    short plusendtype
) {

    mins = chrono::high_resolution_clock::now();
    //get some related objects
    Compartment* C = cc->getCompartment();
    Cylinder* c = cc->getCylinder();

    Filament* f = (Filament*)(c->getParent());
    short filType = f->getType();
    //add monomers to cylinder
    for(int i = 0; i < cc->getSize(); i++) {
        CMonomer* m = new CMonomer(filType);
        initCMonomer(m, filType, C, chemDataBackup_);
        cc->addCMonomer(m);

        if(find(SysParams::Chemistry().bindingSites[filType].begin(),
                SysParams::Chemistry().bindingSites[filType].end(), i)
           !=  SysParams::Chemistry().bindingSites[filType].end()) {

            //add callback to all binding sites
            UpdateBrancherBindingCallback bcallback(c, i);

            Species* bs = cc->getCMonomer(i)->speciesBound(
                    SysParams::Chemistry().brancherBoundIndex[filType]);
            bs->connect(bcallback);

            UpdateLinkerBindingCallback lcallback(c, i);

            Species* ls = cc->getCMonomer(i)->speciesBound(
                    SysParams::Chemistry().linkerBoundIndex[filType]);
            ls->connect(lcallback);

            UpdateMotorBindingCallback mcallback(c, i);

            Species* ms = cc->getCMonomer(i)->speciesBound(
                    SysParams::Chemistry().motorBoundIndex[filType]);
            ms->connect(mcallback);
        }
    }

    mine = chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_time1(mine - mins);
    tchemmanager1 += elapsed_time1.count();


    //get last ccylinder
    CCylinder* lastcc = nullptr;

    //extension of front
    if(extensionFront) {
        mins = chrono::high_resolution_clock::now();
        lastcc = f->getCylinderVector().back()->getCCylinder();
        for(auto &r : _filRxnTemplates[filType]) r->addReaction(lastcc, cc);
        mine = chrono::high_resolution_clock::now();
        chrono::duration<floatingpoint> elapsed_time2(mine - mins);
        tchemmanager2 += elapsed_time2.count();
    }
        //extension of back
    else if(extensionBack) {
        mins = chrono::high_resolution_clock::now();
        lastcc = f->getCylinderVector().front()->getCCylinder();
        for(auto &r : _filRxnTemplates[filType]) r->addReaction(cc, lastcc);
        mine = chrono::high_resolution_clock::now();
        chrono::duration<floatingpoint> elapsed_time2(mine - mins);
        tchemmanager2 += elapsed_time2.count();
    }

        //Base case, initialization
    else if (initialization) {
        mins = chrono::high_resolution_clock::now();
        //Check if this is the first cylinder
        if(!f->getCylinderVector().empty()) {

            //remove plus end from last, add to this.
            lastcc = f->getCylinderVector().back()->getCCylinder();

            if(SysParams::RUNSTATE){
                //Turn off PlusEnd Species in the previous cylinder  and set it as
                // SpeciesFilament
                CMonomer* m1 = lastcc->getCMonomer(lastcc->getSize() - 1);
                m1->speciesPlusEnd(0)->down();
                //fill last cylinder with default filament value
                m1->speciesFilament(0)->up();

                for(auto j : SysParams::Chemistry().bindingIndices[filType])
                    m1->speciesBound(j)->up();
                //Set the end of current cylinder to be Plus End
                CMonomer* m2 = cc->getCMonomer(cc->getSize() - 1);
                m2->speciesPlusEnd(0)->up();
                //fill new cylinder with default filament value
                for(int i = 0; i < cc->getSize() - 1; i++) {
                    cc->getCMonomer(i)->speciesFilament(0)->up();

                    for(auto j : SysParams::Chemistry().bindingIndices[filType])
                        cc->getCMonomer(i)->speciesBound(j)->up();
                }
            }
            else{/*RESTARTPHASE*/
                int start = firstmonomer;
                int end = lastmonomer;

                if(minusendstatus){
                    CMonomer *m2 = cc->getCMonomer(firstmonomer);
                    //minusendtype should have valid ( not -1 ) values if minusendstatus
                    // is non-negative.
                    m2->speciesMinusEnd(minusendtype)->up();
                    start = start+1;
                }
                if(plusendstatus) {
                    CMonomer *m2 = cc->getCMonomer(lastmonomer);
	                //plusendtype should have valid ( not -1 ) values if plusendstatus
	                // is non-negative.
                    m2->speciesPlusEnd(plusendtype)->up();
                    end = end -1;
                }

                //fill new cylinder with default filament value (start-end)
                for(int i = start; i <= end; i++) {
                    cc->getCMonomer(i)->speciesFilament(0)->up();

                    for(auto j : SysParams::Chemistry().bindingIndices[filType])
                        cc->getCMonomer(i)->speciesBound(j)->up();
                }
            }
            //Add cross cylinder reactions.
            for(auto &r : _filRxnTemplates[filType]) r->addReaction(lastcc, cc);
        }
            //this is first one
        else {
            if(SysParams::RUNSTATE){
                //set back and front
                CMonomer* m1 = cc->getCMonomer(cc->getSize() - 1);
                m1->speciesPlusEnd(0)->up();

                CMonomer* m2 = cc->getCMonomer(0);
                m2->speciesMinusEnd(0)->up();
                //fill with default filament value
                for(int i = 1; i < cc->getSize() - 1; i++) {
                    cc->getCMonomer(i)->speciesFilament(0)->up();

                    for(auto j : SysParams::Chemistry().bindingIndices[filType])
                        cc->getCMonomer(i)->speciesBound(j)->up();
                }
            }
            else {
	            int start = firstmonomer;
	            int end = lastmonomer;

                if(minusendstatus){
                    CMonomer *m2 = cc->getCMonomer(firstmonomer);
                    //minusendtype should have valid ( not -1 ) values if minusendstatus
                    // is non-negative.
                    m2->speciesMinusEnd(minusendtype)->up();
                    start = start+1;
                }
                if(plusendstatus) {
                    CMonomer *m2 = cc->getCMonomer(lastmonomer);
	                //plusendtype should have valid ( not -1 ) values if plusendstatus
	                // is non-negative.
	                m2->speciesPlusEnd(plusendtype)->up();
                    end = end -1;
                }

	            //fill new cylinder with default filament value (start-end)
	            for(int i = start; i <= end; i++) {
		            cc->getCMonomer(i)->speciesFilament(0)->up();

		            for(auto j : SysParams::Chemistry().bindingIndices[filType])
			            cc->getCMonomer(i)->speciesBound(j)->up();
	            }
            }
        }
        mine = chrono::high_resolution_clock::now();
        chrono::duration<floatingpoint> elapsed_time3(mine - mins);
        tchemmanager3 += elapsed_time3.count();

    }


    mins = chrono::high_resolution_clock::now();
    //Add all reaction templates to this cylinder
    for(auto &r : _filRxnTemplates[filType]) { r->addReaction(cc); }

    mine = chrono::high_resolution_clock::now();
    chrono::duration<floatingpoint> elapsed_time4(mine - mins);
    tchemmanager4 += elapsed_time4.count();
}

floatingpoint ChemManager::tchemmanager1 = 0.0;
floatingpoint ChemManager::tchemmanager2 = 0.0;
floatingpoint ChemManager::tchemmanager3 = 0.0;
floatingpoint ChemManager::tchemmanager4 = 0.0;

} // namespace medyan

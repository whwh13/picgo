#ifndef MEDYAN_Structure_Output_OChemistry_hpp
#define MEDYAN_Structure_Output_OChemistry_hpp

#include <type_traits>

#include "Structure/OutputStruct.hpp"
#include "Structure/SubSystem.h"
#include "Util/Io/H5.hpp"

namespace medyan {

inline auto getGlobalSpeciesNames(const SimulConfig& conf) {
    std::vector<std::string> res;

    // All diffusing and bulk species.
    for(auto& sd : conf.chemistryData.speciesDiffusing) {
        res.push_back(sd.name);
    }
    for(auto& sb : conf.chemistryData.speciesBulk) {
        res.push_back(sb.name);
    }

    // Membrane diffusing species.
    for(int si = 0; si < conf.chemistryData.speciesMembraneDiffusing.size(); ++si) {
        res.push_back(conf.chemistryData.speciesMembraneDiffusing[si].name);
    }

    // All filament related species.
    for(int filType = 0; filType < conf.chemParams.numFilaments; ++filType) {
        for(auto& s : conf.chemistryData.speciesFilament[filType]) {
            res.push_back(s);
        }
        for(auto& s : conf.chemistryData.speciesPlusEnd[filType]) {
            res.push_back(s);
        }
        for(auto& s : conf.chemistryData.speciesMinusEnd[filType]) {
            res.push_back(s);
        }
        for(auto& s : conf.chemistryData.speciesLinker[filType]) {
            res.push_back(s);
        }
        for(auto& s : conf.chemistryData.speciesMotor[filType]) {
            res.push_back(s);
        }
        for(auto& s : conf.chemistryData.speciesBrancher[filType]) {
            res.push_back(s);
        }
    }

    return res;
}
inline auto getGlobalSpeciesCopyNumbers(const SubSystem& sys, const SimulConfig& conf) {

    std::vector<Size> res;

    auto& grid = *sys.getCompartmentGrid();

    // All diffusing and bulk species.
    for(auto& sd : conf.chemistryData.speciesDiffusing) {
        res.push_back(grid.countDiffusingSpecies(sd.name));
    }
    for(auto& sb : conf.chemistryData.speciesBulk) {
        res.push_back(grid.countBulkSpecies(sb.name));
    }

    // Membrane diffusing species.
    for(int si = 0; si < conf.chemistryData.speciesMembraneDiffusing.size(); ++si) {
        Size count = 0;
        for(auto& v : sys.vertices) {
            count += v.cVertex.species.findSpeciesByIndex(si)->getRSpecies().getTrueN();
        }
        res.push_back(count);
    }

    // All filament related species.
    for(int filType = 0; filType < conf.chemParams.numFilaments; ++filType) {
        for(auto& s : conf.chemistryData.speciesFilament[filType]) {
            res.push_back(Filament::countSpecies(filType, s));
        }
        for(auto& s : conf.chemistryData.speciesPlusEnd[filType]) {
            res.push_back(Filament::countSpecies(filType, s));
        }
        for(auto& s : conf.chemistryData.speciesMinusEnd[filType]) {
            res.push_back(Filament::countSpecies(filType, s));
        }
        for(auto& s : conf.chemistryData.speciesLinker[filType]) {
            res.push_back(Linker::countSpecies(s));
        }
        for(auto& s : conf.chemistryData.speciesMotor[filType]) {
            res.push_back(MotorGhost::countSpecies(s));
        }
        for(auto& s : conf.chemistryData.speciesBrancher[filType]) {
            res.push_back(BranchingPoint::countSpecies(s));
        }
    }

    return res;
}

inline void extract(OutputStructChemistry& outChemistry, const SubSystem& sys, const SimulConfig& conf) {
    const auto& grid = *sys.getCompartmentGrid();
    const auto nc = grid.getCompartments().size();

    // Diffusing species.
    {
        const Size ns = conf.chemistryData.speciesDiffusing.size();
        outChemistry.diffusingSpeciesCopyNumbers.resize(std::array<Size, 4>{ns, grid.shape[0], grid.shape[1], grid.shape[2]});

        for(Index k = 0; k < grid.shape[2]; ++k) {
            for(Index j = 0; j < grid.shape[1]; ++j) {
                for(Index i = 0; i < grid.shape[0]; ++i) {
                    const auto& comp = grid.getCompartment(i, j, k);
                    for(Index s = 0; s < ns; ++s) {
                        const auto& species = *comp.findSpeciesByMolecule(SpeciesNamesDB::stringToInt(conf.chemistryData.speciesDiffusing[s].name));
                        outChemistry.diffusingSpeciesCopyNumbers(s, i, j, k) = species.getRSpecies().getTrueN();
                    }
                }
            }
        }
    }

    // Globally reported species.
    outChemistry.globalSpeciesCopyNumbers = getGlobalSpeciesCopyNumbers(sys, conf);
}

inline void write(h5::Group& grpSnapshot, const OutputStructChemistry& outChemistry) {
    // Diffusing species.
    h5::writeDataSet(grpSnapshot, "diffusingSpeciesCopyNumbers", outChemistry.diffusingSpeciesCopyNumbers);

    // Globally reported species.
    h5::writeDataSet(grpSnapshot, "globalSpeciesCopyNumbers", outChemistry.globalSpeciesCopyNumbers);
}

inline void read(OutputStructChemistry& outChemistry, const h5::Group& grpSnapshot) {
    // Diffusing species.
    h5::readDataSet(outChemistry.diffusingSpeciesCopyNumbers, grpSnapshot, "diffusingSpeciesCopyNumbers");

    // Globally reported species.
    h5::readDataSet(outChemistry.globalSpeciesCopyNumbers, grpSnapshot, "globalSpeciesCopyNumbers");
}

} // namespace medyan

#endif

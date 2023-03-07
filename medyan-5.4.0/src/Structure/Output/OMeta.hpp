#ifndef MEDYAN_Structure_Output_OMeta_hpp
#define MEDYAN_Structure_Output_OMeta_hpp

#include "Structure/Output/OBubble.hpp"
#include "Structure/Output/OChemistry.hpp"
#include "Structure/Output/OFilament.hpp"
#include "Structure/Output/OLinker.hpp"
#include "Structure/Output/OMembrane.hpp"

namespace medyan {

inline void extract(OutputStructMeta& outMeta, const SimulConfig& conf, const ForceFieldManager& ffm) {
    extract(outMeta.membraneMeta, conf.membraneMeshChemistryInfo);
    extract(outMeta.filamentMeta, conf.mechParams.globalFilamentModel);

    outMeta.simulConfig = conf;

    // Energies.
    outMeta.energyNames = { "total" };
    auto otherEnergyNames = ffm.getIndividualEnergyNames();
    outMeta.energyNames.insert(outMeta.energyNames.end(), otherEnergyNames.begin(), otherEnergyNames.end());

    // Chemistry.
    //----------------------------------
    {
        // Diffusing species names.
        outMeta.diffusingSpeciesNames.clear();
        for(auto& sd : conf.chemistryData.speciesDiffusing) {
            outMeta.diffusingSpeciesNames.push_back(sd.name);
        }
        // Global species names.
        outMeta.globalSpeciesNames = getGlobalSpeciesNames(conf);
    }
}

inline void write(h5::Group& grpHeader, const OutputStructMeta& outMeta) {
    auto grpMeta = grpHeader.createGroup("meta");
    write(grpMeta, outMeta.membraneMeta);
    write(grpMeta, outMeta.filamentMeta);

    // Write simulation config in strings.
    {
        const auto systemInputStr = [&] {
            std::ostringstream oss;
            outputTokenList(oss, buildTokens(outMeta.simulConfig.value(), systemParser()));
            return oss.str();
        }();
        const auto chemInputStr = [&] {
            std::ostringstream oss;
            outputTokenList(oss, buildTokens(outMeta.simulConfig.value(), chemDataParser()));
            return oss.str();
        }();

        h5::writeDataSet(grpMeta, "system", systemInputStr);
        h5::writeDataSet(grpMeta, "chem", chemInputStr);
    }

    h5::writeDataSet(grpMeta, "energyNames", outMeta.energyNames);

    // Chemistry.
    h5::writeDataSet(grpMeta, "diffusingSpeciesNames", outMeta.diffusingSpeciesNames);
    h5::writeDataSet(grpMeta, "globalSpeciesNames", outMeta.globalSpeciesNames);
}

inline void read(OutputStructMeta& outMeta, const h5::Group& grpHeader) {
    auto grpMeta = grpHeader.getGroup("meta");
    read(outMeta.membraneMeta, grpMeta);
    read(outMeta.filamentMeta, grpMeta);

    // Read simulation config from strings.
    {
        auto systemInputStr = h5::readDataSet<std::string>(grpMeta, "system");
        auto chemInputStr   = h5::readDataSet<std::string>(grpMeta, "chem");

        outMeta.simulConfig.emplace();
        parseKeyValueList(
            *outMeta.simulConfig,
            lexTokenList(sExprTokenize(systemInputStr)),
            systemParser()
        );
        parseKeyValueList(
            *outMeta.simulConfig,
            lexTokenList(sExprTokenize(chemInputStr)),
            chemDataParser()
        );
    }

    h5::readDataSet(outMeta.energyNames, grpMeta, "energyNames");

    // Chemistry.
    h5::readDataSet(outMeta.diffusingSpeciesNames, grpMeta, "diffusingSpeciesNames");
    h5::readDataSet(outMeta.globalSpeciesNames, grpMeta, "globalSpeciesNames");
}

} // namespace medyan

#endif

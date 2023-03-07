#ifndef MEDYAN_Structure_Output_OFilament_hpp
#define MEDYAN_Structure_Output_OFilament_hpp

#include "Structure/OutputStruct.hpp"
#include "Structure/SubSystem.h"
#include "Util/Io/H5.hpp"

namespace medyan {

// Filament meta.
inline void extract(OutputStructFilamentMeta& outFilMeta, FilamentModel globalFilamentModel) {
    outFilMeta.globalFilamentModel = globalFilamentModel;
}
inline void write(h5::Group& grpFilMeta, const OutputStructFilamentMeta& outFilMeta) {
    const auto globalFilamentModelStr = toString(outFilMeta.globalFilamentModel);
    h5::writeDataSet(grpFilMeta, "globalFilamentModel", globalFilamentModelStr);
}
inline void read(OutputStructFilamentMeta& outFilMeta, const h5::Group& grpFilMeta) {
    auto globalFilamentModelStr = h5::readDataSet<std::string>(grpFilMeta, "globalFilamentModel");
    parse(outFilMeta.globalFilamentModel, globalFilamentModelStr);
}


inline void extract(OutputStructFilament& outFil, FilamentModel globalFilamentModel, const Filament& fil) {
    auto& cylvec = fil.getCylinderVector();

    outFil.id = fil.getId();
    outFil.type = fil.getType();
    outFil.numBeads = cylvec.size() + 1;
    outFil.deltaMinusEnd = fil.getDeltaMinusEnd();
    outFil.deltaPlusEnd = fil.getDeltaPlusEnd();

    if(globalFilamentModel == FilamentModel::beadCylinder) {
        outFil.rawCoords.resize(3, outFil.numBeads);
        const auto addBeadCoord = [&](const Bead& b, int bi) {
            for(int dim = 0; dim < 3; ++dim) {
                outFil.rawCoords(dim, bi) = b.coordinate()[dim];
            }
        };
        for(int ci = 0; ci < cylvec.size(); ++ci) {
            addBeadCoord(*cylvec[ci]->getFirstBead(), ci);
        }
        addBeadCoord(*cylvec.back()->getSecondBead(), cylvec.size());
    }
}
inline void extract(std::vector<OutputStructFilament>& outFils, FilamentModel globalFilamentModel, const SubSystem& sys) {
    outFils.clear();
    outFils.reserve(Filament::getFilaments().size());
    for(auto pf : Filament::getFilaments()) {
        OutputStructFilament outFil;
        extract(outFil, globalFilamentModel, *pf);
        outFils.push_back(std::move(outFil));
    }
}

inline void write(h5::Group& grpFils, const OutputStructFilament& outFil, int filIndex) {
    // Write a filament.
    //----------------------------------
    h5::writeDataSet(grpFils, to_string(filIndex), outFil.rawCoords);
}
inline void write(h5::Group& grpSnapshot, const std::vector<OutputStructFilament>& outFils) {
    using namespace h5;

    Group grpFils = grpSnapshot.createGroup("filaments");

    // Store number of filaments.
    h5::writeDataSet<std::int64_t>(grpFils, "count", outFils.size());

    // Store filament data.
    for(int i = 0; i < outFils.size(); ++i) {
        write(grpFils, outFils[i], i);
    }
}

inline void read(OutputStructFilament& outFil, const h5::Group& grpFils, int filIndex) {
    // Read a filament.
    //----------------------------------
    h5::readDataSet(outFil.rawCoords, grpFils, to_string(filIndex));
}
inline void read(std::vector<OutputStructFilament>& outFils, const h5::Group& grpSnapshot) {
    using namespace h5;

    Group grpFils = grpSnapshot.getGroup("filaments");

    // Read number of filaments.
    auto count = h5::readDataSet<std::int64_t>(grpFils, "count");

    // Read filament data.
    outFils.resize(count);
    for(int i = 0; i < count; ++i) {
        read(outFils[i], grpFils, i);
    }
}

} // namespace medyan

#endif

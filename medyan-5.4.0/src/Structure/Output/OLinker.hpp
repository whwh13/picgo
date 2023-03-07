#ifndef MEDYAN_Structure_Output_OLinker_hpp
#define MEDYAN_Structure_Output_OLinker_hpp

#include "Structure/OutputStruct.hpp"
#include "Structure/SubSystem.h"
#include "Util/Io/H5.hpp"

namespace medyan {

inline void extract(OutputStructLinker& outLinker, const Linker& linker) {
    outLinker.id = linker.getId();
    outLinker.type = "linker";
    outLinker.subtype = linker.getType();

    // Store coordinates
    const auto coord1 = midPointCoordinate(
        linker.getFirstCylinder()->getFirstBead()->vcoordinate(),
        linker.getFirstCylinder()->getSecondBead()->vcoordinate(),
        linker.getFirstCylinder()->adjustedrelativeposition(linker.getFirstPosition())
    );
    const auto coord2 = midPointCoordinate(
        linker.getSecondCylinder()->getFirstBead()->vcoordinate(),
        linker.getSecondCylinder()->getSecondBead()->vcoordinate(),
        linker.getSecondCylinder()->adjustedrelativeposition(linker.getSecondPosition())
    );
    outLinker.rawCoords(0, 0) = coord1[0];
    outLinker.rawCoords(1, 0) = coord1[1];
    outLinker.rawCoords(2, 0) = coord1[2];
    outLinker.rawCoords(0, 1) = coord2[0];
    outLinker.rawCoords(1, 1) = coord2[1];
    outLinker.rawCoords(2, 1) = coord2[2];
}
inline void extract(OutputStructLinker& outLinker, const MotorGhost& linker) {
    outLinker.id = linker.getId();
    outLinker.type = "motor";
    outLinker.subtype = linker.getType();

    // Store coordinates
    const auto coord1 = midPointCoordinate(
        linker.getFirstCylinder()->getFirstBead()->vcoordinate(),
        linker.getFirstCylinder()->getSecondBead()->vcoordinate(),
        linker.getFirstCylinder()->adjustedrelativeposition(linker.getFirstPosition())
    );
    const auto coord2 = midPointCoordinate(
        linker.getSecondCylinder()->getFirstBead()->vcoordinate(),
        linker.getSecondCylinder()->getSecondBead()->vcoordinate(),
        linker.getSecondCylinder()->adjustedrelativeposition(linker.getSecondPosition())
    );
    outLinker.rawCoords(0, 0) = coord1[0];
    outLinker.rawCoords(1, 0) = coord1[1];
    outLinker.rawCoords(2, 0) = coord1[2];
    outLinker.rawCoords(0, 1) = coord2[0];
    outLinker.rawCoords(1, 1) = coord2[1];
    outLinker.rawCoords(2, 1) = coord2[2];
}
inline void extract(OutputStructLinker& outLinker, const BranchingPoint& linker) {
    outLinker.id = linker.getId();
    outLinker.type = "brancher";
    outLinker.subtype = linker.getType();

    // Store coordinates
    const auto coord1 = linker.coordinate;
    const auto coord2 = linker.getSecondCylinder()->getFirstBead()->coordinate();
    outLinker.rawCoords(0, 0) = coord1[0];
    outLinker.rawCoords(1, 0) = coord1[1];
    outLinker.rawCoords(2, 0) = coord1[2];
    outLinker.rawCoords(0, 1) = coord2[0];
    outLinker.rawCoords(1, 1) = coord2[1];
    outLinker.rawCoords(2, 1) = coord2[2];
}
inline void extract(std::vector<OutputStructLinker>& outLinkers, const SubSystem& sys) {
    outLinkers.clear();
    for(auto pl : Linker::getLinkers()) {
        OutputStructLinker outLinker;
        extract(outLinker, *pl);
        outLinkers.push_back(std::move(outLinker));
    }
    for(auto pl : MotorGhost::getMotorGhosts()) {
        OutputStructLinker outLinker;
        extract(outLinker, *pl);
        outLinkers.push_back(std::move(outLinker));
    }
    for(auto pl : BranchingPoint::getBranchingPoints()) {
        OutputStructLinker outLinker;
        extract(outLinker, *pl);
        outLinkers.push_back(std::move(outLinker));
    }
}

inline void write(h5::Group& grpLinkers, const OutputStructLinker& outLinker, int linkerIndex) {
    auto grpLinker = grpLinkers.createGroup(toString(linkerIndex));
    h5::writeDataSet<std::int64_t>(grpLinker, "id",      outLinker.id);
    h5::writeDataSet<std::string> (grpLinker, "type",    outLinker.type);
    h5::writeDataSet<std::int64_t>(grpLinker, "subtype", outLinker.subtype);

    h5::writeDataSet              (grpLinker, "coords", outLinker.rawCoords);
}
inline void write(h5::Group& grpSnapshot, const std::vector<OutputStructLinker>& outLinkers) {
    using namespace h5;

    Group grpLinkers = grpSnapshot.createGroup("linkers");

    // Write number of linkers.
    h5::writeDataSet<std::int64_t>(grpLinkers, "count", outLinkers.size());

    // Write linker data.
    for(int i = 0; i < outLinkers.size(); ++i) {
        write(grpLinkers, outLinkers[i], i);
    }
}

inline void read(OutputStructLinker& outLinker, const h5::Group& grpLinkers, int linkerIndex) {
    auto grpLinker = grpLinkers.getGroup(toString(linkerIndex));
    outLinker.id      = h5::readDataSet<std::int64_t>(grpLinker, "id");
    outLinker.type    = h5::readDataSet<std::string> (grpLinker, "type");
    outLinker.subtype = h5::readDataSet<std::int64_t>(grpLinker, "subtype");

    h5::readDataSet(outLinker.rawCoords, grpLinker, "coords");
}
inline void read(std::vector<OutputStructLinker>& outLinkers, const h5::Group& grpSnapshot) {
    using namespace h5;

    Group grpLinkers = grpSnapshot.getGroup("linkers");

    // Read number of linkers.
    auto count = h5::readDataSet<std::int64_t>(grpLinkers, "count");

    // Read linker data.
    outLinkers.resize(count);
    for(int i = 0; i < count; ++i) {
        read(outLinkers[i], grpLinkers, i);
    }
}

} // namespace medyan

#endif

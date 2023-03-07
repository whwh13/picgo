#ifndef MEDYAN_Structure_Output_OSnapshot_hpp
#define MEDYAN_Structure_Output_OSnapshot_hpp

#include "Structure/Output/OBubble.hpp"
#include "Structure/Output/OChemistry.hpp"
#include "Structure/Output/OFilament.hpp"
#include "Structure/Output/OLinker.hpp"
#include "Structure/Output/OMembrane.hpp"

namespace medyan {

inline void extract(OutputStructSnapshot& outSnapshot, const SubSystem& sys, const SimulConfig& conf, int frame) {
    outSnapshot.snapshot = frame;
    outSnapshot.simulationTime = sys.tau();
    extract(outSnapshot.filamentStruct, conf.mechParams.globalFilamentModel, sys);
    extract(outSnapshot.linkerStruct, sys);
    extract(outSnapshot.bubbleStruct, sys);
    extract(outSnapshot.membraneStruct, conf.membraneMeshChemistryInfo, sys);

    // Chemistry.
    extract(outSnapshot.chemistry, sys, conf);

    // Energies.
    {
        auto& ereport = sys.prevMinResult.energiesAfter;
        outSnapshot.energies.resize(ereport.individual.size() + 1);
        outSnapshot.energies[0] = ereport.total;
        for(int i = 0; i < ereport.individual.size(); ++i) {
            outSnapshot.energies[i + 1] = ereport.individual[i].energy;
        }
    }
}

// The function that appends to the output function.
inline void write(h5::Group& grpSnapshots, const OutputStructSnapshot& outSnapshot) {
    using namespace std;
    using namespace h5;

    Group grpSnapshot = grpSnapshots.createGroup(to_string(outSnapshot.snapshot));

    h5::writeDataSet<double>(grpSnapshot, "time", outSnapshot.simulationTime);

    // Write to file.
    write(grpSnapshot, outSnapshot.filamentStruct);
    write(grpSnapshot, outSnapshot.linkerStruct);
    write(grpSnapshot, outSnapshot.bubbleStruct);
    write(grpSnapshot, outSnapshot.membraneStruct);

    // Chemistry.
    write(grpSnapshot, outSnapshot.chemistry);

    // Energies.
    h5::writeDataSet(grpSnapshot, "energies", outSnapshot.energies);
}

inline void read(OutputStructSnapshot& outSnapshot, const h5::Group& grpSnapshots, Index snapshot) {
    using namespace h5;
    Group grpSnapshot = grpSnapshots.getGroup(toString(snapshot));
    outSnapshot.snapshot = snapshot;

    outSnapshot.simulationTime = h5::readDataSet<double>(grpSnapshot, "time");

    read(outSnapshot.filamentStruct, grpSnapshot);
    read(outSnapshot.linkerStruct, grpSnapshot);
    read(outSnapshot.bubbleStruct, grpSnapshot);
    read(outSnapshot.membraneStruct, grpSnapshot);

    // Chemistry.
    read(outSnapshot.chemistry, grpSnapshot);

    // Energies.
    h5::readDataSet(outSnapshot.energies, grpSnapshot, "energies");
}

} // namespace medyan

#endif

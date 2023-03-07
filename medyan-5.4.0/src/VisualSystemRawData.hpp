#ifndef MEDYAN_VisualSystemRawData_hpp
#define MEDYAN_VisualSystemRawData_hpp

#include <array>
#include <cstdint>
#include <mutex>
#include <optional>
#include <vector>

#include "Mechanics/ForceField/ForceFieldManager.h"
#include "Structure/OutputStruct.hpp"
#include "Structure/SubSystem.h"
#include "Util/Math/Vec.hpp"
#include "Visual/FrameData.hpp"

namespace medyan::visual {

namespace raw_data_cat {

    using Type = std::uint_fast8_t;

    constexpr Type none            = 0;
    constexpr Type beadPosition    = 1 << 0;
    constexpr Type beadConnection  = 1 << 1;
    constexpr Type compartment     = 1 << 2;
    constexpr Type concentration   = 1 << 3;

} // namespace raw_data_cat

// SystemRawData is a stateful data container, used for storing system data
// that is useful for visualization.
//
// The raw data is updated in the simulation thread with system data.
// The data can be consumed so that the update state is reset.
struct SystemRawData {

    // Synchronization and states
    std::mutex me;

    raw_data_cat::Type updated = raw_data_cat::none; // Should be reset by the consumer

    // Raw meta data.
    std::optional<OutputStructMeta> outMeta;

    // Processed meta data.
    DisplayTypeMap             displayTypeMap;

    // Data
    DisplayFrame               frameData;

};


// The global variable for shared data
inline SystemRawData sdfv;


// This function should be called after system initialization, but before copySystemData is ever called.
void copySystemMetaData(
    SystemRawData& data,
    const SubSystem& sys,
    const ForceFieldManager& ffm,
    const SimulConfig& conf,
    const CommandLineConfig& cmdConfig
);
// Copy system meta data to the global shared data.
inline void copySystemMetaData(
    const SubSystem& sys,
    const ForceFieldManager& ffm,
    const SimulConfig& conf,
    const CommandLineConfig& cmdConfig
) {
    return copySystemMetaData(sdfv, sys, ffm, conf, cmdConfig);
}


// Function to copy the system data to raw data
//
// Note:
//   - This function should only be called on the simulation thread.
//   - This function will not perform any action is GUI is not turned on.
//
// Input:
//   - data: the raw data object
//   - updated: the part of the system that's needs to be updated
//   - ignoreDataInUse: if set to true, copy would be skipped if the data is
//     being used by other threads
//
// Returns whether the copy is actually successfully made.
bool copySystemData(
    SystemRawData&     data,
    const SubSystem&   sys,
    const SimulConfig& conf,
    const CommandLineConfig& cmdConfig,
    raw_data_cat::Type updated,
    bool               ignoreDataInUse = true
);

// Copy system data to the global shared data.
inline bool copySystemData(
    const SubSystem&   sys,
    const SimulConfig& conf,
    const CommandLineConfig& cmdConfig,
    raw_data_cat::Type updated = 0,
    bool               ignoreDataInUse = true
) {
    return copySystemData(sdfv, sys, conf, cmdConfig, updated, ignoreDataInUse);
}



} // namespace medyan::visual

#endif

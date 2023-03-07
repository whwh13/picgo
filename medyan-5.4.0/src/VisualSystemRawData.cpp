#include "VisualSystemRawData.hpp"

#include "Output.h"
#include "SysParams.h"

namespace medyan::visual {

void copySystemMetaData(
    SystemRawData& data,
    const SubSystem& sys,
    const ForceFieldManager& ffm,
    const SimulConfig& conf,
    const CommandLineConfig& cmdConfig
) {
    // Do not copy data if GUI is not built or turned on.
    if(!(builtGui && cmdConfig.guiEnabled)) {
        return;
    }

    // Acquire the lock.
    std::scoped_lock lk(data.me);

    // Copy the meta data.
    data.outMeta.emplace();
    extract(*data.outMeta, conf, ffm);
}

bool copySystemData(
    SystemRawData&     data,
    const SubSystem&   sys,
    const SimulConfig& conf,
    const CommandLineConfig& cmdConfig,
    raw_data_cat::Type updated,
    bool               ignoreDataInUse
) {
    // Do not copy data if GUI is not built or turned on.
    if(!(builtGui && cmdConfig.guiEnabled)) {
        return false;
    }

    std::unique_lock lk(data.me, std::try_to_lock);
    if(!lk.owns_lock()) {
        if(ignoreDataInUse) {
            // The data is in use, skip copying
            return false;
        } else {
            // Acquire the lock
            lk.lock();
        }
    }

    // Do not copy data if the meta data is not set.
    if(!data.outMeta.has_value()) {
        log::warn("copySystemData: meta data has not been set.");
        return false;
    }

    // Extract system data.
    // Maybe future: use "updated" to filter out the data that's not needed.
    OutputStructSnapshot outSnapshot;
    extract(outSnapshot, sys, conf, 0); // Frame not used, so set to 0.

    // Set frame data.
    data.frameData = readOneFrameDataFromOutput(*data.outMeta, outSnapshot, data.displayTypeMap);

    // Save updated.
    data.updated |= updated;

    return true;
}

} // namespace medyan::visual

#ifndef MEDYAN_Structure_Output_OBubble_hpp
#define MEDYAN_Structure_Output_OBubble_hpp

#include "Structure/OutputStruct.hpp"
#include "Structure/SubSystem.h"
#include "Util/Io/H5.hpp"

namespace medyan {

inline void extract(OutputStructBubble& outBubble, const Bubble& bubble) {
    outBubble.id = bubble.getId();
    outBubble.type = bubble.getType();
    outBubble.coords[0] = bubble.coord[0];
    outBubble.coords[1] = bubble.coord[1];
    outBubble.coords[2] = bubble.coord[2];
    outBubble.radius = bubble.getRadius();
}
inline void extract(std::vector<OutputStructBubble>& outBubbles, const SubSystem& sys) {
    outBubbles.clear();
    for(auto& b : sys.bubbles) {
        auto& outBubble = outBubbles.emplace_back();
        extract(outBubble, b);
    }
}

inline void write(h5::Group& grpBubbles, const OutputStructBubble& outBubble, int bubbleIndex) {
    auto grpBubble = grpBubbles.createGroup(toString(bubbleIndex));
    h5::writeDataSet<std::int64_t>(grpBubble, "id",      outBubble.id);
    h5::writeDataSet<std::int64_t>(grpBubble, "type",    outBubble.type);
    h5::writeDataSet              (grpBubble, "coords",  outBubble.coords);
    h5::writeDataSet              (grpBubble, "radius",  outBubble.radius);
}
inline void write(h5::Group& grpSnapshot, const std::vector<OutputStructBubble>& outBubbles) {
    auto grpBubbles = grpSnapshot.createGroup("bubbles");

    // Number of bubbles.
    h5::writeDataSet<std::int64_t>(grpBubbles, "count", outBubbles.size());

    // Bubble data.
    for(int i = 0; i < outBubbles.size(); ++i) {
        write(grpBubbles, outBubbles[i], i);
    }
}

inline void read(OutputStructBubble& outBubble, const h5::Group& grpBubbles, int bubbleIndex) {
    auto grpBubble = grpBubbles.getGroup(toString(bubbleIndex));
    outBubble.id   = h5::readDataSet<std::int64_t>(grpBubble, "id");
    outBubble.type = h5::readDataSet<std::int64_t>(grpBubble, "type");
    h5::readDataSet(outBubble.coords, grpBubble, "coords");
    h5::readDataSet(outBubble.radius, grpBubble, "radius");
}
inline void read(std::vector<OutputStructBubble>& outBubbles, const h5::Group& grpSnapshot) {
    auto grpBubbles = grpSnapshot.getGroup("bubbles");

    // Number of bubbles.
    auto count = h5::readDataSet<std::int64_t>(grpBubbles, "count");

    // Bubble data.
    outBubbles.resize(count);
    for(int i = 0; i < count; ++i) {
        read(outBubbles[i], grpBubbles, i);
    }
}

} // namespace medyan

#endif

#ifndef MEDYAN_Visualization_hpp
#define MEDYAN_Visualization_hpp

// This file should be included when using visualization in main codes.

#include "Util/Io/Log.hpp"

#ifdef NO_GUI
// GUI functions will be excluded, allowing for fewer dependencies in compilation.
namespace medyan {
    // Print warning messages on GUI running attempts
    void guiRunDisabledWarning() {
        log::warn("GUI is not built into this version. Possibly because NO_GUI was specified in compilation.");
    }

    void guiRunRealtime()   { guiRunDisabledWarning(); }
    void guiRunTrajectory() { guiRunDisabledWarning(); }

} // namespace medyan

#else // #ifdef NO_GUI
// NO_GUI is not defined. GUI functions are enabled.
#include "Visual/Window.hpp"

namespace medyan {
    void guiRunInitialMode(medyan::visual::DisplayMode displayMode) {
        medyan::visual::VisualDisplay visualDisplay(displayMode);
        visualDisplay.run();
    }

    void guiRunRealtime()   { guiRunInitialMode(medyan::visual::DisplayMode::realtime); }
    void guiRunTrajectory() { guiRunInitialMode(medyan::visual::DisplayMode::trajectory); }

} // namespace medyan

#endif // #ifdef NO_GUI


#endif

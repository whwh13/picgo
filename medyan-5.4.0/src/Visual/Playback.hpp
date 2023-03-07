#ifndef MEDYAN_Visual_Playback_hpp
#define MEDYAN_Visual_Playback_hpp

#include <algorithm>

#include "Visual/DisplaySettings.hpp"
#include "Visual/DisplayStates.hpp"

namespace medyan::visual {

// Update max number of frames of the playback
inline void playbackUpdateMaxFrames(
    DisplayStates& states
) {
    int maxFrames = 0;
    for(const auto& traj : states.trajectoryDataStates.trajectories) {
        maxFrames = std::max<int>(maxFrames, traj.data.frames.size());
    }

    states.playback.maxFrame = max(0, maxFrames - 1);
}


// check whether a new frame should be issued
template< typename NewFrameFunc >
inline void playbackNewFrame(
    const TrajectoryPlaybackSettings& settings,
    TrajectoryPlaybackStates&         states,
    float                             glfwTime,
    NewFrameFunc&&                    func
) {
    if(states.offscreenRender.has_value()) {
        if(states.currentFrame < states.offscreenRender->frameRangeHi) {
            ++states.currentFrame;
        }
        else {
            states.offscreenRender.reset();
            log::info("Offscreen render complete.");
        }
    }
    else if(states.isPlaying) {
        if(glfwTime >= states.lastGlfwTime + 1.0f / settings.fps) {
            if(states.currentFrame < states.maxFrame) ++states.currentFrame;
        }
    }

    if(states.currentFrame != states.previousFrame) {

        // Execute functions for the new frame
        func();

        states.lastGlfwTime = glfwTime;
        states.previousFrame = states.currentFrame;
    }
}

// Check whether a new frame is issued.
// If yes, it will mark all mesh data to be updated.
inline void playbackCheckTrajectory(
    const DisplaySettings& settings,
    DisplayStates&         states,
    float                  glfwTime
) {
    playbackNewFrame(
        settings.playback,
        states.playback,
        glfwTime,
        [&] {
            for(auto& traj : states.trajectoryDataStates.trajectories) {
                for(auto& profileData : traj.profileData) {
                    profileData.shouldUpdateMeshData = true;
                }
            }
        }
    );
}

} // namespace medyan::visual

#endif

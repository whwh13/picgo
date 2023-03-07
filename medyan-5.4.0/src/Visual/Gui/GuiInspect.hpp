#ifndef MEDYAN_Visual_Gui_GuiInspect_hpp
#define MEDYAN_Visual_Gui_GuiInspect_hpp

#include <numeric>

#include <implot.h>

#include "Visual/DisplaySettings.hpp"
#include "Visual/DisplayStates.hpp"
#include "Visual/Gui/GuiAux.hpp"

namespace medyan::visual {

inline void guiInspectWindow(
    DisplaySettings& displaySettings,
    DisplayStates  & displayStates
) {
    if(!displaySettings.gui.inspectWindow) return;

    auto& trajs = displayStates.trajectoryDataStates.trajectories;
    const auto curFrame = displayStates.playback.currentFrame;
    const auto maxFrame = displayStates.playback.maxFrame;

    ImGui::SetNextWindowSize(ImVec2(700, 500), ImGuiCond_FirstUseEver);

    if (ImGui::Begin("inspect", &displaySettings.gui.inspectWindow)) {
        if(displaySettings.displayMode == DisplayMode::realtime) {
            ImGui::Text("Realtime energy inspection is not yet supported.");
        }
        else {
            // Plot energies as a function of frame.
            if(ImPlot::BeginPlot("energies", "frame", "energy (zJ)")) {
                // Prepare x-axis.
                static std::vector<double> xs;
                if(xs.size() < maxFrame + 1) {
                    xs.resize(maxFrame + 1);
                    for(int i = 0; i <= maxFrame; ++i) {
                        xs[i] = i;
                    }
                }

                // Prepare data.
                for(int i = 0; i < trajs.size(); ++i) {
                    auto& traj = trajs[i];
                    if(traj.displayMasterSwitch) {
                        const auto numEnergies = traj.data.energyNames.size();
                        traj.energyPlotSwitch.resize(numEnergies);
                        for(int ei = 0; ei < numEnergies; ++ei) {
                            if(traj.energyPlotSwitch[ei]) {
                                char label[32];
                                std::snprintf(label, sizeof(label), "%d-%s", i, traj.data.energyNames[ei].c_str());
                                ImPlot::PlotLine(label, xs.data(), traj.data.energyValues.col(ei).data(), traj.data.energyValues.rows());
                            }
                        }
                    }
                }

                // Plot current frame.
                ImPlot::PlotVLines("current frame", &curFrame, 1);

                ImPlot::EndPlot();
            }

            // Select energies for inspection.
            ImGui::Text("Select energies to inspect:");
            for(int i = 0; i < trajs.size(); ++i) {
                auto& traj = trajs[i];
                if(traj.displayMasterSwitch) {
                    char trajName[32];
                    std::snprintf(trajName, sizeof(trajName), "traj %d", i);
                    if(ImGui::TreeNode(trajName)) {
                        // Print energy names and values as selectables.
                        const auto numEnergies = traj.data.energyNames.size();
                        traj.energyPlotSwitch.resize(numEnergies);
                        for(int ei = 0; ei < numEnergies; ++ei) {
                            ImGui::Selectable(traj.data.energyNames[ei].c_str(), reinterpret_cast<bool*>(&traj.energyPlotSwitch[ei]));
                            ImGui::SameLine(300);
                            const auto nrows = traj.data.energyValues.rows();
                            if(curFrame < nrows) {
                                ImGui::Text("%f", traj.data.energyValues(curFrame, ei));
                            }
                            else {
                                ImGui::Text("-");
                            }
                        }

                        ImGui::TreePop();
                    }
                }
            }
        }
    }

    // End of window.
    ImGui::End();
}

} // namespace medyan::visual

#endif

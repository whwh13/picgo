#ifndef MEDYAN_Visual_Gui_Gui_hpp
#define MEDYAN_Visual_Gui_Gui_hpp

#include <cstdio>
#include <type_traits>

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <imgui_internal.h>
#include <imgui_stdlib.h>
#include <implot.h>
#include <nfd.h>

#include "Visual/DisplaySettings.hpp"
#include "Visual/DisplayStates.hpp"
#include "Visual/FrameData.hpp"
#include "Visual/Gui/GuiAux.hpp"
#include "Visual/Gui/GuiInspect.hpp"
#include "Visual/Gui/GuiKeyMapping.hpp"
#include "Visual/Gui/GuiProfile.hpp"

namespace medyan::visual {


// Note:
//   - This must be used while OpenGL and GLFW environments are live
class ImguiGuard {

private:
    ImGuiContext*  context_ = nullptr;
    ImPlotContext* plotContext_ = nullptr;

public:
    ImguiGuard(GLFWwindow* window) {
        // Setup Dear ImGui context
        IMGUI_CHECKVERSION();
        context_     = ImGui::CreateContext();
        plotContext_ = ImPlot::CreateContext();
        ImGuiIO& io = ImGui::GetIO(); (void)io;
        //io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
        //io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls

        // Setup Dear ImGui style
        ImGui::StyleColorsDark();
        //ImGui::StyleColorsClassic();

        // Setup Platform/Renderer bindings
        ImGui_ImplGlfw_InitForOpenGL(window, true);
        ImGui_ImplOpenGL3_Init("#version 330 core");

        // Load Fonts
        // - If no fonts are loaded, dear imgui will use the default font. You can also load multiple fonts and use ImGui::PushFont()/PopFont() to select them.
        // - AddFontFromFileTTF() will return the ImFont* so you can store it if you need to select the font among multiple.
        // - If the file cannot be loaded, the function will return NULL. Please handle those errors in your application (e.g. use an assertion, or display an error and quit).
        // - The fonts will be rasterized at a given size (w/ oversampling) and stored into a texture when calling ImFontAtlas::Build()/GetTexDataAsXXXX(), which ImGui_ImplXXXX_NewFrame below will call.
        // - Read 'docs/FONTS.md' for more instructions and details.
        // - Remember that in C/C++ if you want to include a backslash \ in a string literal you need to write a double backslash \\ !
        //io.Fonts->AddFontDefault();
        //io.Fonts->AddFontFromFileTTF("../../misc/fonts/Roboto-Medium.ttf", 16.0f);
        //io.Fonts->AddFontFromFileTTF("../../misc/fonts/Cousine-Regular.ttf", 15.0f);
        //io.Fonts->AddFontFromFileTTF("../../misc/fonts/DroidSans.ttf", 16.0f);
        //io.Fonts->AddFontFromFileTTF("../../misc/fonts/ProggyTiny.ttf", 10.0f);
        //ImFont* font = io.Fonts->AddFontFromFileTTF("c:\\Windows\\Fonts\\ArialUni.ttf", 18.0f, NULL, io.Fonts->GetGlyphRangesJapanese());
        //IM_ASSERT(font != NULL);
    }

    ~ImguiGuard() {
        // Cleanup
        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImGui::DestroyContext(context_);
        ImPlot::DestroyContext(plotContext_);
    }
};


class ImguiDisableHelperGuard {
    bool yesPleaseReallyDisableItForReal_;

public:
    ImguiDisableHelperGuard(bool forReal) : yesPleaseReallyDisableItForReal_(forReal) {
        if(yesPleaseReallyDisableItForReal_) {
            ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
            ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
        }
    }

    ~ImguiDisableHelperGuard() {
        if(yesPleaseReallyDisableItForReal_) {
            ImGui::PopItemFlag();
            ImGui::PopStyleVar();
        }
    }
};



//-----------------------------------------------------------------------------
// Main GUI functions
//-----------------------------------------------------------------------------

inline void guiHelpWindow(DisplaySettings& displaySettings) {
    if(!displaySettings.gui.helpWindow) return;

    ImGui::SetNextWindowSize(ImVec2(360, 400), ImGuiCond_FirstUseEver);

    if (ImGui::Begin("help", &displaySettings.gui.helpWindow)) {

        if(ImGui::CollapsingHeader("controls", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::Text("sliders: ctrl + click to manually input values.");
            ImGui::Text("key controls:");
            ImGui::BulletText("%s: take snapshot of current screen", display(displaySettings.mainView.control.keyMapping.keyBindingList[underlying(KeyCallbackFunction::takeSnapshot)]).c_str());
            ImGui::BulletText("%s: toggle gui on/off", display(displaySettings.mainView.control.keyMapping.keyBindingList[underlying(KeyCallbackFunction::toggleGui)]).c_str());
            ImGui::BulletText("%s: toggle play/pause", display(displaySettings.mainView.control.keyMapping.keyBindingList[underlying(KeyCallbackFunction::togglePlayPause)]).c_str());
            ImGui::BulletText("w/a/s/d: move camera horizontally");
        }
    }


    // End of help window, no matter Begin returns true or false.
    ImGui::End();
}



inline void guiProfileWindow(
    DisplaySettings& displaySettings,
    DisplayStates&   displayStates
) {
    if(!displaySettings.gui.profileWindow) return;

    ImGui::SetNextWindowSize(ImVec2(750, 300), ImGuiCond_FirstUseEver);

    if (ImGui::Begin("profile", &displaySettings.gui.profileWindow)) {

        static int selectedTrajectoryIndex = -1; // used only in trajectory mode. -1 means not selected
        static int selectedProfileIndex = -1; // used only when profiles are listed. -1 means not selected

        // left pane, for trajectory list
        {
            ImGui::BeginChild("trajectory pane", ImVec2(200, 0), true);

            if(displaySettings.displayMode == DisplayMode::trajectory) {
                if(ImGui::CollapsingHeader("new file", ImGuiTreeNodeFlags_Framed | ImGuiTreeNodeFlags_DefaultOpen)) {

                    static DisplayTrajectoryFileSettings newFile;

                    const auto fileSelection = [&](std::filesystem::path& file, std::filesystem::path& lastPath, const char* buttonName) {
                        if(ImGui::Button(buttonName)) {
                            guiAuxSelectFile(file, lastPath);
                            lastPath = file.parent_path();
                        }
                        if(!file.empty()) {
                            ImGui::SameLine();
                            char clearButtonLabelId[128];
                            std::snprintf(clearButtonLabelId, 128, "x##clear %s", buttonName);
                            if(ImGui::Button(clearButtonLabelId)) {
                                file.clear();
                            }

                            ImGui::TextWrapped(file.string().c_str());
                        }
                    };

                    fileSelection(newFile.trajSnapshot, newFile.lastTrajPath, "snapshot.traj");

                    // The "add" button
                    {
                        ImguiDisableHelperGuard guard(displayStates.sync.busy.load());

                        if(ImGui::Button("add")) {
                            pushAnAsyncTask(
                                displayStates.sync,
                                // file information is copied but not referenced
                                [&, file = newFile] {
                                    backgroundTaskReadTrajectory(displayStates.sync, file);
                                }
                            );
                        }
                    }
                }
            }

            if(displaySettings.displayMode == DisplayMode::realtime) {
                // always select realtime
                ImGui::Selectable("realtime", true);
            }
            else {
                const int numTrajectories = displayStates.trajectoryDataStates.trajectories.size();
                for (int i = 0; i < numTrajectories; ++i)
                {
                    auto& thisTraj = displayStates.trajectoryDataStates.trajectories[i];

                    char label[64];
                    std::snprintf(label, 64, "%d %s", i, thisTraj.displayMasterSwitch ? "" : "disabled");
                    if (ImGui::Selectable(label, selectedTrajectoryIndex == i)) {
                        selectedTrajectoryIndex = i;
                    }

                    if(ImGui::IsItemHovered() && ImGui::IsMouseClicked(ImGuiMouseButton_Right)) {
                        thisTraj.displayMasterSwitch ^= true;
                        if(thisTraj.displayMasterSwitch) {
                            // mark mesh data to be updated
                            for(auto& profileData : thisTraj.profileData) {
                                profileData.shouldUpdateMeshData = true;
                            }
                        }
                    }
                }
            }

            ImGui::EndChild();
        }
        ImGui::SameLine();

        // middle pane, for trajectory configuration and profiles
        {
            ImGui::BeginChild("config pane", ImVec2(250, 0), true);

            // The function to display all profiles
            const auto displayProfiles = [&](std::vector< ProfileWithMeshData >& vecProfileData) {

                // add profiles
                if(ImGui::Button("add profile")) {
                    vecProfileData.emplace_back();
                }

                const int numProfiles = vecProfileData.size();

                for(int i = 0; i < numProfiles; ++i) {
                    char label[128];
                    std::snprintf(label, 128, "%d %s", i, elementProfileDisplayName(vecProfileData[i].profile.index()));
                    if(ImGui::Selectable(label, selectedProfileIndex == i)) {
                        selectedProfileIndex = i;
                    }
                }
            };

            if(displaySettings.displayMode == DisplayMode::trajectory) {
                auto& trajectories = displayStates.trajectoryDataStates.trajectories;
                selectedTrajectoryIndex = std::min(selectedTrajectoryIndex, (int)trajectories.size() - 1);

                if(selectedTrajectoryIndex >= 0) {
                    // files
                    if(ImGui::CollapsingHeader("file")) {
                        ImGui::TextWrapped("snapshot file: %s", trajectories[selectedTrajectoryIndex].inputs.trajSnapshot.string().c_str());
                        if(ImGui::Button("close")) {
                            // close this file
                            trajectories.erase(trajectories.begin() + selectedTrajectoryIndex);
                        }
                    }

                    // profiles
                    if(ImGui::CollapsingHeader("profiles", ImGuiTreeNodeFlags_DefaultOpen)) {
                        if(selectedTrajectoryIndex < displayStates.trajectoryDataStates.trajectories.size()) {
                            displayProfiles(displayStates.trajectoryDataStates.trajectories[selectedTrajectoryIndex].profileData);
                        }
                    }
                }
            }
            else {
                if(ImGui::CollapsingHeader("profiles", ImGuiTreeNodeFlags_DefaultOpen)) {
                    displayProfiles(displayStates.realtimeDataStates.profileData);
                }
            }

            ImGui::EndChild();
        }
        ImGui::SameLine();

        // right pane, for profile configuration
        {
            ImGui::BeginChild("profile pane", ImVec2(0, 0), true);

            if(displaySettings.displayMode == DisplayMode::trajectory) {
                const bool valid =
                    selectedTrajectoryIndex >= 0 && selectedTrajectoryIndex < displayStates.trajectoryDataStates.trajectories.size() &&
                    selectedProfileIndex    >= 0 && selectedProfileIndex    < displayStates.trajectoryDataStates.trajectories[selectedTrajectoryIndex].profileData.size();

                if(valid) {
                    auto& traj = displayStates.trajectoryDataStates.trajectories[selectedTrajectoryIndex];
                    if(ImGui::Button("remove profile")) {
                        traj.profileData.erase(traj.profileData.begin() + selectedProfileIndex);
                    }
                    else {
                        auto& profileData = traj.profileData[selectedProfileIndex];
                        const auto frame = profileData.frameInfo.has_value() ? displayStates.playback.currentFrame : -1;
                        if(guiActiveProfileConfig(profileData.profile, traj.data.displayTypeMap, frame)) {
                            profileData.shouldUpdateMeshData = true;
                        }
                    }
                }
            }
            else {
                const bool valid =
                    selectedProfileIndex >= 0 && selectedProfileIndex < displayStates.realtimeDataStates.profileData.size();

                if(valid) {
                    if(ImGui::Button("remove profile")) {
                        displayStates.realtimeDataStates.profileData.erase(displayStates.realtimeDataStates.profileData.begin() + selectedProfileIndex);
                    }
                    else {
                        auto& profileData = displayStates.realtimeDataStates.profileData[selectedProfileIndex];
                        {
                            std::scoped_lock lock(sdfv.me);
                            if(guiActiveProfileConfig(profileData.profile, sdfv.displayTypeMap)) {
                                profileData.shouldUpdateMeshData = true;
                            }
                        }
                    }
                }
            }

            ImGui::EndChild();
        }
    }

    // End of help window, no matter Begin returns true or false.
    ImGui::End();
}

inline void guiTrajectoryControlPanel(
    DisplaySettings& displaySettings,
    DisplayStates&   displayStates
) {
    // Playback
    ImGui::Checkbox("play", &displayStates.playback.isPlaying);
    ImGui::SliderFloat("speed", &displaySettings.playback.fps, 1, 24, "%.1f");
    ImGui::SliderInt("playback", &displayStates.playback.currentFrame, 0, displayStates.playback.maxFrame);
}

inline void guiViewSettings(ObjectViewSettings& viewSettings, ObjectViewStates& viewStates, const DisplaySettings& settings, const DisplayStates& states) {
    guiAuxColorPicker4Popup("background", viewSettings.canvas.bgColor.data());
    // Depth cueing
    ImGui::Checkbox("depth cueing", &viewSettings.canvas.depthCueing);
    if(viewSettings.canvas.depthCueing) {
        const float minDiffDepth = 0.1;
        const bool dMinChanged = ImGui::SliderFloat("depth cueing min", &viewSettings.canvas.depthCueingDepthMin, 0.0, 0.8, "%.2f");
        if(dMinChanged) {
            viewSettings.canvas.depthCueingDepthMax = std::max(viewSettings.canvas.depthCueingDepthMax, viewSettings.canvas.depthCueingDepthMin + minDiffDepth);
        }
        const bool dMaxChanged = ImGui::SliderFloat("depth cueing max", &viewSettings.canvas.depthCueingDepthMax, 0.2, 1.0, "%.2f");
        if(dMaxChanged) {
            viewSettings.canvas.depthCueingDepthMin = std::min(viewSettings.canvas.depthCueingDepthMin, viewSettings.canvas.depthCueingDepthMax - minDiffDepth);
        }
    }

    if(ImGui::TreeNode("transformation")) {
        auto& view = viewSettings.camera.view;
        auto& proj = viewSettings.camera.projection;

        if(ImGui::Button("auto adjust camera")) {
            autoAdjustCamera(viewSettings.camera, settings, states);
        }

        // View settings.
        const bool distanceChanged = ImGui::SliderFloat("target distance", &view.distance, 100.0, 20000.0, "%.1f");

        auto cameraPos = view.position();
        ImGui::Text(
            "camera position: (%.1f, %.1f, %.1f)",
            cameraPos[0],
            cameraPos[1],
            cameraPos[2]
        );
        ImGui::Text(
            "camera target: (%.1f, %.1f, %.1f)",
            view.target[0],
            view.target[1],
            view.target[2]
        );

        // Projection settings.
        guiAuxEnumComboBox("projection", proj.type);

        const float minDiffZNearFar = 10;
        const bool zNearChanged = ImGui::SliderFloat("z near", &proj.zNear, 0.1, 2000.0, "%.1f");
        if(zNearChanged) {
            proj.zFar = std::max(proj.zFar, proj.zNear + minDiffZNearFar);
        }
        const bool zFarChanged = ImGui::SliderFloat("z far", &proj.zFar, 100.0, 25000.0, "%.1f");
        if(zFarChanged) {
            proj.zNear = std::min(proj.zNear, proj.zFar - minDiffZNearFar);
        }

        if(proj.type == ObjectViewSettings::Camera::Projection::Type::perspective) {
            ImGui::Text("fov: %.2f", proj.fov);
        } else {
            ImGui::Text("scale: %.1f", proj.scale);
        }

        ImGui::SliderFloat("camera key speed", &viewSettings.control.cameraKeyPositionPerFrame, 50.0, 2000.0, "%.1f");

        guiAuxEnumComboBox("mouse mode", viewSettings.control.cameraMouseMode);

        ImGui::TreePop();
    }

    if(ImGui::TreeNode("snapshot")) {
        guiAuxEnumComboBox("snapshot resolution", viewSettings.control.snapshotResolution);
        if(viewSettings.control.snapshotResolution == ObjectViewSettings::Control::SnapshotResolution::scaleWithScreen) {
            ImGui::SliderFloat("snapshot scale", &viewSettings.control.snapshotScale, 0.1f, 10.0f, "%.1f");
            ImGui::Checkbox("undo scaling orthographic", &viewSettings.control.snapshotUndoScaleOrtho);
        }
        else {
            ImGui::SliderInt("snapshot width", &viewSettings.control.snapshotWidth, 200, 1920);
            ImGui::SliderInt("snapshot height", &viewSettings.control.snapshotHeight, 150, 1080);
        }

        if(ImGui::Button("snapshot directory")) {
            guiAuxSelectFile(viewSettings.control.snapshotSavePath, viewSettings.control.snapshotSavePath.parent_path(), true);
        }
        ImGui::TextWrapped("%s", viewSettings.control.snapshotSavePath.string().c_str());

        ImGui::TreePop();
    }

    if(ImGui::TreeNode("keyboard")) {
        guiComponentKeyMapping(viewSettings.control.keyMapping, viewStates.control);

        ImGui::TreePop();
    }

}

inline void guiMainWindow(
    DisplaySettings& displaySettings,
    DisplayStates  & displayStates
) {
    const bool busy = displayStates.sync.busy.load();

    // Exceptionally add an extra assert here for people confused about initial Dear ImGui setup
    // Most ImGui functions would normally just crash if the context is missing.
    IM_ASSERT(ImGui::GetCurrentContext() != NULL && "Missing dear imgui context.");

    ImGuiWindowFlags windowFlags =
        ImGuiWindowFlags_MenuBar;
        // | ImGuiWindowFlags_NoCollapse;

    // We specify a default position/size in case there's no data in the .ini file.
    // We only do it to make the demo applications a little more welcoming, but typically this isn't required.
    ImGui::SetNextWindowPos(ImVec2(0, 0), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(400, 600), ImGuiCond_FirstUseEver);

    // Main body of the Demo window starts here.
    if (!ImGui::Begin("medyan control", nullptr, windowFlags))
    {
        // Early out if the window is collapsed, as an optimization.
        ImGui::End();
        return;
    }

    // Most "big" widgets share a common width settings by default. See 'Demo->Layout->Widgets Width' for details.

    // e.g. Use 2/3 of the space for widgets and 1/3 for labels (default)
    // ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.65f);

    // e.g. Leave a fixed amount of width for labels (by passing a negative value), the rest goes to widgets.
    // ImGui::PushItemWidth(ImGui::GetFontSize() * -12);

    // Menu Bar
    if (ImGui::BeginMenuBar())
    {
        if (ImGui::BeginMenu("file")) {
            ImGui::MenuItem("(empty menu)", NULL, false, false);
            ImGui::EndMenu();
        }

        if(ImGui::BeginMenu("window")) {
            ImGui::MenuItem("profile", nullptr, &displaySettings.gui.profileWindow, true);
            ImGui::MenuItem("inspect", nullptr, &displaySettings.gui.inspectWindow, true);
            ImGui::EndMenu();
        }

        if(ImGui::BeginMenu("help"))
        {
            ImGui::MenuItem("help window", nullptr, &displaySettings.gui.helpWindow, true);
            ImGui::EndMenu();
        }
        ImGui::EndMenuBar();
    }

    // Main display settings
    guiAuxEnumComboBox(
        "display mode",
        displaySettings.displayMode,
        [&](DisplayMode modeOld, DisplayMode modeNew) {
            switch(modeNew) {
                case DisplayMode::trajectory:
                    if(modeOld != modeNew && displayStates.trajectoryDataStates.trajectories.empty()) {
                        pushAnAsyncTask(
                            displayStates.sync,
                            [&] {
                                backgroundTaskReadTrajectory(displayStates.sync, DisplayTrajectoryFileSettings {});
                            }
                        );
                    }
                    break;
            }
        },
        busy ? ImGuiSelectableFlags_Disabled : 0
    );

    if(displaySettings.displayMode == DisplayMode::trajectory) {
        guiTrajectoryControlPanel(displaySettings, displayStates);
    }

    ImGui::Spacing();


    if (ImGui::CollapsingHeader("info", ImGuiTreeNodeFlags_DefaultOpen)) {

        // Function to print view object information
        const auto printView = [](const ObjectViewSettings& viewSettings) {
            ImGui::Text("view");
            ImGui::BulletText(
                "window size: (%d, %d)",
                viewSettings.canvas.width,
                viewSettings.canvas.height
            );
        };

        // busy information
        if(busy) {
            ImGui::Text("busy...");
            ImGui::Separator();
        }
        // fps information
        ImGui::Text("fps: %.1f", displayStates.timing.fps);
        ImGui::Separator();
        // main view information
        printView(displaySettings.mainView);

    }

    if(ImGui::CollapsingHeader("settings", ImGuiTreeNodeFlags_DefaultOpen)) {

        // main view
        guiViewSettings(displaySettings.mainView, displayStates.mainView, displaySettings, displayStates);
    }


    // End of main window
    ImGui::End();
}

// Note:
//   - This should only be used in GLFW main loop
inline void imguiLoopRender(
    DisplaySettings& displaySettings,
    DisplayStates  & displayStates
) {
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    if(displaySettings.gui.enabled) {
        guiMainWindow(displaySettings, displayStates);
        guiProfileWindow(displaySettings, displayStates);
        guiInspectWindow(displaySettings, displayStates);
        guiHelpWindow(displaySettings);
    }

    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

} // namespace medyan::visual

#endif

#ifndef MEDYAN_Visual_Window_Hpp
#define MEDYAN_Visual_Window_Hpp

#include <algorithm>
#include <array>
#include <functional>
#include <iostream> // cout, endl
#include <stdexcept> // runtime_error
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_WRITE_STATIC
#include "stb_image_write.h"

#include "Util/Environment.hpp"
#include "Util/Io/Log.hpp"
#include "Visual/Common.hpp"
#include "Visual/Control.hpp"
#include "Visual/Gui/Gui.hpp"
#include "Visual/Playback.hpp"
#include "Visual/Shader.hpp"
#include "Visual/ShaderSrc.hpp"
#include "Visual/VisualElement.hpp"
#include "VisualSystemRawData.hpp"


// For best portability, the window signal handling could only be done from the
// main thread (due to MacOS Cocoa framework).

namespace medyan::visual {


inline void glfwError(int id, const char* description) {
    log::error(description);
    throw std::runtime_error("Error in GLFW environment");
}



// The RAII object for managing visualization context and window
//
// Note:
//   - Only one object of this type can be created at a time.
//   - The window signal handling should only be done from the main thread (due
//     to MacOS Cocoa framework).
class VisualContext {
public:
    struct OffscreenBufferGuard {
        bool bufferBuilt = false;
        unsigned offscreenFbo = 0;
        unsigned offscreenRbo = 0;
        unsigned offscreenColorRbo = 0;

        // destructor will deallocate the frame buffers
        ~OffscreenBufferGuard() {
            unbuild();
        }

        // build a new buffer
        void build(int width, int height) {
            if(bufferBuilt) unbuild();

            bufferBuilt = true;

            glGenFramebuffers(1, &offscreenFbo);
            glBindFramebuffer(GL_FRAMEBUFFER, offscreenFbo);

            glGenRenderbuffers(1, &offscreenRbo);
            glBindRenderbuffer(GL_RENDERBUFFER, offscreenRbo);
            glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, width, height);

            glGenRenderbuffers(1, &offscreenColorRbo);
            glBindRenderbuffer(GL_RENDERBUFFER, offscreenColorRbo);
            glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA8, width, height);

            glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, offscreenRbo);
            glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, offscreenColorRbo);

            if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
                log::error("Framebuffer is not complete.");
            }
        }

        void unbuild() {
            if(bufferBuilt) {
                glDeleteFramebuffers(1, &offscreenFbo);
                glDeleteRenderbuffers(1, &offscreenRbo);
                glDeleteRenderbuffers(1, &offscreenColorRbo);
                bufferBuilt = false;
            }
        }
    };



    // Display settings and states
    DisplaySettings displaySettings;
    DisplayStates   displayStates;


    VisualContext(
        GLFWframebuffersizefun framebufferSizeCallback,
        GLFWcursorposfun       cursorPosCallback,
        GLFWmousebuttonfun     mouseButtonCallback,
        GLFWscrollfun          scrollCallback,
        GLFWkeyfun             keyCallback
    ) {
        // GLFW initializing
        log::debug("Initializing GLFW");
        glfwSetErrorCallback(&glfwError);
        glfwInit();
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
        // The following line is generally useless, but it is fun to see a transparent window.
        // glfwWindowHint(GLFW_TRANSPARENT_FRAMEBUFFER, GL_TRUE);

        // Window initializing
        window_ = glfwCreateWindow(displaySettings.mainView.canvas.width, displaySettings.mainView.canvas.height, "MEDYAN", NULL, NULL);
        if(window_ == NULL) {
            log::error("Failed to create GLFW window");
            glfwTerminate();
            return;
        }
        glfwMakeContextCurrent(window_);

        // GLAD initializing
        log::debug("initializing GLAD");
        if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
            log::error("Failed to initialize GLAD");
            return;
        }

        // Obtain the actual framebuffer size.
        glfwGetFramebufferSize(window_, &displaySettings.mainView.canvas.width, &displaySettings.mainView.canvas.height);
        glViewport(0, 0, displaySettings.mainView.canvas.width, displaySettings.mainView.canvas.height);

        // Initialize key mapping
        //-----------------------------
        displaySettings.mainView.control.keyMapping = KeyMapping::createDefault();

        // Set window callbacks
        //-----------------------------
        glfwSetFramebufferSizeCallback(window_, framebufferSizeCallback);
        glfwSetCursorPosCallback(window_, cursorPosCallback);
        glfwSetMouseButtonCallback(window_, mouseButtonCallback);
        glfwSetScrollCallback(window_, scrollCallback);
        glfwSetKeyCallback(window_, keyCallback);

    } // VisualContext()

    ~VisualContext() {
        glfwTerminate();
    }

    auto window() const { return window_; }

    // Helper function to process window inputs
    void processInput() {

        auto& camera = displaySettings.mainView.camera;

        const float currentTime = glfwGetTime();
        const float deltaTime = currentTime - displayStates.timing.glfwTimeLastFrame;
        displayStates.timing.update(currentTime);
        const float cameraMove = displaySettings.mainView.control.cameraKeyPositionPerFrame * deltaTime;

        if (glfwGetKey(window_, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
            // Do nothing
        }

        if(glfwGetKey(window_, GLFW_KEY_W) == GLFW_PRESS) {
            const auto change = -cameraMove * camera.view.backward();
            camera.view.target += change;
        }
        if(glfwGetKey(window_, GLFW_KEY_S) == GLFW_PRESS) {
            const auto change = cameraMove * camera.view.backward();
            camera.view.target += change;
        }
        if(glfwGetKey(window_, GLFW_KEY_A) == GLFW_PRESS) {
            const auto change = -cameraMove * camera.view.right();
            camera.view.target += change;
        }
        if(glfwGetKey(window_, GLFW_KEY_D) == GLFW_PRESS) {
            const auto change = cameraMove * camera.view.right();
            camera.view.target += change;
        }
    }

private:

    // Member variables
    GLFWwindow* window_;

}; // VisualContext


// The RAII object for all the rendering process
struct VisualDisplay {
    // Callbacks bound to key press.
    //---------------------------------
    static void keyCallbackTakeSnapshot(DisplaySettings*, DisplayStates* pstates) {
        pstates->mainView.control.snapshotRenderingNextFrame = true;
    }
    static void keyCallbackTakeSnapshotsAll(DisplaySettings* psettings, DisplayStates* pstates) {
        if(
            auto& offscreenRender = pstates->playback.offscreenRender;
            psettings->displayMode == DisplayMode::trajectory &&
            !offscreenRender.has_value()
        ) {
            offscreenRender.emplace(
                TrajectoryPlaybackStates::OffscreenRender {
                    0,
                    pstates->playback.maxFrame
                }
            );
            // Reset the playhead to be the previous frame of the minimum of the range of frames,
            // even if the frame does not exist.
            pstates->playback.currentFrame = offscreenRender->frameRangeLo - 1;
        }
    }
    static void keyCallbackToggleGui(DisplaySettings* psettings, DisplayStates*) {
        psettings->gui.enabled = !psettings->gui.enabled;
    }
    static void keyCallbackTogglePlayPause(DisplaySettings*, DisplayStates* pstates) {
        pstates->playback.isPlaying = !pstates->playback.isPlaying;
    }
    static void keyCallbackControlMouseSetRotate(DisplaySettings* psettings, DisplayStates*) {
        psettings->mainView.control.cameraMouseMode = ObjectViewSettings::Control::CameraMouseMode::rotate;
    }
    static void keyCallbackControlMouseSetPan(DisplaySettings* psettings, DisplayStates*) {
        psettings->mainView.control.cameraMouseMode = ObjectViewSettings::Control::CameraMouseMode::pan;
    }

    // Static key callback function list
    inline static const std::array< std::function<void(DisplaySettings*, DisplayStates*)>, underlying(KeyCallbackFunction::last_) > keyCallbackList {
        keyCallbackTakeSnapshot,
        keyCallbackTakeSnapshotsAll,
        keyCallbackToggleGui,
        keyCallbackTogglePlayPause,
        keyCallbackControlMouseSetRotate,
        keyCallbackControlMouseSetPan,
    };


    // GLFW window callbacks.
    //---------------------------------
    static void framebufferSizeCallback(GLFWwindow* window, int width, int height) {
        auto& vd = *static_cast< VisualDisplay* >(glfwGetWindowUserPointer(window));
        auto& canvas = vd.vc.displaySettings.mainView.canvas;
        canvas.width = width;
        canvas.height = height;
        glViewport(0, 0, width, height);
    }
    static void cursorPosCallback(GLFWwindow* window, double xpos, double ypos) {
        auto& vd = *static_cast< VisualDisplay* >(glfwGetWindowUserPointer(window));
        auto& controlStates = vd.vc.displayStates.mainView.control;
        auto& camera = vd.vc.displaySettings.mainView.camera;
        const auto& control = vd.vc.displaySettings.mainView.control;
        auto& io = ImGui::GetIO();

        const int mouseState = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
        if(mouseState == GLFW_PRESS) {
            if(controlStates.mouseLeftAlreadyPressed) {
                // Already in the dragging process.
                if(controlStates.mouseLeftDragStartInScene) {
                    mouseDragCameraView(
                        vd.vc.displaySettings.mainView.camera.view,
                        vd.vc.displayStates.mainView.control.mouseLastX,
                        vd.vc.displayStates.mainView.control.mouseLastY,
                        xpos,
                        ypos,
                        vd.vc.displaySettings.mainView.control,
                        vd.vc.displaySettings.mainView.camera.projection,
                        vd.vc.displaySettings.mainView.canvas
                    );
                }

            } else {
                // Starts mouse dragging.
                controlStates.mouseLeftAlreadyPressed = true;
                controlStates.mouseLeftDragStartInScene = !io.WantCaptureMouse; // True if mouse is not in any GUI window.
            }
            controlStates.mouseLastX = xpos;
            controlStates.mouseLastY = ypos;
        } else {
            controlStates.mouseLeftAlreadyPressed = false;
        }
    }
    static void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
        auto& vd = *static_cast< VisualDisplay* >(glfwGetWindowUserPointer(window));

        if (action == GLFW_PRESS) {
            // Any button click should clear key mapping selection.
            vd.vc.displayStates.mainView.control.selectedKeyCallbackFunctionIndex = -1;
        }
    }
    static void scrollCallback(GLFWwindow* window, double xoffset, double yoffset) {
        auto& vd = *static_cast< VisualDisplay* >(glfwGetWindowUserPointer(window));
        auto& proj = vd.vc.displaySettings.mainView.camera.projection;
        auto& fov = proj.fov;
        auto& scale = proj.scale;

        auto& io = ImGui::GetIO();
        if(!io.WantCaptureMouse) {
            fov -= 0.02 * yoffset;
            fov = std::clamp(fov, 0.01f, 3.00f);

            scale -= 0.05 * yoffset;
            scale = std::clamp(scale, 0.05f, 20.0f);
        }
    }
    static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
        auto& vd = *static_cast< VisualDisplay* >(glfwGetWindowUserPointer(window));

        KeyComb keyComb {
            key,
            bool(mods & GLFW_MOD_CONTROL),
            bool(mods & GLFW_MOD_ALT),
            bool(mods & GLFW_MOD_SHIFT),
        };
        auto& keyMapping = vd.vc.displaySettings.mainView.control.keyMapping;
        auto& selectedKeyFuncIndex = vd.vc.displayStates.mainView.control.selectedKeyCallbackFunctionIndex;

        if(action == GLFW_PRESS) {

            // Only perform the actual function outside key mapping mode.
            if(selectedKeyFuncIndex == -1) {
                const auto& keyCallbackMap = keyMapping.keyCallbackMap;
                auto it = keyCallbackMap.find(keyComb);
                if(it != keyCallbackMap.end()) {
                    // Invoke key callback
                    keyCallbackList[underlying(it->second)](&vd.vc.displaySettings, &vd.vc.displayStates);
                }
            }
        }

        if(action == GLFW_RELEASE) {
            // If in key mapping mode, update key mapping and clear selection.
            if(selectedKeyFuncIndex != -1) {
                keyMapping.update(keyComb, static_cast<KeyCallbackFunction>(selectedKeyFuncIndex));
                selectedKeyFuncIndex = -1;
            }
        }

    }


    // The overall opengl context. Must be at top
    VisualContext vc;

    // The ImGui context. Must be below the opengl context
    ImguiGuard guard;

    // Shader objects
    Shader shaderSurface { shader::VertexElementLight, shader::FragElementLight };
    Shader shaderLine    { shader::VertexElementLine,  shader::FragElementLine  };


    VisualDisplay(DisplayMode displayMode = DisplayMode::realtime) :
        vc(framebufferSizeCallback, cursorPosCallback, mouseButtonCallback, scrollCallback, keyCallback),
        guard(vc.window())
    {
        // Configure global opengl state
        glEnable(GL_DEPTH_TEST);

        // Configure GLFW window user.
        glfwSetWindowUserPointer(vc.window(), this);

        // Set initial display mode
        vc.displaySettings.displayMode = displayMode;

        // Setup realtime display profiles
        vc.displayStates.realtimeDataStates.profileData = makeDefaultElementProfileData();

    } // VisualDisplay()

    void run() {

        while (!glfwWindowShouldClose(vc.window())) {
            // input
            vc.processInput();

            // check for update for mesh data
            if(vc.displaySettings.displayMode == DisplayMode::trajectory) {
                appendTrajectoryDataFromLoadDock(vc.displayStates);
                playbackUpdateMaxFrames(vc.displayStates);

                // check whether updates are needed by new frame.
                playbackCheckTrajectory(vc.displaySettings, vc.displayStates, glfwGetTime());

                // update mesh data if needed.
                updateMeshDataForAllTrajectories(vc.displayStates.trajectoryDataStates, vc.displayStates.playback.currentFrame);
            } else {
                // update mesh data if needed.
                updateRealtimeMeshData(vc.displayStates.realtimeDataStates.profileData, sdfv);
            }

            auto& mainViewSettings = vc.displaySettings.mainView;
            auto& mainViewStates   = vc.displayStates.mainView;

            // select frame buffer: on screen display / offscreen snapshot
            const bool offscreen =
                mainViewStates.control.snapshotRenderingNextFrame ||
                vc.displayStates.playback.offscreenRender.has_value();
            mainViewStates.control.snapshotRenderingNextFrame = false;

            auto width  = mainViewSettings.canvas.width;
            auto height = mainViewSettings.canvas.height;
            auto projectionWidth  = mainViewSettings.canvas.width;
            auto projectionHeight = mainViewSettings.canvas.height;
            if(offscreen) {
                std::tie(width, height) = mainViewSettings.control.snapshotSize(width, height);
                if(
                    mainViewSettings.camera.projection.type == ObjectViewSettings::Camera::Projection::Type::orthographic &&
                    mainViewSettings.control.snapshotResolution == ObjectViewSettings::Control::SnapshotResolution::scaleWithScreen &&
                    mainViewSettings.control.snapshotUndoScaleOrtho
                ) {
                    // The projection size stay as is.
                    // so do nothing
                } else {
                    projectionWidth = width;
                    projectionHeight = height;
                }
            }

            VisualContext::OffscreenBufferGuard offscreenBufferGuard;
            if(offscreen) {
                offscreenBufferGuard.build(width, height);
                LOG(STEP) << "Rendering offscreen...";
            }
            glBindFramebuffer(GL_FRAMEBUFFER, offscreen ? offscreenBufferGuard.offscreenFbo : 0);

            glViewport(0, 0, width, height);
            {
                const auto& c = mainViewSettings.canvas.bgColor;
                glClearColor(c[0], c[1], c[2], c[3]);
            }
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            // transform
            const auto model          = glm::mat4(1.0f);
            const auto modelInvTrans3 = glm::transpose(glm::inverse(model));
            const auto view           = mainViewSettings.camera.view.view();
            const auto projection     = mainViewSettings.camera.projection.proj(projectionWidth, projectionHeight);

            // Start to render surface
            glUseProgram(shaderSurface.id());

            shaderSurface.setMat4("projection",     projection);
            shaderSurface.setMat4("model",          model);
            shaderSurface.setMat3("modelInvTrans3", modelInvTrans3);
            shaderSurface.setMat4("view",           view);

            shaderSurface.setBool("depthCueing",    mainViewSettings.canvas.depthCueing);
            shaderSurface.setFloat("depthCueingMin", mainViewSettings.canvas.depthCueingDepthMin);
            shaderSurface.setFloat("depthCueingMax", mainViewSettings.canvas.depthCueingDepthMax);
            shaderSurface.setVec3("CameraPos",      mainViewSettings.camera.view.position());

            // Set lights.
            for(int li = 0; li < ObjectViewSettings::Lighting::numDirLights; ++li) {
                constexpr int sz = 32;
                char name[sz];
                std::snprintf(name, sz, "dirLights[%d].direction", li);
                shaderSurface.setVec3(name, mainViewSettings.lighting.dirLights[li].direction);
                std::snprintf(name, sz, "dirLights[%d].ambient", li);
                shaderSurface.setVec3(name, mainViewSettings.lighting.dirLights[li].ambient);
                std::snprintf(name, sz, "dirLights[%d].diffuse", li);
                shaderSurface.setVec3(name, mainViewSettings.lighting.dirLights[li].diffuse);
                std::snprintf(name, sz, "dirLights[%d].specular", li);
                shaderSurface.setVec3(name, mainViewSettings.lighting.dirLights[li].specular);
            }
            for(int li = 0; li < ObjectViewSettings::Lighting::numPointLights; ++li) {
                constexpr int sz = 32;
                char name[sz];
                std::snprintf(name, sz, "pointLights[%d].position", li);
                shaderSurface.setVec3(name, mainViewSettings.lighting.pointLights[li].position);
                std::snprintf(name, sz, "pointLights[%d].ambient", li);
                shaderSurface.setVec3(name, mainViewSettings.lighting.pointLights[li].ambient);
                std::snprintf(name, sz, "pointLights[%d].diffuse", li);
                shaderSurface.setVec3(name, mainViewSettings.lighting.pointLights[li].diffuse);
                std::snprintf(name, sz, "pointLights[%d].specular", li);
                shaderSurface.setVec3(name, mainViewSettings.lighting.pointLights[li].specular);
                std::snprintf(name, sz, "pointLights[%d].constant", li);
                shaderSurface.setFloat(name, mainViewSettings.lighting.pointLights[li].constant);
                std::snprintf(name, sz, "pointLights[%d].linear", li);
                shaderSurface.setFloat(name, mainViewSettings.lighting.pointLights[li].linear);
                std::snprintf(name, sz, "pointLights[%d].quadratic", li);
                shaderSurface.setFloat(name, mainViewSettings.lighting.pointLights[li].quadratic);
            }


            const auto drawSurfaceProfileData = [&](auto& eachProfileData) {
                if(displayGeometryType(eachProfileData.profile) == DisplayGeometryType::surface) {

                    std::visit(
                        [&](const auto& profile) {
                            draw(eachProfileData.data, std::nullopt, shaderSurface, profile.displaySettings);
                        },
                        eachProfileData.profile
                    );
                }
            };
            if(vc.displaySettings.displayMode == DisplayMode::trajectory) {
                for(auto& traj : vc.displayStates.trajectoryDataStates.trajectories) {
                    if(traj.displayMasterSwitch) {
                        for(auto& profileData : traj.profileData) {
                            drawSurfaceProfileData(profileData);
                        }
                    }
                }
            }
            else {
                for(auto& profileData : vc.displayStates.realtimeDataStates.profileData) {
                    drawSurfaceProfileData(profileData);
                }
            }


            // start to render line
            glUseProgram(shaderLine.id());

            shaderLine.setMat4("projection",     projection);
            shaderLine.setMat4("model",          model);
            shaderLine.setMat3("modelInvTrans3", modelInvTrans3);
            shaderLine.setMat4("view",           view);

            shaderLine.setBool("depthCueing",    mainViewSettings.canvas.depthCueing);
            shaderLine.setFloat("depthCueingMin", mainViewSettings.canvas.depthCueingDepthMin);
            shaderLine.setFloat("depthCueingMax", mainViewSettings.canvas.depthCueingDepthMax);

            const auto drawLineProfileData = [&](auto& eachProfileData) {
                if(displayGeometryType(eachProfileData.profile) == DisplayGeometryType::line) {

                    std::visit(
                        [&](const auto& profile) {
                            draw(eachProfileData.data, std::nullopt, shaderLine, profile.displaySettings);
                        },
                        eachProfileData.profile
                    );
                }
            };
            if(vc.displaySettings.displayMode == DisplayMode::trajectory) {
                for(auto& traj : vc.displayStates.trajectoryDataStates.trajectories) {
                    if(traj.displayMasterSwitch) {
                        for(auto& profileData : traj.profileData) {
                            drawLineProfileData(profileData);
                        }
                    }
                }
            }
            else {
                for(auto& profileData : vc.displayStates.realtimeDataStates.profileData) {
                    drawLineProfileData(profileData);
                }
            }


            glBindVertexArray(0);

            // Output offscreen render results
            if(offscreen) {
                std::vector< std::uint8_t > data(width * height * 4); // 8-bit color
                glReadBuffer(GL_COLOR_ATTACHMENT0);
                glReadPixels(0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, data.data());

                int frameIndex = -1;
                if(vc.displayStates.playback.offscreenRender.has_value()) {
                    frameIndex = vc.displayStates.playback.currentFrame;
                }
                const auto snapshotPngFile = mainViewSettings.control.snapshotPngFile(frameIndex);
                stbi_write_png(snapshotPngFile.string().c_str(), width, height, 4, data.data(), 4 * width);
                LOG(INFO) << "Snapshot saved to " << snapshotPngFile;
            }

            // Update GUI
            imguiLoopRender(vc.displaySettings, vc.displayStates);
            // Perform async tasks
            dispatchAnAsyncTask(vc.displayStates.sync);

            // check
            glfwSwapBuffers(vc.window());
            glfwPollEvents();

        } // End main loop

    } // void run() const
}; // VisualDisplay

} // namespace medyan::visual

#endif

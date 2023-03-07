#ifndef MEDYAN_Visual_DisplaySettings_hpp
#define MEDYAN_Visual_DisplaySettings_hpp

#include <array>
#include <cstdio>
#include <filesystem>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "Util/Io/Log.hpp"
#include "Visual/KeyMapping.hpp"

namespace medyan::visual {

//-------------------------------------
// Global display mode
//-------------------------------------

enum class DisplayMode {
    realtime,        // Display system data from memory during simulation
    trajectory,      // Read and display trajectory from output files

    // Reflection hack for counting elements
    last_
};

constexpr auto text(DisplayMode mode) {
    switch(mode) {
        case DisplayMode::realtime:   return "realtime";
        case DisplayMode::trajectory: return "trajectory";
        default:                      return "";
    }
}

//-------------------------------------
// GUI settings
//-------------------------------------

struct GuiSettings {
    // the master switch
    bool enabled = false;

    // windows switches
    bool helpWindow = false;
    bool profileWindow = true;
    bool inspectWindow = false;
};

//-------------------------------------
// Graphics settings
//-------------------------------------

struct ObjectViewSettings {

    struct Canvas {
        std::array< float, 4 > bgColor { 0.0f, 0.0f, 0.0f, 0.0f };
        bool  depthCueing = false;
        float depthCueingDepthMin = 0;
        float depthCueingDepthMax = 1;

        // The framebuffer size, not necessarily the window size.
        int width = 1200;
        int height = 800;
    };

    struct Lighting {
        struct DirLight {
            glm::vec3 direction { 0.0f, 0.0f, 1.0f };
            glm::vec3 ambient { 0.3f, 0.3f, 0.3f };
            glm::vec3 diffuse { 0.3f, 0.3f, 0.3f };
            glm::vec3 specular { 0.1f, 0.1f, 0.1f };
        };
        struct PointLight {
            glm::vec3 position {};
            glm::vec3 ambient { 0.0f, 0.0f, 0.0f };
            glm::vec3 diffuse { 0.4f, 0.4f, 0.4f };
            glm::vec3 specular { 0.8f, 0.8f, 0.8f };
            float constant = 1.0f;
            float linear = 1.4e-4f;
            float quadratic = 7.2e-8f;
        };

        // Currently, the number of lights is fixed (check the shader source).
        static constexpr int numDirLights = 2;
        static constexpr int numPointLights = 4;

        std::array< DirLight,   numDirLights >   dirLights {{
            DirLight { {  1.0f,  1.0f,  1.0f }, },
            DirLight { { -1.0f, -1.0f, -1.0f }, },
        }};
        std::array< PointLight, numPointLights > pointLights {{
            PointLight { { -500.0f, -500.0f, -500.0f }, },
            PointLight { { -500.0f, 3500.0f, 3500.0f }, },
            PointLight { { 3500.0f, -500.0f, 3500.0f }, },
            PointLight { { 3500.0f, 3500.0f, -500.0f }, },
        }};
    };

    // The Camera struct stores
    // - The view transformation, containing the camera position and orientation, and
    // - The projection transformation, containing the field of view, aspect ratio, and near and far clipping planes.
    struct Camera {
        struct View {
            glm::mat4 view() const {
                auto trans1 = glm::translate(glm::mat4(1.0f), -target);
                auto rot = glm::mat4_cast(viewQuat);
                auto trans2 = glm::translate(glm::mat4(1.0f), glm::vec3(0, 0, -distance));
                return trans2 * rot * trans1;
            }
            glm::vec3 position() const {
                return backward() * distance + target;
            }
            glm::vec3 right() const {
                return glm::conjugate(viewQuat) * glm::vec3(1, 0, 0);
            }
            glm::vec3 up() const {
                return glm::conjugate(viewQuat) * glm::vec3(0, 1, 0);
            }
            glm::vec3 backward() const {
                return glm::conjugate(viewQuat) * glm::vec3(0, 0, 1);
            }

            // Positions
            glm::quat viewQuat = glm::quat(1.0f, 0.0f, 0.0f, 0.0f); // Should always be normalized.
            float     distance = 2000.f;
            glm::vec3 target   = glm::vec3(0.0f, 0.0f, 0.0f);
        };

        struct Projection {
            enum class Type { orthographic, perspective, last_ };

            Type  type       = Type::perspective;
            float fov        = glm::radians(45.0f); // perspective, vertical field of view
            float scale      = 1.0f;                // orthographic
            float zNear      = 10.0f;
            float zFar       = 5000.0f;

            // Helper functions

            auto projOrtho(float width, float height) const {
                return glm::ortho(-width * 0.5f * scale, width * 0.5f * scale, -height * 0.5f * scale, height * 0.5f * scale, zNear, zFar);
            }
            auto projPerspective(float width, float height) const {
                return glm::perspective(fov, height == 0 ? 1.0f : width / height, zNear, zFar);
            }
            auto proj(float width, float height) const {
                return type == Type::orthographic ?
                    projOrtho      (width, height):
                    projPerspective(width, height);
            }
        };

        View view;
        Projection projection;
    };

    struct Control {
        enum class CameraMouseMode {
            none,
            rotate,
            pan,
            last_
        };
        enum class SnapshotResolution {
            scaleWithScreen,
            manual,
            last_
        };

        // Key mapping
        //-----------------------------
        KeyMapping keyMapping;

        // camera
        //-----------------------------
        float cameraKeyPositionPerFrame = 500.0f;

        CameraMouseMode cameraMouseMode = CameraMouseMode::rotate;
        // When rotating camera, the camera-target distance divided by the following value gives the radius of the rotation sphere.
        // When the mouse moves a little, it should appear that the center of front surface of the rotation sphere is dragged by the mouse.
        // Let k be the following divisor, then the angle per pixel (app) can be computed as:
        // - perspective:   app = (k - 1) * vertical_fov / height_in_pixel
        // - orthographic:  app = scale * (k / camera_target_distance)
        // k should be in range (1, +infinity), where 1 indicates smallest rotation speed (zero if perspective) and +infinity indicates infinite rotation speed.
        float cameraRotateSphereDistanceDivisor = 2.5f;
        float cameraPanPositionPerCursorPixel = 2.0f;      // camera move distance per pixel (pan)

        // snapshot saving section
        //-----------------------------
        std::filesystem::path snapshotSavePath = ".";
        std::string           snapshotBaseName = "snapshot";
        SnapshotResolution    snapshotResolution = SnapshotResolution::scaleWithScreen;
        float                 snapshotScale = 1.0f;   // Used with "scale with screen"
        // Used with "scale with screen".
        // If set to true, it will undo the scaling when calculating orthographic projection.
        bool                  snapshotUndoScaleOrtho = true;
        int                   snapshotWidth = 1920;   // Used with "manual"
        int                   snapshotHeight = 1080;   // Used with "manual"

        // This function gives the snapshot size in integers (width, height).
        auto snapshotSize(int screenViewWidth, int screenViewHeight) const {
            int width, height;
            if(snapshotResolution == SnapshotResolution::scaleWithScreen) {
                width = screenViewWidth * snapshotScale;
                height = screenViewHeight * snapshotScale;
            } else {
                width = snapshotWidth;
                height = snapshotHeight;
            }

            return std::pair { width, height };
        }

        // Finds the snapshot file name.
        //
        // The file name will be displayed as "baseName-index.ext", where the
        // index will be displayed with padding zeros.
        //
        // The index should be -1 or non-negative, where -1 will not append
        // "-index" after the base name.
        //
        // One can use
        //   ffmpeg -framerate {frameRate} -i {baseName}-%06d.png -c:v libx264 -vf "scale=ceil(iw/2)*2:ceil(ih/2)*2,format=yuv420p" -crf 18 {output}.mp4
        // to stitch the images to make a video.
        auto snapshotPngFile(int index = -1) const {
            if(index == -1) {
                return snapshotSavePath / (snapshotBaseName + ".png");
            }
            else {
                if(index < 0 || index > 999999) {
                    log::error("Snapshot index {} out of range.", index);
                    throw std::runtime_error("Snapshot index out of range");
                }
                // 11 chars after base name: "-012345.png"
                char appendStr[12] {};
                std::snprintf(appendStr, 12, "-%06d.png", index);
                return snapshotSavePath / (snapshotBaseName + appendStr);
            }
        }
    };

    // Canvas
    Canvas canvas;

    Lighting lighting;

    // View
    Camera camera;

    // Control
    Control control;
};

constexpr auto text(ObjectViewSettings::Camera::Projection::Type proj) {
    switch(proj) {
        case ObjectViewSettings::Camera::Projection::Type::orthographic: return "orthographic";
        case ObjectViewSettings::Camera::Projection::Type::perspective:  return "perspective";
        default:                                                 return "";
    }
}
constexpr auto text(ObjectViewSettings::Control::CameraMouseMode mode) {
    switch(mode) {
        case ObjectViewSettings::Control::CameraMouseMode::none:   return "none";
        case ObjectViewSettings::Control::CameraMouseMode::rotate: return "rotate";
        case ObjectViewSettings::Control::CameraMouseMode::pan:    return "pan";
        default:                                                   return "";
    }
}
constexpr auto text(ObjectViewSettings::Control::SnapshotResolution value) {
    switch(value) {
        case ObjectViewSettings::Control::SnapshotResolution::scaleWithScreen: return "scale with screen";
        case ObjectViewSettings::Control::SnapshotResolution::manual:          return "manual";
        default:                                                               return "";
    }
}

//-------------------------------------
// Playback settings
//-------------------------------------

struct TrajectoryPlaybackSettings {
    float fps = 6;

};

//-------------------------------------
// All display settings
//-------------------------------------

struct DisplaySettings {
    DisplayMode displayMode = DisplayMode::realtime;

    // GUI
    GuiSettings gui;

    // Main display
    ObjectViewSettings mainView;

    // Trajectory playback
    TrajectoryPlaybackSettings playback;
};

} // namespace medyan::visual

#endif

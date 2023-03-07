#ifndef MEDYAN_Visual_Control_hpp
#define MEDYAN_Visual_Control_hpp

#include "Visual/DisplaySettings.hpp"
#include "Visual/DisplayStates.hpp"

namespace medyan::visual {

// mouse camera drag

inline void mouseDragRotateCameraAroundTarget(
    ObjectViewSettings::Camera::View& view,
    double mouseLastX,
    double mouseLastY,
    double mouseX,
    double mouseY,
    float  cameraRotateAnglePerCursorPixel
) {

    const float dx = mouseX - mouseLastX;
    const float dy = mouseLastY - mouseY;
    const float d = std::sqrt(dx * dx + dy * dy);
    if(d > 0) {
        const float angle = d * cameraRotateAnglePerCursorPixel;
        const auto axis = glm::vec3(-dy / d, dx / d, 0.0f);
        view.viewQuat = glm::normalize(glm::angleAxis(angle, axis) * view.viewQuat);
    }
}

inline void mouseDragPanCamera(
    ObjectViewSettings::Camera::View& view,
    double mouseLastX,
    double mouseLastY,
    double mouseX,
    double mouseY,
    float  cameraPanPositionPerCursorPixel
) {
    const auto objectDelta = (
        view.right() * float(mouseX - mouseLastX) +
        view.up()    * float(mouseLastY - mouseY)
    ) * cameraPanPositionPerCursorPixel;

    view.target   -= objectDelta;
}

inline void mouseDragCameraView(
    ObjectViewSettings::Camera::View& view,
    double mouseLastX,
    double mouseLastY,
    double mouseX,
    double mouseY,
    const ObjectViewSettings::Control& control,
    const ObjectViewSettings::Camera::Projection& proj,
    const ObjectViewSettings::Canvas& canvas
) {
    using MM = ObjectViewSettings::Control::CameraMouseMode;

    switch(control.cameraMouseMode) {
        case MM::rotate:
            mouseDragRotateCameraAroundTarget(
                view,
                mouseLastX, mouseLastY, mouseX, mouseY,
                proj.type == ObjectViewSettings::Camera::Projection::Type::perspective
                    ? (control.cameraRotateSphereDistanceDivisor - 1) * proj.fov / canvas.height
                    : control.cameraRotateSphereDistanceDivisor * proj.scale / view.distance
            );
            break;

        case MM::pan:
            mouseDragPanCamera(
                view,
                mouseLastX, mouseLastY, mouseX, mouseY,
                control.cameraPanPositionPerCursorPixel
            );
            break;

        default:
            break;
    }
}


// Automatically adjust the camera parameters to fit the current scene.
inline void autoAdjustCamera(
    ObjectViewSettings::Camera& camera,
    const DisplaySettings& settings,
    const DisplayStates& states
) {
    auto& view = camera.view;
    auto& proj = camera.projection;

    // Find bounding box of all loaded mesh data, except for trajectories that are completely disabled.
    // Using the bounding box to obtain an approximation of the bounding sphere, because the latter is non-trivial.
    std::optional<std::array<float, 6>> bounds;

    if(settings.displayMode == DisplayMode::realtime) {
        for(auto& profileData : states.realtimeDataStates.profileData) {
            extend3DExtents(bounds, profileData.data);
        }
    }
    else {
        for(auto& traj : states.trajectoryDataStates.trajectories) {
            if(traj.displayMasterSwitch) {
                for(auto& profileData : traj.profileData) {
                    extend3DExtents(bounds, profileData.data);
                }
            }
        }
    }

    if(bounds.has_value()) {
        // Obtain an approximation of the bounding sphere.
        const auto& b = *bounds;
        const auto cx = (b[0] + b[1]) / 2;
        const auto cy = (b[2] + b[3]) / 2;
        const auto cz = (b[4] + b[5]) / 2;
        const auto radius = std::sqrt(
            (b[1] - b[0]) * (b[1] - b[0]) +
            (b[3] - b[2]) * (b[3] - b[2]) +
            (b[5] - b[4]) * (b[5] - b[4])
        ) / 2;

        // Adjust camera parameters.
        view.target = glm::vec3(cx, cy, cz);

        auto defaultProj = ObjectViewSettings::Camera::Projection {};
        proj.fov = defaultProj.fov;
        const float verticalRoomFactor = 1.05f; // To allow for extra room between vertical boundaries and the bounding sphere.
        view.distance = radius * verticalRoomFactor / std::sin(proj.fov / 2);

        const float boundRoomFactor = 1.05f; // To allow for extra room between boundaries and the bounding sphere.
        const auto minScreenSize = std::min(settings.mainView.canvas.width, settings.mainView.canvas.height);
        proj.scale = 2 * radius * boundRoomFactor / minScreenSize;

        const float frontBackRoomFactor = 1.2f; // To allow for extra room between front and back boundaries.
        proj.zNear = std::max(view.distance / 100, view.distance - radius * frontBackRoomFactor);
        proj.zFar  = view.distance + radius * frontBackRoomFactor;
    }
    else {
        // No visible objects.
        ObjectViewSettings::Camera::View defaultView {};
        view.target = defaultView.target;
        view.distance = defaultView.distance;

        ObjectViewSettings::Camera::Projection defaultProj {};
        proj.fov = defaultProj.fov;
        proj.scale = defaultProj.scale;
        proj.zNear = defaultProj.zNear;
        proj.zFar  = defaultProj.zFar;
    }
}

} // namespace medyan::visual

#endif

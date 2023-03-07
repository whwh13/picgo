#ifndef MEDYAN_Visual_Gui_GuiAux_hpp
#define MEDYAN_Visual_Gui_GuiAux_hpp

#include <algorithm> // clamp
#include <cstdio>
#include <filesystem>
#include <type_traits>

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <imgui_internal.h>
#include <imgui_stdlib.h>
#include <nfd.h>

#include "Util/Io/Log.hpp"
#include "Util/Math/ColorMap.hpp"

namespace medyan::visual {


//-----------------------------------------------------------------------------
// Auxiliary functions
//-----------------------------------------------------------------------------

// Returns whether the color is changed.
template<
    int colorDof,
    std::enable_if_t< colorDof == 3 || colorDof == 4 >* = nullptr
>
inline bool guiAuxColorPickerPopup(
    const char*        strId, 
    float*             pColor
) {
    bool changed = false;

    const bool bgColorPopup = ImGui::ColorButton(
        strId,
        ImVec4(
            pColor[0],
            pColor[1],
            pColor[2],
            colorDof == 3 ? 1.0f : pColor[3]
        ),
        0
    );
    ImGui::SameLine();
    ImGui::Text(strId);

    if(bgColorPopup) {
        ImGui::OpenPopup(strId);
    }

    if(ImGui::BeginPopup(strId)) {
        if(colorDof == 3) {
            changed = ImGui::ColorPicker3(strId, pColor, ImGuiColorEditFlags_PickerHueWheel);
        } else {
            changed = ImGui::ColorPicker4(strId, pColor, ImGuiColorEditFlags_PickerHueWheel);
        }

        ImGui::EndPopup();
    }

    return changed;
}

inline bool guiAuxColorPicker3Popup(const char* strId, float* pColor) {
    return guiAuxColorPickerPopup<3>(strId, pColor);
}
inline bool guiAuxColorPicker4Popup(const char* strId, float* pColor) {
    return guiAuxColorPickerPopup<4>(strId, pColor);
}

// Function to build combo box automatically for an enumeration type.
//
// Returns whether a new value is selected.
//
// Notes:
//   - The elements in the enum type must be automatically valued (ie the
//     values are automatically 0, 1, 2, ...).
//   - The type must have "last_" as the last element.
//   - A "text(Enum)" function must be implemented to display the elements.
template<
    typename Enum,
    typename Reselect,              // function void(Enum old, Enum new) to execute when selected
    std::enable_if_t<
        std::is_enum_v< Enum > &&
        std::is_invocable_r_v< void, Reselect, Enum, Enum > // Reselect: (Enum, Enum) -> void
    >* = nullptr                    // type requirements
>
inline bool guiAuxEnumComboBox(
    const char*          name,
    Enum&                value,
    Reselect&&           reselect,
    ImGuiSelectableFlags flags = 0
) {
    bool changed = false;

    if(ImGui::BeginCombo(name, text(value), 0)) {
        for (int i = 0; i < underlying(Enum::last_); ++i) {

            const Enum valueI = static_cast<Enum>(i);

            const bool isSelected = (value == valueI);
            if (ImGui::Selectable(text(valueI), isSelected, flags)) {
                const auto oldValue = value;
                value = valueI;
                reselect(oldValue, valueI);

                if(!isSelected) {
                    // selected one item that was previously not selected
                    changed = true;
                }
            }

            // Set the initial focus when opening the combo (scrolling + keyboard navigation focus)
            if (isSelected) {
                ImGui::SetItemDefaultFocus();
            }
        }
        ImGui::EndCombo();
    }

    return changed;
}
// Where Reselect function is a no-op.
template<
    typename Enum,
    std::enable_if_t< std::is_enum_v< Enum > >* = nullptr     // type requirements
>
inline bool guiAuxEnumComboBox(
    const char*          name,
    Enum&                value,
    ImGuiSelectableFlags flags = 0
) {
    return guiAuxEnumComboBox(name, value, [](Enum, Enum) {}, flags);
}

// Tristate checkbox
//
// State: 0 -> unchecked, 1 -> checked, -1 -> intermediate
// Clicking the checkbox will never set the state to -1.
//
// adapted from https://github.com/ocornut/imgui/issues/2644#issuecomment-507023896
inline bool guiAuxCheckboxTristate(const char* label, int* pState) {
    bool clicked = false;
    if (*pState == -1)
    {
        ImGui::PushItemFlag(ImGuiItemFlags_MixedValue, true);
        bool b = false;
        clicked = ImGui::Checkbox(label, &b);
        if (clicked)
            *pState = 1;
        ImGui::PopItemFlag();
    }
    else
    {
        bool b = (*pState != 0);
        clicked = ImGui::Checkbox(label, &b);
        if (clicked)
            *pState = (int)b;
    }
    return clicked;
}

// check box to toggle flags
//
// Returns whether the flag is changed.
template< typename UInt >
inline bool guiAuxCheckboxFlags(const char* label, UInt& flag, UInt setTo) {
    const UInt oldFlag     = flag;
    const bool atLeast1Set = flag & setTo;
    const bool allSet      = (flag & setTo) == setTo;
    int        state       = allSet ? 1 : (atLeast1Set ? -1 : 0);

    if(guiAuxCheckboxTristate(label, &state)) {
        if(state == 1) {
            flag |= setTo;
        } else if(state == 0) {
            flag &= ~setTo;
        } // cannot be -1
    }

    return flag != oldFlag;
}

// select file
inline void guiAuxSelectFile(std::filesystem::path& file, const std::filesystem::path& defaultPath, bool selectDirectory = false) {
    nfdchar_t* filepath = nullptr;
    auto result = selectDirectory
        ? NFD_PickFolder(defaultPath.string().c_str(), &filepath)
        : NFD_OpenDialog(nullptr, defaultPath.string().c_str(), &filepath);

    if ( result == NFD_OKAY ) {
        file = filepath;
        std::free(filepath);
    }
    else if ( result == NFD_CANCEL ) {
        // do nothing
    }
    else {
        log::error("Error: {}", NFD_GetError());
    }
}

// Draw horizontal color bar.
inline void guiAuxDrawHorizontalColorBar(const char* label, const colormap::DynamicColorMap<float>& colorMap) {
    // Auxiliary function to convert ColorRgb to ImU32.
    const auto toImU32 = [](const colormap::ColorRgb<float>& color) {
        return ImGui::GetColorU32(IM_COL32(
            std::clamp(static_cast<int>(color[0] * 256), 0, 255),
            std::clamp(static_cast<int>(color[1] * 256), 0, 255),
            std::clamp(static_cast<int>(color[2] * 256), 0, 255),
            255
        ));
    };

    ImGui::PushItemWidth(-ImGui::GetFontSize() * 10);
    ImDrawList* drawList = ImGui::GetWindowDrawList();

    // Draw gradients.
    const ImVec2 barSize = ImVec2(ImGui::CalcItemWidth(), ImGui::GetFrameHeight());
    const ImVec2 p0 = ImGui::GetCursorScreenPos();
    const Size numAnchors = colorMap.interpList.size();

    if(numAnchors > 1) {
        ImU32 colorLast = toImU32(colorMap.interpList[0]);
        for(int i = 1; i < numAnchors; ++i) {
            const ImVec2 p1 = ImVec2(p0.x + barSize.x * (i-1) / (numAnchors-1), p0.y);
            const ImVec2 p2 = ImVec2(p0.x + barSize.x * i / (numAnchors-1), p0.y + barSize.y);
            const ImU32 colorNew = toImU32(colorMap.interpList[i]);
            drawList->AddRectFilledMultiColor(p1, p2, colorLast, colorNew, colorNew, colorLast);

            colorLast = colorNew;
        }
    }

    char buf[64] {};
    std::snprintf(buf, sizeof(buf), "##%s", label);
    ImGui::InvisibleButton(buf, barSize);
    ImGui::PopItemWidth();
}

} // namespace medyan::visual

#endif

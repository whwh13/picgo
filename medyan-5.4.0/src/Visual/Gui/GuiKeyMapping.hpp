#ifndef MEDYAN_Visual_Gui_GuiKeyMapping_hpp
#define MEDYAN_Visual_Gui_GuiKeyMapping_hpp

#include <iterator> // std::size

#include <imgui.h>

#include "Visual/DisplaySettings.hpp"
#include "Visual/DisplayStates.hpp"
#include "Visual/Gui/GuiAux.hpp"

namespace medyan::visual {

inline void guiComponentKeyMapping(KeyMapping& keyMapping, ObjectViewStates::Control& controlStates) {

    if(ImGui::Button("Reset default")) {
        keyMapping = KeyMapping::createDefault();
    }

    if (ImGui::BeginTable("split", 2, ImGuiTableFlags_BordersOuter | ImGuiTableFlags_Resizable))
    {

        // Iterate placeholder objects (all the same data)
        for (int i = 0; i < underlying(KeyCallbackFunction::last_); ++i) {

            auto keyfunci = static_cast<KeyCallbackFunction>(i);
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::AlignTextToFramePadding();
            ImGui::Text(text(keyfunci));

            ImGui::TableSetColumnIndex(1);
            std::string selectableLabel = display(keyMapping.keyBindingList[i]) + "##" + std::to_string(i);
            if(ImGui::Selectable(selectableLabel.c_str(), controlStates.selectedKeyCallbackFunctionIndex == i)) {
                controlStates.selectedKeyCallbackFunctionIndex = i;
            }

        }
        ImGui::EndTable();
    }

}

} // namespace medyan::visual

#endif

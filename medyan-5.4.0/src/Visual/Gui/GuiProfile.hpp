#ifndef MEDYAN_Visual_Gui_GuiProfile_hpp
#define MEDYAN_Visual_Gui_GuiProfile_hpp

#include "Visual/DisplaySettings.hpp"
#include "Visual/DisplayStates.hpp"
#include "Visual/Gui/GuiAux.hpp"

namespace medyan::visual {

// Functions for profile configuration
// Returns whether the profile is changed.
inline bool guiGeometryDisplaySettings(SurfaceDisplaySettings& settings) {
    bool changed = false;

    changed |= ImGui::Checkbox("enabled", &settings.enabled);
    changed |= guiAuxEnumComboBox("polygon", settings.polygonMode);

    return changed;
}
inline bool guiGeometryDisplaySettings(LineDisplaySettings& settings) {
    bool changed = false;

    changed |= ImGui::Checkbox("enabled", &settings.enabled);

    return changed;
}


// Frame index -1 indicates using the current frame with an unknown frame index.
inline bool guiActiveProfileConfig(MembraneProfile& profile, const DisplayTypeMap& displayTypeMap, Index frame = -1) {
    bool changed = false;

    changed |= guiAuxEnumComboBox("element mode", profile.displaySettings.elementMode);

    if(profile.displaySettings.elementMode == MembraneDisplaySettings::ElementMode::edge) {
        changed |= ImGui::SliderFloat("edge extrude radius", &profile.displaySettings.edgeExtrudeRadius, 1.0f, 10.0f, "%.1f");
        changed |= ImGui::SliderInt("edge extrude sides", &profile.displaySettings.edgeExtrudeSides, 4, 20);
    }

    changed |= guiGeometryDisplaySettings(profile.displaySettings.surface);

    changed |= guiAuxEnumComboBox("color mode", profile.displaySettings.colorMode);
    if(profile.displaySettings.colorMode == MembraneDisplaySettings::ColorMode::attribute) {
        auto& attrIndex = profile.displaySettings.attributeIndex;
        // Clamp attribute index.
        if(attrIndex < -1 || attrIndex >= displayTypeMap.membraneVertexAttributeNames.size()) {
            attrIndex = -1;
        }
        const auto& attrName = [&](Index index) -> const char* {
            if(index == -1) return "none";
            return displayTypeMap.membraneVertexAttributeNames[index].c_str();
        };

        // Select attribute name.
        if(ImGui::BeginCombo("attribute", attrName(attrIndex), 0)) {
            for(Index i = -1; i < (Size)displayTypeMap.membraneVertexAttributeNames.size(); ++i) {

                const bool isSelected = attrIndex == i;
                if(ImGui::Selectable(attrName(i), isSelected, 0)) {
                    if(!isSelected) {
                        // A new attribute is chosen.
                        attrIndex = i;

                        // Because we selected an attribute that was not previously selected.
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

        // Attribute value ranges.
        changed |= ImGui::Checkbox("auto attrib range", &profile.displaySettings.autoAttributeRange);
        if(frame < (Size)profile.displayStates.attribRanges.size()) {
            ImGui::Text("attrib range (%.3g, %.3g)",
                profile.displayStates.attribRanges[frame == -1 ? 0 : frame].min,
                profile.displayStates.attribRanges[frame == -1 ? 0 : frame].max
            );
        } else {
            ImGui::Text("attrib range (-, -)");
        }

        if(!profile.displaySettings.autoAttributeRange) {
            const bool attribMinChanged = ImGui::DragFloat("min attrib", &profile.displaySettings.manualAttribMin, 0.01f);
            if(attribMinChanged) {
                changed = true;
                profile.displaySettings.manualAttribMax = std::max(profile.displaySettings.manualAttribMin, profile.displaySettings.manualAttribMax);
            }
            const bool attribMaxChanged = ImGui::DragFloat("max attrib", &profile.displaySettings.manualAttribMax, 0.01f);
            if(attribMaxChanged) {
                changed = true;
                profile.displaySettings.manualAttribMin = std::min(profile.displaySettings.manualAttribMin, profile.displaySettings.manualAttribMax);
            }
        }

        // Choose color map.
        guiAuxDrawHorizontalColorBar("colorbar display", profile.displaySettings.colorMap);
        ImGui::SameLine();
        if(ImGui::BeginCombo("##colormap", "colormap", 0)) {
            struct DynamicColorMapInfo {
                const char* name;
                const colormap::DynamicColorMap<float>* pmap;
            };
            const DynamicColorMapInfo colorMapInfo[] {
                { "bwr", &colormap::bwrf },
                { "jet", &colormap::jetf },
            };

            for(Index i = 0; i < std::size(colorMapInfo); ++i) {
                if(ImGui::Selectable(colorMapInfo[i].name, false, 0)) {
                    profile.displaySettings.colorMap = *colorMapInfo[i].pmap;
                    changed = true;
                }
            }
            ImGui::EndCombo();
        }
    } else {
        changed |= guiAuxColorPicker3Popup("fixed color", profile.displaySettings.colorFixed.value.data());
    }

    return changed;
}
inline bool guiActiveProfileConfig(FilamentProfile& profile, const DisplayTypeMap& displayTypeMap, Index frame = -1) {
    bool changed = false;

    // Selector.
    // type
    {
        // on/off checkbox
        bool filterOn = profile.selector.type.has_value();
        const bool filterOnChanged = ImGui::Checkbox("filter type##filter type check", &filterOn);
        changed |= filterOnChanged;

        if(filterOnChanged) {
            if(filterOn) profile.selector.type.emplace();
            else         profile.selector.type.reset();
        }

        // show text box
        if(filterOn) {
            changed |= ImGui::InputInt("filter type##filter type text", &*profile.selector.type);
        }
    }

    // path mode specific
    changed |= guiAuxEnumComboBox("path mode", profile.displaySettings.pathMode);

    switch(profile.displaySettings.pathMode) {
        case FilamentDisplaySettings::PathMode::line:
            break;
        case FilamentDisplaySettings::PathMode::extrude:
            changed |= ImGui::SliderFloat("extrude radius", &profile.displaySettings.pathExtrudeRadius, 1.0f, 24.0f, "%.1f");
            changed |= ImGui::SliderInt("extrude sides", &profile.displaySettings.pathExtrudeSides, 4, 20);
            break;
        case FilamentDisplaySettings::PathMode::bead:
            changed |= ImGui::SliderFloat("bead radius", &profile.displaySettings.beadRadius, 1.0f, 50.0f, "%.1f");
            break;
    }

    if(displayGeometryType(profile.displaySettings) == DisplayGeometryType::line) {
        changed |= guiGeometryDisplaySettings(profile.displaySettings.line);
    } else {
        changed |= guiGeometryDisplaySettings(profile.displaySettings.surface);
    }

    changed |= guiAuxColorPicker3Popup("fixed color", profile.displaySettings.colorFixed.value.data());

    return changed;
}
inline bool guiActiveProfileConfig(LinkerProfile& profile, const DisplayTypeMap& displayTypeMap, Index frame = -1) {
    bool changed = false;

    // Selectors.
    // name
    {
        // on/off checkbox
        bool filterOn = profile.selector.name.has_value();
        const bool filterOnChanged = ImGui::Checkbox("filter name##filter name check", &filterOn);
        changed |= filterOnChanged;

        if(filterOnChanged) {
            if(filterOn) profile.selector.name.emplace();
            else         profile.selector.name.reset();
        }

        // show text box
        if(filterOn) {
            changed |= ImGui::InputText("filter name##filter name text", &*profile.selector.name);
        }
    }
    // subtype
    {
        // on/off checkbox
        bool filterOn = profile.selector.subtype.has_value();
        const bool filterOnChanged = ImGui::Checkbox("filter subtype##filter subtype check", &filterOn);
        changed |= filterOnChanged;

        if(filterOnChanged) {
            if(filterOn) profile.selector.subtype.emplace();
            else         profile.selector.subtype.reset();
        }

        // show text box
        if(filterOn) {
            changed |= ImGui::InputInt("filter subtype##filter subtype text", &*profile.selector.subtype);
        }
    }

    // path mode specific
    changed |= guiAuxEnumComboBox("path mode", profile.displaySettings.pathMode);

    switch(profile.displaySettings.pathMode) {
        case LinkerDisplaySettings::PathMode::line:
            break;
        case LinkerDisplaySettings::PathMode::extrude:
            changed |= ImGui::SliderFloat("extrude radius", &profile.displaySettings.pathExtrudeRadius, 1.0f, 24.0f, "%.1f");
            changed |= ImGui::SliderInt("extrude sides", &profile.displaySettings.pathExtrudeSides, 4, 20);
            break;
    }

    if(displayGeometryType(profile.displaySettings) == DisplayGeometryType::line) {
        changed |= guiGeometryDisplaySettings(profile.displaySettings.line);
    } else {
        changed |= guiGeometryDisplaySettings(profile.displaySettings.surface);
    }

    // color
    changed |= guiAuxColorPicker3Popup("fixed color", profile.displaySettings.colorFixed.value.data());

    return changed;
}
inline bool guiActiveProfileConfig(BubbleProfile& profile, const DisplayTypeMap& displayTypeMap, Index frame = -1) {
    bool changed = false;

    changed |= ImGui::SliderInt("segs longitude", &profile.displaySettings.sphereLongitudeSegs, 6, 30);
    changed |= ImGui::SliderInt("segs latitude", &profile.displaySettings.sphereLatitudeSegs, 3, 15);

    changed |= guiGeometryDisplaySettings(profile.displaySettings.surface);

    changed |= guiAuxColorPicker3Popup("fixed color", profile.displaySettings.colorFixed.value.data());

    return changed;
}
inline bool guiActiveProfileConfig(AuxLineProfile& profile, const DisplayTypeMap& displayTypeMap, Index frame = -1) {
    bool changed = false;

    changed |= guiAuxCheckboxFlags("compartment border", profile.displaySettings.flag, AuxLineDisplaySettings::targetCompartmentBorder);
    changed |= guiAuxCheckboxFlags("compartment grid", profile.displaySettings.flag, AuxLineDisplaySettings::targetCompartmentAll);

    changed |= guiGeometryDisplaySettings(profile.displaySettings.line);

    return changed;
}

// Returns whether the profile is changed
inline bool guiActiveProfileConfig(ElementProfile& profile, const DisplayTypeMap& displayTypeMap, Index frame = -1) {
    bool changed = false;

    // Combo box to select target
    if(ImGui::BeginCombo("target", elementProfileDisplayName(profile.index()), 0)) {
        for (int i = 0; i < std::variant_size_v< ElementProfile >; ++i) {

            const bool isSelected = (profile.index() == i);
            if (ImGui::Selectable(elementProfileDisplayName(i), isSelected, 0)) {
                if(!isSelected) {
                    // A new profile is chosen
                    profile = makeProfileWithIndex(i);

                    // Because we selected a profile that was not previously selected
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

    // Display actual settings
    changed |= std::visit([&](auto& actualProfile) { return guiActiveProfileConfig(actualProfile, displayTypeMap, frame); }, profile);

    return changed;
}

} // namespace medyan::visual

#endif

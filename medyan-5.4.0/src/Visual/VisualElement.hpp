#ifndef MEDYAN_Visual_VisualElement_Hpp
#define MEDYAN_Visual_VisualElement_Hpp

#include <cstdint>
#include <memory> // shared_ptr
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "utility.h" // Overload
#include "Visual/Common.hpp"
#include "Visual/MeshData.hpp"
#include "Visual/Shader.hpp"

namespace medyan::visual {

//-------------------------------------
// Element selectors and profiles
//-------------------------------------

struct MembraneSelector {
};
struct MembraneProfile {
    // selects all membranes
    MembraneSelector selector;

    // apply this setting
    MembraneDisplaySettings displaySettings;

    // States.
    MembraneDisplayStates displayStates;
};
inline auto select(
    const std::vector< MembraneFrame >& membranes,
    const MembraneSelector& selector
) {
    std::vector< const MembraneFrame* > res;
    for(auto& membrane : membranes) {
        res.push_back(&membrane);
    }
    return res;
}

struct FilamentSelector {
    std::optional< int > type;
};
struct FilamentProfile {
    // selects all filaments
    FilamentSelector selector;

    // apply this setting
    FilamentDisplaySettings displaySettings;

};
inline auto select(
    const std::vector< FilamentFrame >& filaments,
    const FilamentSelector& selector
) {
    std::vector< const FilamentFrame* > res;
    for(auto& filament : filaments) {
        if(!selector.type.has_value() || *selector.type == filament.type) {
            res.push_back(&filament);
        }
    }
    return res;
}


struct LinkerSelector {
    std::optional< std::string > name;
    std::optional< int > subtype;
};
struct LinkerProfile {
    // selector
    LinkerSelector selector;

    // setting
    LinkerDisplaySettings displaySettings;
};
inline auto select(
    const std::vector< LinkerFrame >& linkers,
    const LinkerSelector& selector,
    const DisplayTypeMap& typeMap
) {
    std::vector< const LinkerFrame* > res;
    for(auto& linker : linkers) {
        if(
            (!selector.name.has_value() || *selector.name == typeMap.linkerTypeName.at(linker.type)) &&
            (!selector.subtype.has_value() || *selector.subtype == linker.subtype)
        ) {
            res.push_back(&linker);
        }
    }
    return res;
}


struct BubbleSelector {
};
struct BubbleProfile {
    // Selector.
    BubbleSelector selector;

    // Settings.
    SphereDisplaySettings displaySettings;
};
inline auto select(
    const std::vector< BubbleFrame >& bubbles,
    const BubbleSelector& selector
) {
    std::vector< const BubbleFrame* > res;
    for(auto& bubble : bubbles) {
        res.push_back(&bubble);
    }
    return res;
}

struct AuxLineProfile {
    AuxLineDisplaySettings displaySettings;
};


using ElementProfile = std::variant<
    MembraneProfile,
    FilamentProfile,
    LinkerProfile,
    BubbleProfile,
    AuxLineProfile
>;

template< typename ProfileType >
constexpr auto profileTypeDisplayName() {
    if constexpr(std::is_same_v< ProfileType, MembraneProfile >)      { return "membrane"; }
    else if constexpr(std::is_same_v< ProfileType, FilamentProfile >) { return "filament"; }
    else if constexpr(std::is_same_v< ProfileType, LinkerProfile >)   { return "linker"; }
    else if constexpr(std::is_same_v< ProfileType, BubbleProfile >)   { return "bubble"; }
    else if constexpr(std::is_same_v< ProfileType, AuxLineProfile >)  { return "auxiliary"; }
    else { return ""; }
}

// Function to get name of variant type at runtime
template< int current = 0 >
inline const char* elementProfileDisplayName(int index) {
    if constexpr(current >= std::variant_size_v< ElementProfile >) {
        throw std::runtime_error("Element profile index " + std::to_string(current + index) + " out of bounds");
    }
    else {
        return index == 0 ?
            profileTypeDisplayName< std::variant_alternative_t< current, ElementProfile > >() :
            elementProfileDisplayName< current + 1 >(index - 1);
    }
}

// Function to make profile with certain variant index at runtime.
//
// Adapted from https://stackoverflow.com/a/60567091/7120360
template< int current = 0 >
inline ElementProfile makeProfileWithIndex(int index) {
    if constexpr(current >= std::variant_size_v< ElementProfile >) {
        throw std::runtime_error("Element profile index " + std::to_string(current + index) + " out of bounds");
    }
    else {
        return index == 0 ?
            ElementProfile { std::in_place_index< current > } :
            makeProfileWithIndex< current + 1 >(index - 1);
    }
}

struct MeshDataFrameInfo {
    int start = 0;
    int size = 0;
};
struct ProfileWithMeshData {
    ElementProfile                                   profile;
    MeshData                                         data;
    std::optional< std::vector< MeshDataFrameInfo >> frameInfo;

    // Whether to inform the visualization to update mesh state.
    // If true, mesh data is outdated, and should be updated.
    //
    // It should be set to true by whoever changes the profile or source data.
    // Examples:
    //   - Profile modification
    //   - New frame in playback
    //   - Updated source data (file load / realtime system data update)
    //
    // It should be set to false by the mesh data creation/update function.
    bool shouldUpdateMeshData = true;
};


//-------------------------------------
// Functions for creating mesh data
//-------------------------------------



inline auto makeDefaultElementProfileData() {
    std::vector< ProfileWithMeshData > profiles;

    profiles.push_back(ProfileWithMeshData { MembraneProfile {} });

    profiles.push_back(ProfileWithMeshData { FilamentProfile {} });

    {
        LinkerProfile lp;
        lp.selector.name = "linker";
        profiles.push_back(ProfileWithMeshData { std::move(lp) } );
    }
    {
        LinkerProfile lp;
        lp.selector.name = "motor";
        lp.displaySettings.colorFixed = Vec3f { 0.1f, 0.1f, 0.99f };
        profiles.push_back(ProfileWithMeshData { std::move(lp) } );
    }
    {
        LinkerProfile lp;
        lp.selector.name = "brancher";
        lp.displaySettings.colorFixed = Vec3f { 0.95f, 0.8f, 0.1f };
        lp.displaySettings.pathExtrudeRadius = 12.0f;
        profiles.push_back(ProfileWithMeshData { std::move(lp) } );
    }

    profiles.push_back(ProfileWithMeshData { BubbleProfile {} });

    profiles.push_back(ProfileWithMeshData { AuxLineProfile {
        { AuxLineDisplaySettings::targetCompartmentBorder }
    }});

    return profiles;
}


inline auto displayGeometryType(const ElementProfile& profile) {
    return std::visit([](const auto& actualProfile) { return displayGeometryType(actualProfile.displaySettings); }, profile);
}

inline auto createMeshData(
    const DisplayFrame&   frameData,
    const DisplayTypeMap& typeMap,
    ElementProfile&       profile
) {
    return std::visit(
        Overload {
            [&](MembraneProfile& membraneProfile) {
                return createMembraneMeshData(
                    membraneProfile.displayStates,
                    select(frameData.membranes, membraneProfile.selector),
                    membraneProfile.displaySettings
                );
            },
            [&](const FilamentProfile& filamentProfile) {
                return createFilamentMeshData(
                    select(frameData.filaments, filamentProfile.selector),
                    filamentProfile.displaySettings
                );
            },
            [&](const LinkerProfile& linkerProfile) {
                return createLinkerMeshData(
                    select(frameData.linkers, linkerProfile.selector, typeMap),
                    linkerProfile.displaySettings
                );
            },
            [&](const BubbleProfile& bubbleProfile) {
                return createBubbleMeshData(
                    select(frameData.bubbles, bubbleProfile.selector),
                    bubbleProfile.displaySettings
                );
            },
            [&](const AuxLineProfile& auxLineProfile) {
                return createAuxLineMeshData(
                    frameData,
                    auxLineProfile.displaySettings
                );
            },
        },
        profile
    );
}


// Returns
//   - The new mesh data created, containing all the frames
//   - The frame infomation
inline auto createMeshDataAllFrames(
    ElementProfile&       profile,
    const DisplayData&    displayData
) {
    MeshData meshData;
    std::vector< MeshDataFrameInfo > frameInfo;

    meshData.descriptor = displayGeometryType(profile) == DisplayGeometryType::surface ? 
        meshDataDescriptorSurface :
        meshDataDescriptorLine;

    const int numFrames = displayData.frames.size();
    frameInfo.reserve(numFrames);
    for(int i = 0; i < numFrames; ++i) {
        const int sizeCur = meshData.data.size();
        const int sizeInc = std::visit(
            Overload {
                [&](MembraneProfile& prof) {
                    return appendMembraneMeshData(
                        meshData,
                        prof.displayStates,
                        select(displayData.frames[i].membranes, prof.selector),
                        prof.displaySettings
                    );
                },
                [&](const FilamentProfile& prof) {
                    return appendFilamentMeshData(
                        meshData,
                        select(displayData.frames[i].filaments, prof.selector),
                        prof.displaySettings
                    );
                },
                [&](const LinkerProfile& prof) {
                    return appendLinkerMeshData(
                        meshData,
                        select(displayData.frames[i].linkers, prof.selector, displayData.displayTypeMap),
                        prof.displaySettings
                    );
                },
                [&](const BubbleProfile& prof) {
                    return appendBubbleMeshData(
                        meshData,
                        select(displayData.frames[i].bubbles, prof.selector),
                        prof.displaySettings
                    );
                },
                [&](const AuxLineProfile& prof) {
                    return appendAuxLineMeshData(
                        meshData,
                        displayData.frames[i],
                        prof.displaySettings
                    );
                },
            },
            profile
        );

        frameInfo.push_back({ sizeCur, sizeInc });
    }

    return std::tuple { meshData, frameInfo };
}


//-------------------------------------
// Functions for drawing
//-------------------------------------

// Auxiliary function to replace buffer data in OpenGL
template< typename T >
inline void replaceBuffer(GLenum target, const std::vector<T>& source, int start, int size) {
    GLint prevSize;
    glGetBufferParameteriv(target, GL_BUFFER_SIZE, &prevSize);

    const auto newSize = sizeof(T) * size;

    if(newSize > prevSize) {
        glBufferData(target, newSize, source.data() + start, GL_DYNAMIC_DRAW);
    } else {
        glBufferSubData(target, 0, newSize, source.data() + start);
    }
}
template< typename T >
inline void replaceBuffer(GLenum target, const std::vector<T>& source) {
    replaceBuffer(target, source, 0, source.size());
}


// Draw elements using the mesh data.
//
// Notes:
//   - This function must be used in an OpenGL context
inline void draw(
    MeshData& meshData,
    std::optional< MeshDataFrameInfo > frameInfo,
    const GlVertexBufferManager& vertexBufferManager,
    GLenum polygonMode,
    GLenum elementMode
) {
    glBindVertexArray(vertexBufferManager.vao());

    glPolygonMode(GL_FRONT_AND_BACK, polygonMode);

    // Update data
    if(meshData.updated || frameInfo.has_value()) {
        glBindBuffer(GL_ARRAY_BUFFER, vertexBufferManager.vbo());
        if(frameInfo.has_value()) {
            replaceBuffer(GL_ARRAY_BUFFER, meshData.data, frameInfo->start, frameInfo->size);
        } else {
            replaceBuffer(GL_ARRAY_BUFFER, meshData.data);
        }
        meshData.updated = false;
    }

    // Draw
    const int numStrides = (frameInfo.has_value() ? frameInfo->size : meshData.data.size()) / meshData.descriptor.strideSize;
    glDrawArrays(elementMode, 0, numStrides);
    // glDrawElements(ve->state.eleMode, ve->state.vertexIndices.size(), GL_UNSIGNED_INT, (void*)0);

}
template<
    typename GeometryDisplaySettings,
    std::enable_if_t<
        std::is_same_v< GeometryDisplaySettings, SurfaceDisplaySettings > ||
        std::is_same_v< GeometryDisplaySettings, LineDisplaySettings >
    >* = nullptr
>
inline void draw(
    MeshData&                          meshData,
    std::optional< MeshDataFrameInfo > frameInfo,
    const Shader&                      shader,
    const GeometryDisplaySettings&     settings
) {
    if(settings.enabled) {
        // setting uniform
        if constexpr(std::is_same_v< GeometryDisplaySettings, SurfaceDisplaySettings >) {
            shader.setVec3("material.specular", convertToGlm(settings.colorSpecular));
            shader.setFloat("material.shininess", settings.colorShininess);
        }

        draw(meshData, frameInfo, settings.vertexBufferManager, polygonModeGl(settings), elementModeGl(settings));
    }
}

inline void draw(MeshData& meshData, std::optional< MeshDataFrameInfo > frameInfo, const Shader& shader, const MembraneDisplaySettings& settings) {
    draw(meshData, frameInfo, shader, settings.surface);
}
inline void draw(MeshData& meshData, std::optional< MeshDataFrameInfo > frameInfo, const Shader& shader, const FilamentDisplaySettings& settings) {
    if(settings.pathMode == FilamentDisplaySettings::PathMode::line) {
        draw(meshData, frameInfo, shader, settings.line);
    } else {
        draw(meshData, frameInfo, shader, settings.surface);
    }
}
inline void draw(MeshData& meshData, std::optional< MeshDataFrameInfo > frameInfo, const Shader& shader, const LinkerDisplaySettings& settings) {
    if(settings.pathMode == LinkerDisplaySettings::PathMode::line) {
        draw(meshData, frameInfo, shader, settings.line);
    } else {
        draw(meshData, frameInfo, shader, settings.surface);
    }
}
inline void draw(MeshData& meshData, std::optional< MeshDataFrameInfo > frameInfo, const Shader& shader, const SphereDisplaySettings& settings) {
    draw(meshData, frameInfo, shader, settings.surface);
}
inline void draw(MeshData& meshData, std::optional< MeshDataFrameInfo > frameInfo, const Shader& shader, const AuxLineDisplaySettings& settings) {
    draw(meshData, frameInfo, shader, settings.line);
}



} // namespace medyan::visual

#endif

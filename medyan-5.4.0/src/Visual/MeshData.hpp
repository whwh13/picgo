#ifndef MEDYAN_Visual_MeshData_hpp
#define MEDYAN_Visual_MeshData_hpp

// This file defines the data structure for mesh data stored in the memory,
// which can then be mapped into GPU memory for drawing.

#include <array>
#include <iterator>
#include <numeric>
#include <tuple>

#include <Eigen/Geometry>

#include "Util/Math/ColorMap.hpp"
#include "Visual/Common.hpp"
#include "Visual/FrameData.hpp"
#include "Visual/Geometry/PathExtrude.hpp"
#include "Visual/Geometry/Sphere.hpp"
#include "Visual/Shader.hpp"

namespace medyan::visual {

//--------------------------------------
// Data structures
//--------------------------------------

struct MeshDataDescriptor {
    int strideSize = 9;

    int positionStart = 0;
    int positionSize = 3;
    int normalStart = 3;
    int normalSize = 3;
    int colorStart = 6;
    int colorSize = 3;
};
// preset descriptors
inline constexpr MeshDataDescriptor meshDataDescriptorSurface { 9, 0, 3, 3, 3, 6, 3 };
inline constexpr MeshDataDescriptor meshDataDescriptorLine    { 6, 0, 3, 3, 0, 3, 3 };

class GlVertexBufferManager {
private:
    // vao, vbo, ebo
    GLuint vao_ = 0;
    GLuint vbo_ = 0;
    // GLuint ebo_;

    bool bufferMade_ = false;

public:

    auto vao() const { return vao_; }
    auto vbo() const { return vbo_; }

    auto bufferMade() const { return bufferMade_; }

    // ctor and dtor
    GlVertexBufferManager() = default;
    GlVertexBufferManager(MeshDataDescriptor desc) {
        makeBuffer(desc);
    }
    GlVertexBufferManager(GlVertexBufferManager&& rhs) {
        swap(*this, rhs);
    }

    ~GlVertexBufferManager() {
        deleteBuffer();
    }

    GlVertexBufferManager& operator=(GlVertexBufferManager&& rhs) {
        if(this != &rhs) {
            // clear this
            deleteBuffer();

            swap(*this, rhs);
        }

        return *this;
    }

    friend void swap(GlVertexBufferManager& lhs, GlVertexBufferManager& rhs) noexcept {
        std::swap(lhs.vao_,        rhs.vao_);
        std::swap(lhs.vbo_,        rhs.vbo_);
        std::swap(lhs.bufferMade_, rhs.bufferMade_);
    }

    void makeBuffer(MeshDataDescriptor desc) {
        if(bufferMade_) deleteBufferImpl_();

        makeBufferImpl_(desc);
        bufferMade_ = true;
    }
    void deleteBuffer() {
        if(bufferMade_) {
            deleteBufferImpl_();
        }
        bufferMade_ = false;
    }

private:
    // Note: this function does not delete previously built buffer
    void makeBufferImpl_(MeshDataDescriptor desc) {
        glGenBuffers(1, &vbo_);
        // glGenBuffers(1, &ebo);
        glGenVertexArrays(1, &vao_);

        glBindVertexArray(vao_);
        glBindBuffer(GL_ARRAY_BUFFER, vbo_);
        // glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo_);

        // Vertex attribute
        //---------------------------------------------------------------------
        GLuint idx = 0;
        const auto addBlock = [&](int blockOffset, int blockSize, int strideSize) {
            if(blockSize > 0) {
                glVertexAttribPointer(
                    idx,
                    blockSize,
                    GL_FLOAT,
                    GL_FALSE,
                    strideSize * sizeof(float),
                    static_cast<const char*>(0) + sizeof(float) * blockOffset
                );
                glEnableVertexAttribArray(idx);
                ++idx;
            }
        };
        // Position
        addBlock(desc.positionStart, desc.positionSize, desc.strideSize);
        addBlock(desc.normalStart,   desc.normalSize,   desc.strideSize);
        addBlock(desc.colorStart,    desc.colorSize,    desc.strideSize);

        // temporarily retarget
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        // glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        glBindVertexArray(0);
    }

    // Note: this function does not check for creation
    void deleteBufferImpl_() {
        glDeleteVertexArrays(1, &vao_);
        glDeleteBuffers(1, &vbo_);
        // glDeleteBuffers(1, &ebo_);
    }
};


struct MeshData {
    MeshDataDescriptor descriptor;
    std::vector< float > data;

    // Auxiliary state to indicate that the data is updated.
    // Can be set to false by the downstream pipeline.
    bool updated = true;
};




//--------------------------------------
// Settings (for display)
//--------------------------------------

enum class DisplayGeometryType { surface, line };

struct SurfaceDisplaySettings {
    enum class PolygonMode { wireframe, fill, last_ };

    bool            enabled = true;
    PolygonMode     polygonMode = PolygonMode::fill;

    Vec3f           colorSpecular { 0.5, 0.5, 0.5 };
    float           colorShininess = 32.0f;

    GlVertexBufferManager vertexBufferManager;

    SurfaceDisplaySettings() : vertexBufferManager(meshDataDescriptorSurface) {}
};

constexpr auto text(SurfaceDisplaySettings::PolygonMode value) {
    switch(value) {
        case SurfaceDisplaySettings::PolygonMode::wireframe: return "wireframe";
        case SurfaceDisplaySettings::PolygonMode::fill:      return "fill";
        default:                                             return "";
    }
}

struct LineDisplaySettings {
    bool enabled = true;

    GlVertexBufferManager vertexBufferManager;

    LineDisplaySettings() : vertexBufferManager(meshDataDescriptorLine) {}
};


//--------------------------------------
// Settings for making mesh data.
// States updated during making mesh data.
//--------------------------------------

struct MembraneDisplaySettings {
    enum class ElementMode { triangle, edge, last_ };
    enum class ColorMode { fixed, attribute, last_ };

    ElementMode elementMode = ElementMode::triangle;

    float edgeExtrudeRadius = 5.0f; // used only when element mode is "edge".
    int   edgeExtrudeSides  = 8;

    ColorMode                        colorMode = ColorMode::fixed;
    Vec3f                            colorFixed { 0.4f, 0.6f, 0.95f };
    Index                            attributeIndex = -1;
    bool                             autoAttributeRange = true;
    float                            manualAttribMin = -1;
    float                            manualAttribMax = 1;
    colormap::DynamicColorMap<float> colorMap = colormap::bwrf;

    SurfaceDisplaySettings surface;
};
constexpr auto text(MembraneDisplaySettings::ElementMode value) {
    switch(value) {
        case MembraneDisplaySettings::ElementMode::triangle: return "triangle";
        case MembraneDisplaySettings::ElementMode::edge:     return "edge";
        default:                                             return "";
    }
}
constexpr auto text(MembraneDisplaySettings::ColorMode value) {
    switch(value) {
        case MembraneDisplaySettings::ColorMode::fixed:     return "fixed";
        case MembraneDisplaySettings::ColorMode::attribute: return "attribute";
        default:                                            return "";
    }
}

struct MembraneDisplayStates {
    struct AttribRange {
        float min = inff;
        float max = -inff;
    };
    // Range of attribute values.
    // If all frames are generated in batch, it is indexed by frame index.
    // Otherwise, it contains one entry for the current frame.
    std::vector<AttribRange> attribRanges;
};


struct FilamentDisplaySettings {
    enum class PathMode { line, extrude, bead, last_ };

    Vec3f colorFixed { 0.95f, 0.1f, 0.15f };

    PathMode pathMode = PathMode::extrude;

    // path extrude parameters
    float  pathExtrudeRadius = 7.5f;
    int    pathExtrudeSides = 10;

    // bead parameters
    float  beadRadius = 12.0f;
    int    beadLongitudeSegs = 10;
    int    beadLatitudeSegs = 5;

    // display settings
    SurfaceDisplaySettings surface;  // used in "extrude" or "bead" path mode
    LineDisplaySettings    line;     // used in "line" path mode

};

constexpr auto text(FilamentDisplaySettings::PathMode value) {
    switch(value) {
        case FilamentDisplaySettings::PathMode::line:    return "line";
        case FilamentDisplaySettings::PathMode::extrude: return "extrude";
        case FilamentDisplaySettings::PathMode::bead:    return "bead";
        default:                                         return "";
    }
}

struct LinkerDisplaySettings {
    enum class PathMode { line, extrude, last_ };

    Vec3f colorFixed { 0.1f, 0.9f, 0.0f };

    PathMode pathMode = PathMode::extrude;

    // path extrude parameters
    float  pathExtrudeRadius = 7.5f;
    int    pathExtrudeSides = 10;

    // display settings
    SurfaceDisplaySettings surface;  // used in "extrude" or "bead" path mode
    LineDisplaySettings    line;     // used in "line" path mode

};

constexpr auto text(LinkerDisplaySettings::PathMode value) {
    switch(value) {
        case LinkerDisplaySettings::PathMode::line:    return "line";
        case LinkerDisplaySettings::PathMode::extrude: return "extrude";
        default:                                       return "";
    }
}

struct SphereDisplaySettings {
    Vec3f colorFixed { 0.95f, 0.8f, 0.1f };

    // Sphere display parameters.
    bool  useManualRadius = false;
    float manualRadius = 12.0f;
    int   sphereLongitudeSegs = 13;
    int   sphereLatitudeSegs = 7;

    // display settings
    SurfaceDisplaySettings surface;
};

struct AuxLineDisplaySettings {
    using Flag = std::uint_fast8_t;
    inline static constexpr Flag targetCompartmentBorder = 1 << 0;
    inline static constexpr Flag targetCompartmentAll    = 1 << 1;

    Flag flag = 0;

    Vec3f colorFixed { 1.0f, 1.0f, 1.0f };

    LineDisplaySettings line;
};



//--------------------------------------
// Functions (mesh data creation)
//--------------------------------------

// Note that the following functions may be categorized as the following:
//   - The actual function appending data of each element to the mesh data.
//   - The function appending data of all selected elements of a frame to existing mesh data.
//   - The function creating a new mesh data solely for the data of all selected elements of a frame.

// For membrane mesh data creation, find the maximum and minimum values of the attribute.
// Returns { attribMin, attribMax }.
inline auto findMinMaxMembraneAttribute(
    const std::vector<const MembraneFrame*>& pMembranes,
    const MembraneDisplaySettings&           membraneSettings
) {
    float attribMin = inff;
    float attribMax = -inff;

    const bool useAttrib = membraneSettings.colorMode == MembraneDisplaySettings::ColorMode::attribute && membraneSettings.attributeIndex >= 0;
    if(useAttrib) {
        for(auto pm : pMembranes) {
            const auto numVertices = pm->vertexAttributes.cols();
            for(Index vi = 0; vi < numVertices; ++vi) {
                const auto attrib = pm->vertexAttributes(membraneSettings.attributeIndex, vi);
                attribMin = std::min(attribMin, attrib);
                attribMax = std::max(attribMax, attrib);
            }
        }
    }

    return std::tuple { attribMin, attribMax };
}

// Returns how many numbers added
inline auto appendMembraneMeshData(
    MeshData&                      meshData,
    float                          actualAttribMin,
    float                          actualAttribMax,
    const MembraneFrame&           membrane,
    const MembraneDisplaySettings& membraneSettings
) {
    using namespace std;

    const Size sizePrev = meshData.data.size();

    // Preprosessing.
    const bool useAttrib = membraneSettings.colorMode == MembraneDisplaySettings::ColorMode::attribute && membraneSettings.attributeIndex >= 0;
    // Find range of the attribute.
    float attribMin = membraneSettings.autoAttributeRange ? actualAttribMin : membraneSettings.manualAttribMin;
    float attribMax = membraneSettings.autoAttributeRange ? actualAttribMax : membraneSettings.manualAttribMax;

    const auto attrib01 = [&](float attrib) -> float {
        if(attribMin < attribMax) {
            return (attrib - attribMin) / (attribMax - attribMin);
        }
        return 0.5f;
    };
    const auto getColor = [&](float attrib) {
        return useAttrib
            ? colormap::color(membraneSettings.colorMap, attrib01(attrib))
            : membraneSettings.colorFixed;
    };

    switch(membraneSettings.elementMode) {
        case MembraneDisplaySettings::ElementMode::triangle:

            for(const auto& t : membrane.triangles) {
                Eigen::Map<const Eigen::Vector3f> c[] {
                    Eigen::Vector3f::Map(&membrane.vertexAttributes(0, t[0])),
                    Eigen::Vector3f::Map(&membrane.vertexAttributes(0, t[1])),
                    Eigen::Vector3f::Map(&membrane.vertexAttributes(0, t[2])),
                };
                const Eigen::Vector3f un = (c[1] - c[0]).cross(c[2] - c[0]).normalized();

                for(int i = 0; i < 3; ++i) {
                    const auto color = getColor(useAttrib ? membrane.vertexAttributes(membraneSettings.attributeIndex, t[i]) : 0);
                    meshData.data.push_back(c[i][0]);
                    meshData.data.push_back(c[i][1]);
                    meshData.data.push_back(c[i][2]);
                    meshData.data.push_back(un[0]);
                    meshData.data.push_back(un[1]);
                    meshData.data.push_back(un[2]);
                    meshData.data.push_back(color[0]);
                    meshData.data.push_back(color[1]);
                    meshData.data.push_back(color[2]);
                }
            }
            break;

        case MembraneDisplaySettings::ElementMode::edge:

            // For each triangle, find all three edges and make a path extrusion for it.
            //
            // Currently, each non-border edge will be generated twice, which is redundant.
            for(const auto& memTriangle : membrane.triangles) {
                for(int edgeIdx = 0; edgeIdx < 3; ++edgeIdx) {
                    array< int, 2 > indices { memTriangle[edgeIdx], memTriangle[(edgeIdx + 1) % 3] };

                    const auto [genVertices, genVertexNormals, genAttribs, genTriInd]
                        = pathExtrudeGenerateWithAttrib<float>(
                            membrane.vertexAttributes,
                            // Function that gets certain attribute.
                            [&](const auto& attributes, Index index) -> float {
                                return useAttrib
                                    ? attributes(membraneSettings.attributeIndex, index)
                                    : 0;
                            },
                            indices,
                            membraneSettings.edgeExtrudeRadius,
                            membraneSettings.edgeExtrudeSides
                        );

                    const int numDisplayTriangles = genTriInd.size();
                    for(int t = 0; t < numDisplayTriangles; ++t) {
                        const auto& triInds = genTriInd[t];

                        for(int i = 0; i < 3; ++i) {
                            const auto coord = genVertices[triInds[i]];
                            const auto un    = genVertexNormals[triInds[i]];
                            const auto color = getColor(genAttribs[triInds[i]]);
                            meshData.data.push_back(coord[0]);
                            meshData.data.push_back(coord[1]);
                            meshData.data.push_back(coord[2]);
                            meshData.data.push_back(un[0]);
                            meshData.data.push_back(un[1]);
                            meshData.data.push_back(un[2]);
                            meshData.data.push_back(color[0]);
                            meshData.data.push_back(color[1]);
                            meshData.data.push_back(color[2]);
                        }
                    }
                }
            }
            break;

    }

    const Size sizeCur = meshData.data.size();
    return sizeCur - sizePrev;
}

// Note:
//   - This function does not set descriptor.
inline auto appendMembraneMeshData(
    MeshData&                         meshData,
    MembraneDisplayStates&            membraneStates,
    std::vector<const MembraneFrame*> pMembranes,
    const MembraneDisplaySettings&    membraneSettings
) {
    // Find range of the attribute.
    const auto [attribMin, attribMax] = findMinMaxMembraneAttribute(pMembranes, membraneSettings);
    membraneStates.attribRanges.push_back({ attribMin, attribMax });

    Size numAdded = 0;
    for(auto pm : pMembranes) {
        numAdded += appendMembraneMeshData(meshData, attribMin, attribMax, *pm, membraneSettings);
    }

    return numAdded;
}

// Generate membrane mesh with a selected range of membranes.
//
// The selected membranes should be replaced with a range object with C++20.
inline auto createMembraneMeshData(
    MembraneDisplayStates&            membraneStates,
    std::vector<const MembraneFrame*> pMembranes,
    const MembraneDisplaySettings&    membraneSettings
) {
    MeshData res;
    res.descriptor = meshDataDescriptorSurface;

    // Find range of the attribute.
    const auto [attribMin, attribMax] = findMinMaxMembraneAttribute(pMembranes, membraneSettings);
    membraneStates.attribRanges.resize(1);
    membraneStates.attribRanges[0] = { attribMin, attribMax };

    int numTriangles = 0;
    for(auto pm : pMembranes) numTriangles += pm->triangles.size();
    res.data.reserve(3 * numTriangles * res.descriptor.strideSize);
    for(auto pm : pMembranes) {
        appendMembraneMeshData(res, attribMin, attribMax, *pm, membraneSettings);
    }

    return res;
}

// Make mesh data for a single filament and append it to the mesh data.
//
// Returns how many numbers added.
inline auto appendFilamentMeshData(
    MeshData&                      meshData,
    const FilamentFrame&           filament,
    const FilamentDisplaySettings& filamentSettings
) {
    using PM = FilamentDisplaySettings::PathMode;
    using namespace std;

    const int sizePrev = meshData.data.size();
    const int numBeads = filament.coords.size();
    switch(filamentSettings.pathMode) {
        case PM::line:
            if(numBeads > 1) {
                const int numSegments = numBeads - 1;

                for(int i = 0; i < numSegments; ++i) {
                    meshData.data.push_back(filament.coords[i][0]);
                    meshData.data.push_back(filament.coords[i][1]);
                    meshData.data.push_back(filament.coords[i][2]);
                    meshData.data.push_back(filamentSettings.colorFixed[0]);
                    meshData.data.push_back(filamentSettings.colorFixed[1]);
                    meshData.data.push_back(filamentSettings.colorFixed[2]);

                    meshData.data.push_back(filament.coords[i+1][0]);
                    meshData.data.push_back(filament.coords[i+1][1]);
                    meshData.data.push_back(filament.coords[i+1][2]);
                    meshData.data.push_back(filamentSettings.colorFixed[0]);
                    meshData.data.push_back(filamentSettings.colorFixed[1]);
                    meshData.data.push_back(filamentSettings.colorFixed[2]);
                }
            }
            break;

        case PM::extrude:

            {
                vector< int > trivialIndices(numBeads);
                iota(begin(trivialIndices), end(trivialIndices), 0);
                const auto [genVertices, genVertexNormals, genTriInd]
                    = pathExtrudeGenerate<float>(
                        filament.coords,
                        trivialIndices,
                        filamentSettings.pathExtrudeRadius,
                        filamentSettings.pathExtrudeSides
                    );

                const int numTriangles = genTriInd.size();
                for(int t = 0; t < numTriangles; ++t) {
                    const auto& triInds = genTriInd[t];
                    const Vec3f coord[] {
                        genVertices[triInds[0]],
                        genVertices[triInds[1]],
                        genVertices[triInds[2]]
                    };
                    const Vec3f un[] {
                        genVertexNormals[triInds[0]],
                        genVertexNormals[triInds[1]],
                        genVertexNormals[triInds[2]]
                    };

                    for(int i = 0; i < 3; ++i) {
                        meshData.data.push_back(coord[i][0]);
                        meshData.data.push_back(coord[i][1]);
                        meshData.data.push_back(coord[i][2]);
                        meshData.data.push_back(un[i][0]);
                        meshData.data.push_back(un[i][1]);
                        meshData.data.push_back(un[i][2]);
                        meshData.data.push_back(filamentSettings.colorFixed[0]);
                        meshData.data.push_back(filamentSettings.colorFixed[1]);
                        meshData.data.push_back(filamentSettings.colorFixed[2]);
                    }
                }
            }
            break;

        case PM::bead:

            {
                const auto sphereGen = SphereUv<float> {
                    filamentSettings.beadRadius,
                    filamentSettings.beadLongitudeSegs,
                    filamentSettings.beadLatitudeSegs
                };
                const auto sphereCache = sphereGen.makeCache();

                for(const auto& bc : filament.coords) {

                    const auto [genVertices, _] = sphereGen.generate(
                        {
                            static_cast<float>(bc[0]),
                            static_cast<float>(bc[1]),
                            static_cast<float>(bc[2])
                        },
                        sphereCache
                    );

                    const int numTriangles = sphereCache.triInd.size();
                    for(int t = 0; t < numTriangles; ++t) {
                        const typename decltype(genVertices)::value_type coord[] {
                            genVertices[sphereCache.triInd[t][0]],
                            genVertices[sphereCache.triInd[t][1]],
                            genVertices[sphereCache.triInd[t][2]]
                        };
                        const auto un = normalizedVector(cross(coord[1] - coord[0], coord[2] - coord[0]));

                        for(int i = 0; i < 3; ++i) {
                            meshData.data.push_back(coord[i][0]);
                            meshData.data.push_back(coord[i][1]);
                            meshData.data.push_back(coord[i][2]);
                            meshData.data.push_back(un[0]);
                            meshData.data.push_back(un[1]);
                            meshData.data.push_back(un[2]);
                            meshData.data.push_back(filamentSettings.colorFixed[0]);
                            meshData.data.push_back(filamentSettings.colorFixed[1]);
                            meshData.data.push_back(filamentSettings.colorFixed[2]);
                        }
                    }
                } // end loop beads in a filament

            }
            break;
    } // end switch

    const Size sizeCur = meshData.data.size();
    return sizeCur - sizePrev;
}

// Note:
//   - This function does not set descriptor.
inline auto appendFilamentMeshData(
    MeshData&                         meshData,
    std::vector<const FilamentFrame*> pFilaments,
    const FilamentDisplaySettings&    filamentSettings
) {
    Size numAdded = 0;
    for(auto pf : pFilaments) {
        numAdded += appendFilamentMeshData(meshData, *pf, filamentSettings);
    }

    return numAdded;
}

// Generate filament mesh with selected filaments
//
// The selected objects should be replaced with a range object with C++20.
inline auto createFilamentMeshData(
    std::vector<const FilamentFrame*> pFilaments,
    const FilamentDisplaySettings&    filamentSettings
) {
    using PM = FilamentDisplaySettings::PathMode;

    MeshData res;

    // reserve space and set descriptor
    switch(filamentSettings.pathMode) {
        case PM::line:
            res.descriptor = meshDataDescriptorLine;
            {
                int numSegments = 0;
                for(auto pf : pFilaments) {
                    if(pf->coords.size() >= 2) {
                        numSegments += pf->coords.size() - 1;
                    }
                }
                res.data.reserve(2 * numSegments * res.descriptor.strideSize);
            }
            break;

        case PM::extrude:
            res.descriptor = meshDataDescriptorSurface;
            {
                int numTriangles = 0;
                for(auto pf : pFilaments) {
                    numTriangles += pathExtrudeEstimateNumTriangles(
                        pf->coords.size(),
                        filamentSettings.pathExtrudeSides
                    );
                }
                res.data.reserve(3 * numTriangles * res.descriptor.strideSize);
            }
            break;

        case PM::bead:
            res.descriptor = meshDataDescriptorSurface;
            {
                const int numTrianglesPerBead = SphereUv<float> {
                    filamentSettings.beadRadius,
                    filamentSettings.beadLongitudeSegs,
                    filamentSettings.beadLatitudeSegs
                }.estimateNumTriangles();

                int numBeads = 0;
                for(auto pf : pFilaments) {
                    numBeads += pf->coords.size();
                }

                res.data.reserve(3 * numTrianglesPerBead * numBeads * res.descriptor.strideSize);
            }
            break;
    }

    for(auto pf : pFilaments) {
        appendFilamentMeshData(res, *pf, filamentSettings);
    }

    return res;
}


// Make mesh data for a single filament and append it to the mesh data.
//
// Returns how many numbers added.
inline auto appendLinkerMeshData(
    MeshData&                    meshData,
    const LinkerFrame&           linker,
    const LinkerDisplaySettings& linkerSettings
) {
    using PM = LinkerDisplaySettings::PathMode;
    using namespace std;

    const int sizePrev = meshData.data.size();
    switch(linkerSettings.pathMode) {
        case PM::line:
            {
                meshData.data.push_back(linker.coords[0][0]);
                meshData.data.push_back(linker.coords[0][1]);
                meshData.data.push_back(linker.coords[0][2]);
                meshData.data.push_back(linkerSettings.colorFixed[0]);
                meshData.data.push_back(linkerSettings.colorFixed[1]);
                meshData.data.push_back(linkerSettings.colorFixed[2]);

                meshData.data.push_back(linker.coords[1][0]);
                meshData.data.push_back(linker.coords[1][1]);
                meshData.data.push_back(linker.coords[1][2]);
                meshData.data.push_back(linkerSettings.colorFixed[0]);
                meshData.data.push_back(linkerSettings.colorFixed[1]);
                meshData.data.push_back(linkerSettings.colorFixed[2]);
            }
            break;

        case PM::extrude:

            {
                vector< int > trivialIndices { 0, 1 };
                const auto [genVertices, genVertexNormals, genTriInd]
                    = pathExtrudeGenerate<float>(
                        linker.coords,
                        trivialIndices,
                        linkerSettings.pathExtrudeRadius,
                        linkerSettings.pathExtrudeSides
                    );

                const int numTriangles = genTriInd.size();
                for(int t = 0; t < numTriangles; ++t) {
                    const auto& triInds = genTriInd[t];
                    const Vec3f coord[] {
                        genVertices[triInds[0]],
                        genVertices[triInds[1]],
                        genVertices[triInds[2]]
                    };
                    const Vec3f un[] {
                        genVertexNormals[triInds[0]],
                        genVertexNormals[triInds[1]],
                        genVertexNormals[triInds[2]]
                    };

                    for(int i = 0; i < 3; ++i) {
                        meshData.data.push_back(coord[i][0]);
                        meshData.data.push_back(coord[i][1]);
                        meshData.data.push_back(coord[i][2]);
                        meshData.data.push_back(un[i][0]);
                        meshData.data.push_back(un[i][1]);
                        meshData.data.push_back(un[i][2]);
                        meshData.data.push_back(linkerSettings.colorFixed[0]);
                        meshData.data.push_back(linkerSettings.colorFixed[1]);
                        meshData.data.push_back(linkerSettings.colorFixed[2]);
                    }
                }
            }
            break;

    } // end switch

    const Size sizeCur = meshData.data.size();
    return sizeCur - sizePrev;
}


// Note:
//   - This function does not set descriptor.
inline auto appendLinkerMeshData(
    MeshData&                       meshData,
    std::vector<const LinkerFrame*> pLinkers,
    const LinkerDisplaySettings&    linkerSettings
) {
    Size numAdded = 0;
    for(auto pl : pLinkers) {
        numAdded += appendLinkerMeshData(meshData, *pl, linkerSettings);
    }

    return numAdded;
}


// Create linkers given the range of linkers
//
// The selected objects should be replaced with a range object with C++20.
inline auto createLinkerMeshData(
    std::vector<const LinkerFrame*> pLinkers,
    const LinkerDisplaySettings&    linkerSettings
) {
    using PM = LinkerDisplaySettings::PathMode;

    MeshData res;

    // reserve space and set descriptor
    switch(linkerSettings.pathMode) {
        case PM::line:
            res.descriptor = meshDataDescriptorLine;
            {
                const int numSegments = pLinkers.size();
                res.data.reserve(2 * numSegments * res.descriptor.strideSize);
            }
            break;

        case PM::extrude:
            res.descriptor = meshDataDescriptorSurface;
            {
                const int numTriangles =
                    pLinkers.size() *
                    pathExtrudeEstimateNumTriangles(2, linkerSettings.pathExtrudeSides);
                res.data.reserve(3 * numTriangles * res.descriptor.strideSize);
            }
            break;
    }

    for(auto pl : pLinkers) {
        appendLinkerMeshData(res, *pl, linkerSettings);
    }

    return res;
}


// Make mesh data for a single bubble and append it to the mesh data.
//
// Returns how many numbers added.
inline auto appendBubbleMeshData(
    MeshData&                      meshData,
    const BubbleFrame&             bubble,
    const SphereDisplaySettings&   sphereSettings
) {
    using PM = FilamentDisplaySettings::PathMode;
    using namespace std;

    const Size sizePrev = meshData.data.size();

    {
        const auto sphereGen = SphereUv<float> {
            sphereSettings.useManualRadius ? sphereSettings.manualRadius : bubble.radius,
            sphereSettings.sphereLongitudeSegs,
            sphereSettings.sphereLatitudeSegs
        };
        const auto sphereCache = sphereGen.makeCache();

        const auto [genVertices, _] = sphereGen.generate(
            {
                static_cast<float>(bubble.coord[0]),
                static_cast<float>(bubble.coord[1]),
                static_cast<float>(bubble.coord[2])
            },
            sphereCache
        );

        const int numTriangles = sphereCache.triInd.size();
        for(int t = 0; t < numTriangles; ++t) {
            const typename decltype(genVertices)::value_type coord[] {
                genVertices[sphereCache.triInd[t][0]],
                genVertices[sphereCache.triInd[t][1]],
                genVertices[sphereCache.triInd[t][2]]
            };
            const auto un = normalizedVector(cross(coord[1] - coord[0], coord[2] - coord[0]));

            for(int i = 0; i < 3; ++i) {
                meshData.data.push_back(coord[i][0]);
                meshData.data.push_back(coord[i][1]);
                meshData.data.push_back(coord[i][2]);
                meshData.data.push_back(un[0]);
                meshData.data.push_back(un[1]);
                meshData.data.push_back(un[2]);
                meshData.data.push_back(sphereSettings.colorFixed[0]);
                meshData.data.push_back(sphereSettings.colorFixed[1]);
                meshData.data.push_back(sphereSettings.colorFixed[2]);
            }
        }

    }

    const Size sizeCur = meshData.data.size();
    return sizeCur - sizePrev;
}

// Note:
//   - This function does not set descriptor.
inline auto appendBubbleMeshData(
    MeshData&                       meshData,
    std::vector<const BubbleFrame*> pBubbles,
    const SphereDisplaySettings&    sphereSettings
) {
    Size numAdded = 0;
    for(auto p : pBubbles) {
        numAdded += appendBubbleMeshData(meshData, *p, sphereSettings);
    }

    return numAdded;
}

// Generate bubble mesh with selected bubbles.
//
// The selected objects should be replaced with a range object with C++20.
inline auto createBubbleMeshData(
    std::vector<const BubbleFrame*> pBubbles,
    const SphereDisplaySettings&    sphereSettings
) {
    using PM = FilamentDisplaySettings::PathMode;

    MeshData res;
    res.descriptor = meshDataDescriptorSurface;

    const int numTrianglesPerPoint = SphereUv<float> {
        1.0f, // Dummy radius.
        sphereSettings.sphereLongitudeSegs,
        sphereSettings.sphereLatitudeSegs
    }.estimateNumTriangles();

    res.data.reserve(3 * numTrianglesPerPoint * pBubbles.size() * res.descriptor.strideSize);

    for(auto p : pBubbles) {
        appendBubbleMeshData(res, *p, sphereSettings);
    }

    return res;
}


// Returns how many numbers added.
inline auto appendAuxLineMeshData(
    MeshData&                       meshData,
    const DisplayFrame&             frameData,
    const AuxLineDisplaySettings&   auxLineSettings
) {
    const int sizePrev = meshData.data.size();

    if(frameData.compartmentInfo.has_value()) {

        // Function to generate auxiliary lines
        //
        // Inputs:
        //   - fixedAxis: the axis parallel to the line to be drawn. range 0, 1, 2.
        //   - dax1:      the number of compartments interval along the next axis after fixedAxis.
        //   - dax2:      the number of compartments interval along the 2nd next axis after fixedAxis.
        //
        // Notes:
        //   - The next axis index are found using the rotation 0 -> 1 -> 2 -> 0.
        const auto genLineArray = [&](size_t fixedAxis, size_t dax1, size_t dax2) {
            const int ax1 = (fixedAxis + 1) % 3;
            const int ax2 = (fixedAxis + 2) % 3;
            for(size_t x1 = 0; x1 <= frameData.compartmentInfo->number[ax1]; x1 += dax1) {
                for(size_t x2 = 0; x2 <= frameData.compartmentInfo->number[ax2]; x2 += dax2) {
                    const auto v1 = frameData.compartmentInfo->size[ax1] * x1;
                    const auto v2 = frameData.compartmentInfo->size[ax2] * x2;
                    Vec3f coord0, coord1;
                    coord0[fixedAxis] = 0;
                    coord0[ax1] = v1;
                    coord0[ax2] = v2;
                    coord0 += frameData.compartmentInfo->offset;
                    coord1[fixedAxis] = frameData.compartmentInfo->size[fixedAxis] * frameData.compartmentInfo->number[fixedAxis];
                    coord1[ax1] = v1;
                    coord1[ax2] = v2;
                    coord1 += frameData.compartmentInfo->offset;

                    meshData.data.push_back(coord0[0]);
                    meshData.data.push_back(coord0[1]);
                    meshData.data.push_back(coord0[2]);
                    meshData.data.push_back(auxLineSettings.colorFixed[0]);
                    meshData.data.push_back(auxLineSettings.colorFixed[1]);
                    meshData.data.push_back(auxLineSettings.colorFixed[2]);
                    meshData.data.push_back(coord1[0]);
                    meshData.data.push_back(coord1[1]);
                    meshData.data.push_back(coord1[2]);
                    meshData.data.push_back(auxLineSettings.colorFixed[0]);
                    meshData.data.push_back(auxLineSettings.colorFixed[1]);
                    meshData.data.push_back(auxLineSettings.colorFixed[2]);
                }
            }
        };

        if(auxLineSettings.flag & AuxLineDisplaySettings::targetCompartmentBorder) {
            genLineArray(0, frameData.compartmentInfo->number[1], frameData.compartmentInfo->number[2]);
            genLineArray(1, frameData.compartmentInfo->number[2], frameData.compartmentInfo->number[0]);
            genLineArray(2, frameData.compartmentInfo->number[0], frameData.compartmentInfo->number[1]);
        }
        if(auxLineSettings.flag & AuxLineDisplaySettings::targetCompartmentAll) {
            genLineArray(0, 1, 1);
            genLineArray(1, 1, 1);
            genLineArray(2, 1, 1);
        }
    }

    const Size sizeCur = meshData.data.size();
    return sizeCur - sizePrev;
}
// Auxiliary lines do not need a selector.
// The actual "selection" happens inside a mesh.

inline auto createAuxLineMeshData(
    const DisplayFrame&             frameData,
    const AuxLineDisplaySettings&   auxLineSettings
) {
    MeshData res;
    res.descriptor = meshDataDescriptorLine;
    appendAuxLineMeshData(res, frameData, auxLineSettings);

    return res;
}

//-------------------------------------
// Functions (mesh data display)
//-------------------------------------


inline auto convertToGlm(const Vec3f& vec) {
    return glm::vec3( vec[0], vec[1], vec[2] );
}


// Auxiliary functions to get polygon mode
inline auto polygonModeGl(const SurfaceDisplaySettings& settings) {
    return settings.polygonMode == SurfaceDisplaySettings::PolygonMode::wireframe ? GL_LINE : GL_FILL;
}
inline auto polygonModeGl(const LineDisplaySettings& settings) {
    return GL_LINE;
}

// Auxiliary functions to get element mode
inline auto elementModeGl(const SurfaceDisplaySettings&) {
    return GL_TRIANGLES;
}
inline auto elementModeGl(const LineDisplaySettings&) {
    return GL_LINES;
}

inline auto displayGeometryType(const MembraneDisplaySettings&) {
    return DisplayGeometryType::surface;
}
inline auto displayGeometryType(const FilamentDisplaySettings& settings) {
    return settings.pathMode == FilamentDisplaySettings::PathMode::line ?
        DisplayGeometryType::line :
        DisplayGeometryType::surface;
}
inline auto displayGeometryType(const LinkerDisplaySettings& settings) {
    return settings.pathMode == LinkerDisplaySettings::PathMode::line ?
        DisplayGeometryType::line :
        DisplayGeometryType::surface;
}
inline auto displayGeometryType(const SphereDisplaySettings& settings) {
    return DisplayGeometryType::surface;
}
inline auto displayGeometryType(const AuxLineDisplaySettings&) {
    return DisplayGeometryType::line;
}


//-------------------------------------
// Functions (mesh data auxiliary information)
//-------------------------------------

// Adjust the previous bounding box to include the new mesh data.
//
// bounds may contain an array of 6 values: { xmin, xmax, ymin, ymax, zmin, zmax }.
// bounds may also be empty, indicating that the bounding box is not yet initialized.
// If the data is empty or does not contain 3D coordinates, bounds will not be modified.
inline void extend3DExtents(std::optional<std::array<float, 6>>& bounds, const MeshData& meshData) {
    if(meshData.data.empty() || meshData.descriptor.positionSize != 3) {
        return;
    }

    if(!bounds.has_value()) {
        bounds = std::array {
            std::numeric_limits<float>::infinity(),
            -std::numeric_limits<float>::infinity(),
            std::numeric_limits<float>::infinity(),
            -std::numeric_limits<float>::infinity(),
            std::numeric_limits<float>::infinity(),
            -std::numeric_limits<float>::infinity()
        };
    }
    float& xmin = (*bounds)[0];
    float& xmax = (*bounds)[1];
    float& ymin = (*bounds)[2];
    float& ymax = (*bounds)[3];
    float& zmin = (*bounds)[4];
    float& zmax = (*bounds)[5];

    for(Index i = 0; i < meshData.data.size(); i += meshData.descriptor.strideSize) {
        const float x = meshData.data[i + meshData.descriptor.positionStart];
        const float y = meshData.data[i + meshData.descriptor.positionStart + 1];
        const float z = meshData.data[i + meshData.descriptor.positionStart + 2];
        xmin = std::min(xmin, x);
        xmax = std::max(xmax, x);
        ymin = std::min(ymin, y);
        ymax = std::max(ymax, y);
        zmin = std::min(zmin, z);
        zmax = std::max(zmax, z);
    }
}


} // namespace medyan::visual

#endif

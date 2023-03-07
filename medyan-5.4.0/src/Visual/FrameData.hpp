#ifndef MEDYAN_Visual_FrameData_hpp
#define MEDYAN_Visual_FrameData_hpp

#include <array>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <optional>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>

#include "Output.h"
#include "Util/Io/H5.hpp"
#include "Util/Math/Vec.hpp"

namespace medyan::visual {


//-------------------------------------
// Data structures
//-------------------------------------

struct MembraneFrame {
    // meta data
    int id = 0;

    // triangular mesh data.
    std::vector< std::array< int, 3 >> triangles;

    // Additional vertex float attributes.
    // Must be of size numAttributes * numVertices.
    Eigen::MatrixXf vertexAttributes;
};

struct FilamentFrame {
    // meta data
    int id = 0;
    int type = 0;

    // line segment data
    std::vector< Vec3f > coords;
};

// Used for line-segment-like elements, such as static linkers and motors.
struct LinkerFrame {
    // meta data
    int id = 0;
    int type = 0;
    int subtype = 0;

    // line data
    std::array< Vec3f, 2 > coords;
};

// Used for bubble-like elements, such as AFMs and MTOCs.
struct BubbleFrame {
    // meta data
    int id = 0;
    int type = 0;

    // point data
    Vec3f coord {};
    float radius = 0;
};


struct DisplayFrame {
    struct CompartmentInfo {
        Vec< 3, int > number { 2, 2, 2 };
        Vec3f         size { 500.0f, 500.0f, 500.0f };
        Vec3f         offset {};
    };

    // meta data
    int serial = 0;
    double simulationTime = 0.0;

    // frame data
    std::vector< MembraneFrame > membranes;
    std::vector< FilamentFrame > filaments;
    std::vector< LinkerFrame >   linkers;
    std::vector< BubbleFrame >   bubbles;

    // other auxiliary optional data
    std::optional< CompartmentInfo > compartmentInfo;
};

struct DisplayTypeMap {
    // linker type map
    std::vector< std::string > linkerTypeName;
    std::unordered_map< std::string, int > linkerTypeMap;

    // Membrane attribute names.
    std::vector< std::string > membraneVertexAttributeNames;
};

struct DisplayData {
    // meta data
    DisplayTypeMap displayTypeMap;

    // all frames
    std::vector< DisplayFrame > frames;

    // Vectorized data.
    //----------------------------------
    // Energies.
    std::vector< std::string > energyNames;
    Eigen::MatrixXd energyValues; // [energy-id][frame]
};


enum class FrameDataOutputFileType {
    unknown,
    traj,
    hdf5,
};
struct DisplayTrajectoryFileSettings {
    std::filesystem::path trajSnapshot = "./snapshot.traj";
    std::filesystem::path lastTrajPath = std::filesystem::current_path();

    FrameDataOutputFileType trajSnapshotFileType() const {
        if (trajSnapshot.extension() == ".traj") {
            return FrameDataOutputFileType::traj;
        } else if (trajSnapshot.extension() == ".h5") {
            return FrameDataOutputFileType::hdf5;
        } else {
            return FrameDataOutputFileType::unknown;
        }
    }
};

//-------------------------------------
// Functions
//-------------------------------------

// Get the displayed linker type from type name.
// If the name does not exist in the map, add a new one.
inline int displayLinkerType(DisplayTypeMap& displayTypeMap, const std::string& name) {
    if(auto it = displayTypeMap.linkerTypeMap.find(name); it == displayTypeMap.linkerTypeMap.end()) {
        // The name was not found
        const int res = displayTypeMap.linkerTypeName.size();
        displayTypeMap.linkerTypeMap.insert({ name, res });
        displayTypeMap.linkerTypeName.push_back(name);
        return res;
    }
    else {
        // The name was found
        return it->second;
    }
}

inline int numFrames(const DisplayData& displayData) { return displayData.frames.size(); }

// Given an open file stream, read one frame of data for display.
//
// Inputs:
// - outMeta: the meta data for the entire output.
// - outSnapshot: the parsed snapshot output read from file.
// - displayTypeMap: the id-name maps of various types.
inline DisplayFrame readOneFrameDataFromOutput(
    const OutputStructMeta&     outMeta,
    const OutputStructSnapshot& outSnapshot,
    DisplayTypeMap&             displayTypeMap
) {
    using namespace std;

    DisplayFrame res;

    // frame number
    res.serial = outSnapshot.snapshot;

    // meta data
    res.simulationTime = outSnapshot.simulationTime;

    // membrane
    {
        // Meta: set membrane attribute names if not set.
        if(displayTypeMap.membraneVertexAttributeNames.empty()) {
            auto& attrnames = displayTypeMap.membraneVertexAttributeNames;
            const Size naf64 = outMeta.membraneMeta.vertexColumnNamesFloat64.size();
            const Size nai64 = outMeta.membraneMeta.vertexColumnNamesInt64.size();
            for(Index i = 0; i < naf64; ++i) {
                attrnames.push_back(outMeta.membraneMeta.vertexColumnNamesFloat64[i]);
            }
            for(Index i = 0; i < nai64; ++i) {
                attrnames.push_back(outMeta.membraneMeta.vertexColumnNamesInt64[i]);
            }
        }

        const int numMembranes = outSnapshot.membraneStruct.size();
        res.membranes.reserve(numMembranes);

        for(const auto& m : outSnapshot.membraneStruct) {
            MembraneFrame mf;

            // Get triangle list.
            mf.triangles.reserve(m.getNumTriangles());
            for(Index ti = 0; ti < m.getNumTriangles(); ++ti) {
                mf.triangles.push_back({
                    static_cast<int>(m.triangleDataInt64(0, ti)),
                    static_cast<int>(m.triangleDataInt64(1, ti)),
                    static_cast<int>(m.triangleDataInt64(2, ti)),
                });
            }

            // Get vertex attributes.
            {
                const Size naf64 = m.vertexDataFloat64.rows();
                const Size nai64 = m.vertexDataInt64.rows();
                mf.vertexAttributes.resize(naf64 + nai64, m.getNumVertices());
                if(naf64 > 0) mf.vertexAttributes.block(0,     0, naf64, m.getNumVertices()) = m.vertexDataFloat64.cast<float>();
                if(nai64 > 0) mf.vertexAttributes.block(naf64, 0, nai64, m.getNumVertices()) = m.vertexDataInt64.cast<float>();

                // Validate.
                if(naf64 + nai64 != displayTypeMap.membraneVertexAttributeNames.size()) {
                    log::error("Membrane vertex attribute names do not match actual data.");
                    throw std::runtime_error("Membrane vertex attribute names do not match actual data.");
                }
            }


            res.membranes.push_back(move(mf));
        }
    }

    // filament
    {
        const int numFilaments = outSnapshot.filamentStruct.size();
        res.filaments.reserve(numFilaments);

        for(const auto& f : outSnapshot.filamentStruct) {
            FilamentFrame ff;

            // get id and type
            ff.id = f.getId();
            ff.type = f.getType();

            // get coordinates
            if(outMeta.filamentMeta.globalFilamentModel == FilamentModel::beadCylinder) {
                const auto numBeads = f.rawCoords.cols();
                ff.coords.reserve(numBeads);
                for(int i = 0; i < numBeads; ++i) {
                    ff.coords.push_back(Vec3f {
                        static_cast<float>(f.rawCoords(0, i)),
                        static_cast<float>(f.rawCoords(1, i)),
                        static_cast<float>(f.rawCoords(2, i)),
                    });
                }
            }

            res.filaments.push_back(move(ff));
        }
    }

    // linkers
    {
        const int numLinkers = outSnapshot.linkerStruct.size();
        res.linkers.reserve(numLinkers);

        for(const auto& l : outSnapshot.linkerStruct) {
            LinkerFrame lf;

            lf.type = displayLinkerType(displayTypeMap, l.type);
            lf.subtype = l.subtype;
            lf.coords[0][0] = l.rawCoords(0, 0);
            lf.coords[0][1] = l.rawCoords(1, 0);
            lf.coords[0][2] = l.rawCoords(2, 0);
            lf.coords[1][0] = l.rawCoords(0, 1);
            lf.coords[1][1] = l.rawCoords(1, 1);
            lf.coords[1][2] = l.rawCoords(2, 1);

            res.linkers.push_back(move(lf));
        }
    }

    // Bubbles.
    {
        const int numBubbles = outSnapshot.bubbleStruct.size();
        res.bubbles.reserve(numBubbles);

        for(const auto& b : outSnapshot.bubbleStruct) {
            BubbleFrame bf;

            bf.id = b.id;
            bf.type = b.type;

            bf.coord[0] = b.coords[0];
            bf.coord[1] = b.coords[1];
            bf.coord[2] = b.coords[2];
            bf.radius = b.radius;

            res.bubbles.push_back(move(bf));
        }
    }

    // Compartment grid.
    if(outMeta.simulConfig.has_value()) {
        auto& geoParams = outMeta.simulConfig->geoParams;

        res.compartmentInfo.emplace();

        res.compartmentInfo->number = {
            geoParams.NX,
            geoParams.NY,
            geoParams.NZ,
        };
        res.compartmentInfo->size = {
            static_cast< float >(geoParams.compartmentSizeX),
            static_cast< float >(geoParams.compartmentSizeY),
            static_cast< float >(geoParams.compartmentSizeZ),
        };
    }

    return res;
}

inline DisplayData readAllFrameDataFromOutput(
    const DisplayTrajectoryFileSettings& inputs
) {
    using namespace std;

    DisplayData res;

    const auto fileType = inputs.trajSnapshotFileType();
    if(fileType == FrameDataOutputFileType::traj) {
        // Read snapshot
        ifstream is(inputs.trajSnapshot);

        int curFrame = 0;

        log::info("Start reading {}", inputs.trajSnapshot);

        OutputStructMeta outMeta;
        outMeta.membraneMeta.vertexColumnNamesFloat64 = { "coord.x", "coord.y", "coord.z" };

        string line;
        while(true) {
            getline(is, line);
            if(!is) break;
            if(line.empty()) continue;

            ++curFrame;
            if (curFrame % 20 == 0) log::info("Frame {}", curFrame);

            istringstream iss(line);
            OutputStructSnapshot outSnapshot;
            outSnapshot.snapshot = curFrame;
            outSnapshot.getFromOutput(is, iss);
            res.frames.push_back(readOneFrameDataFromOutput(outMeta, outSnapshot, res.displayTypeMap));
        }
        log::info("Reading complete. {} frames loaded.", numFrames(res));

    }
    else if(fileType == FrameDataOutputFileType::hdf5) {
        // Read snapshot.
        h5::File file(inputs.trajSnapshot.string(), h5::File::ReadOnly);
        h5::Group groupHeader = file.getGroup("/header");
        h5::Group groupSnapshots = file.getGroup("/snapshots");

        // Read number of frames.
        std::int64_t numFrames = 0;
        h5::readDataSet(numFrames, groupHeader, "count");
        log::info("Start reading {}, {} frames to be loaded.", inputs.trajSnapshot, numFrames);
        // Read meta data.
        OutputStructMeta meta;
        read(meta, groupHeader);
        {
            // Energies.
            res.energyNames = meta.energyNames;
            res.energyValues.resize(numFrames, res.energyNames.size());
        }

        // Read all frames.
        for(Index curFrame = 0; curFrame < numFrames; ++curFrame) {
            if(curFrame % 20 == 0) log::info("Loading frame {}", curFrame);
            OutputStructSnapshot outSnapshot;
            read(outSnapshot, groupSnapshots, curFrame);
            res.frames.push_back(readOneFrameDataFromOutput(meta, outSnapshot, res.displayTypeMap));

            // Energies.
            for(Index i = 0; i < res.energyNames.size(); ++i) {
                res.energyValues(curFrame, i) = outSnapshot.energies[i];
            }
        }
        log::info("Reading complete. {} frames loaded.", numFrames);
    }
    else {
        log::error("Unknown file extension for {}", inputs.trajSnapshot);
        throw runtime_error("Unknown file type.");
    }


    return res;
}


} // namespace medyan::visual

#endif

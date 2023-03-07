#ifndef MEDYAN_Structure_OutputStruct_hpp
#define MEDYAN_Structure_OutputStruct_hpp

#include <array>
#include <iostream>
#include <optional>
#include <sstream>
#include <tuple>
#include <vector>

#include <Eigen/Dense>
#include <xtensor/xtensor.hpp>

#include "Structure/FilamentTypes.hpp"
#include "Util/Math/Vec.hpp"
#include "Util/Parser/StringParser.hpp"
#include "SysParams.h"

namespace medyan {

/*!
 * The OutputStruct class serves as a connection between the real structure in
 * MEDYAN and the output structure. It may contain minimalistic data for
 * converting from system data and converting from/to output literals.
 * 
 * How much required information of each column can be obtained from each row
 * 
 *              | System Data   | Stored Data   | Output Data
 * -------------|---------------|---------------|---------------
 * System Data  | --            | Full          | Full
 * Stored Data  | None          | --            | Full
 * Output Data  | None          | Full          | --
 * 
 * (Except for snapshot serial which is not contained in the system)
 * 
 */


struct OutputStructFilament {
    static constexpr char name[] = "FILAMENT";

    void outputFromStored(std::ostream& os);
    void getFromOutput(std::istream& is, std::istringstream& iss);

    auto getId() const { return id; }
    auto getType() const { return type; }
    int getNumBeads()const { return numBeads; }

    /// Data
    int id = 0;
    int type = 0;
    int numBeads = 0;
    short deltaMinusEnd = 0;
    short deltaPlusEnd = 0;

    Eigen::MatrixXd rawCoords;

};
struct OutputStructFilamentMeta {
    FilamentModel globalFilamentModel = FilamentModel::beadCylinder;
};

struct OutputStructLinker {

    void outputFromStored(std::string_view name, std::ostream& os);
    void getFromOutput(std::string_view type, std::istream& is, std::istringstream& iss);

    auto getType() const { return subtype; }

    /// Data
    std::string type;
    int id = 0;
    int subtype = 0;

    Eigen::Matrix<double, 3, 2> rawCoords;
};


struct OutputStructBubble {
    static constexpr char name[] = "BUBBLE";

    void outputFromStored(std::ostream& os);
    void getFromOutput(std::istream& is, std::istringstream& iss);

    /// Data
    int id = 0;
    int type = 0;

    Eigen::Vector3d coords;
    double radius = 0;

};

struct OutputStructMembrane {

    static constexpr char name[] = "MEMBRANE";

    void outputFromStored(std::ostream& os);
    void getFromOutput(std::istream& is, std::istringstream& iss);

    auto getNumVertices()const { return numVertices; }
    auto getNumTriangles()const { return numTriangles; }

    /// Data
    int type = 0;
    Size numVertices = 0;
    Size numTriangles = 0;
    Size numBorders = 0;

    // The vertex data is stored as a matrix of size (num attributes, num vertices).
    // The columns are specified by the VertexColumnSpec structs.
    Eigen::MatrixXd vertexDataFloat64;
    Eigen::MatrixXi vertexDataInt64;

    // Triangle data.
    Eigen::MatrixXi triangleDataInt64;

    // Packed border vertices.
    // Let m be the number of borders.
    // Let n_i be the number of vertices in border i, for i = 1, 2, ..., m.
    // Let b_ij be the jth vertex of border i, for j = 1, 2, ..., n_i.
    // Then the packed border vertices are stored as an array as follows:
    //   [
    //     n_1,  b_11, b_12, ..., b_{1 n_1},
    //     n_2,  b_21, b_22, ..., b_{2 n_2},
    //     ...
    //     n_m,  b_m1, b_m2, ..., b_{m n_m}
    //   ]
    std::vector< std::int64_t > packedBorderVertices;
};
struct OutputStructMembraneMeta {

    // Vertex attributes are stored in the following order:
    // [
    //   coord.x, coord.y, coord.z,
    //   area., curvature.curv,
    //   conc.<membrane-diffusing-species-name>...,
    // ]
    std::vector<std::string>             vertexColumnNamesFloat64;
    // Vertex attributes are stored in the following order:
    // [
    //   copyNumber.<membrane-diffusing-species-name>...,
    // ]
    std::vector<std::string>             vertexColumnNamesInt64;
};


struct OutputStructChemistry {
    // Diffusing species.
    // - A buffer storing column major 4D data, indexed by (diffusing-species-index, comp-i, comp-j, comp-k).
    xt::xtensor<Size, 4, xt::layout_type::column_major> diffusingSpeciesCopyNumbers;

    // Global species.
    // - The size must be equal to the global species names in the meta data.
    std::vector<Size> globalSpeciesCopyNumbers;
};


struct OutputStructSnapshot {
    /// Bonus data (generated when constructed)
    Index snapshot;

    // data
    double simulationTime = 0;

    // Elements.
    std::vector<OutputStructFilament> filamentStruct;
    std::vector<OutputStructLinker>   linkerStruct;
    std::vector<OutputStructBubble>   bubbleStruct;
    std::vector<OutputStructMembrane> membraneStruct;

    // Chemistry.
    OutputStructChemistry chemistry;

    // Energies.
    // In MEDYAN units, and should correspond to the sequence of energies in the meta data.
    std::vector< double > energies;

    /// Non data

    OutputStructSnapshot(int snapshot = 0): snapshot(snapshot) {}

    void outputFromStored(std::ostream& os);
    void outputFromStoredWithoutChildren(std::ostream& os);
    void getFromOutput(std::istream& is, std::istringstream& iss);
};

struct OutputStructMeta {
    OutputStructMembraneMeta membraneMeta;
    OutputStructFilamentMeta filamentMeta;

    //----------------------------------
    // Other system meta data.
    //----------------------------------

    // Entire simulation configuration.
    std::optional<SimulConfig> simulConfig;

    // Energies.
    std::vector< std::string > energyNames;

    // Diffusing species names.
    std::vector< std::string > diffusingSpeciesNames;

    // Names of species that are reported globally. (May include diffusing species).
    std::vector< std::string > globalSpeciesNames;
};

} // namespace medyan

#endif

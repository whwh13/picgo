
#include <iterator>
#include <string>

#include "common.h"
#include "Structure/OutputStruct.hpp"

namespace medyan {

/******************************************************************************
OutputStruct for Filaments
******************************************************************************/
//@{
constexpr char OutputStructFilament::name[];

void OutputStructFilament::outputFromStored(std::ostream& os) {
    os << name << " "
        << id << " "
        << type << " "
        << numBeads << " "
        << deltaMinusEnd << " "
        << deltaPlusEnd << std::endl;

    for(int i = 0; i < numBeads; ++i) {
        for(int dim = 0; dim < 3; ++dim) {
            os << rawCoords(dim, i) << ' ';
        }
    }
    os << std::endl;
}

void OutputStructFilament::getFromOutput(std::istream& is, std::istringstream& iss) {
    iss >> id
        >> type
        >> numBeads
        >> deltaMinusEnd
        >> deltaPlusEnd;

    std::string nextLine;
    std::getline(is, nextLine);
    std::istringstream newIss(nextLine);

    rawCoords.resize(3, numBeads);
    for(int i = 0; i < numBeads; ++i) {
        for(int dim = 0; dim < 3; ++dim) {
            newIss >> rawCoords(dim, i);
        }
    }
}
//@}

/******************************************************************************
OutputStruct for Linkers
******************************************************************************/

void OutputStructLinker::outputFromStored(std::string_view name, std::ostream& os) {
    os << name << " "
        << id << " "
        << subtype << std::endl;

    for(int n = 0; n < 2; ++n) {
        for(int dim = 0; dim < 3; ++dim) {
            os << rawCoords(dim, n) << ' ';
        }
    }
    os << std::endl;
}

void OutputStructLinker::getFromOutput(std::string_view type, std::istream& is, std::istringstream& iss) {
    this->type = type;
    iss >> id
        >> subtype;

    std::string nextLine;
    std::getline(is, nextLine);
    std::istringstream newIss(nextLine);
    for(int n = 0; n < 2; ++n) {
        for(int dim = 0; dim < 3; ++dim) {
            newIss >> rawCoords(dim, n);
        }
    }
}


/******************************************************************************
OutputStruct for Bubbles
******************************************************************************/
//@{
constexpr char OutputStructBubble::name[];

void OutputStructBubble::outputFromStored(std::ostream& os) {
    os << name << " "
        << id << " "
        << type << std::endl;

    os
        << coords[0] << ' '
        << coords[1] << ' '
        << coords[2] << ' '
        << std::endl;
}

void OutputStructBubble::getFromOutput(std::istream& is, std::istringstream& iss) {
    iss >> id
        >> type;

    std::string nextLine;
    std::getline(is, nextLine);
    std::istringstream newIss(nextLine);
    newIss >> coords[0] >> coords[1] >> coords[2];
}
//@}

/******************************************************************************
OutputStruct for Membranes
******************************************************************************/
//@{
constexpr char OutputStructMembrane::name[];

void OutputStructMembrane::outputFromStored(std::ostream& os) {
    os << name << " "
        << type << " "
        << numVertices << ' '
        << numTriangles << ' '
        << numBorders << '\n';

    // print coordinates.
    for(int i = 0; i < numVertices; ++i) {
        for(int dim = 0; dim < 3; ++dim) {
            os << vertexDataFloat64(dim, i) << ' ';
        }
        os << '\n';
    }

    // print triangles.
    for(int i = 0; i < numTriangles; ++i) {
        for(int j = 0; j < 3; ++j) {
            os << triangleDataInt64(j, i) << ' ';
        }
        os << '\n';
    }

    // print borders.
    Index curBorderCountIndex = 0;
    for(Index i = 0; i < numBorders; ++i) {
        // Number of vertices in the border.
        auto nv = packedBorderVertices[curBorderCountIndex];
        os << nv << ' ';

        // Print vertex indices.
        for(Index vi = 0; vi < nv; ++vi) {
            os << packedBorderVertices[curBorderCountIndex + vi + 1] << ' ';
        }

        // Move to the next border.
        os << '\n';
        curBorderCountIndex += nv + 1;
    }
}

void OutputStructMembrane::getFromOutput(std::istream& is, std::istringstream& iss) {
    iss >> type
        >> numVertices
        >> numTriangles
        >> numBorders;

    vertexDataFloat64.resize(3, numVertices);
    for(Index i = 0; i < numVertices; ++i) {
        std::string nextLine;
        std::getline(is, nextLine);
        std::istringstream newIss(nextLine);

        // Record coordinates.
        for(int dim = 0; dim < 3; ++dim) {
            newIss >> vertexDataFloat64(dim, i);
        }
    }

    triangleDataInt64.resize(3, numTriangles);
    for(Index i = 0; i < numTriangles; ++i) {
        std::string nextLine;
        std::getline(is, nextLine);
        std::istringstream newIss(nextLine);

        // Record indices.
        for(int j = 0; j < 3; ++j) {
            newIss >> triangleDataInt64(j, i);
        }
    }

    packedBorderVertices.clear();
    for(Index i = 0; i < numBorders; ++i) {
        std::string nextLine;
        std::getline(is, nextLine);
        std::istringstream newIss(nextLine);

        // First number is the number of vertices on this border.
        Size nv; newIss >> nv;

        packedBorderVertices.push_back(nv);

        for(Index vi = 0; vi < nv; ++vi) {
            Index v; newIss >> v;
            packedBorderVertices.push_back(v);
        }
    }
}

//@}

/******************************************************************************
OutputStruct for snapshots
******************************************************************************/

void OutputStructSnapshot::outputFromStored(std::ostream& os) {
    outputFromStoredWithoutChildren(os);
    
    for(auto& it: filamentStruct) {
        it.outputFromStored(os);
    }
    for(auto& it: linkerStruct) {
        it.outputFromStored(
            (
                it.type == "linker" ? "LINKER" :
                it.type == "motor"  ? "MOTOR" :
                it.type == "brancher" ? "BRANCHER" :
                "ERROR"
            ),
            os
        );
    }
    for(auto& it: bubbleStruct) {
        it.outputFromStored(os);
    }
    for(auto& it: membraneStruct) {
        it.outputFromStored(os);
    }

    // Note: new children should be added here
}

void OutputStructSnapshot::outputFromStoredWithoutChildren(std::ostream& os) {
    const auto countLinkersByType = [this](std::string_view type) {
        int count = 0;
        for(auto& it: linkerStruct) {
            if(it.type == type) {
                ++count;
            }
        }
        return count;
    };
    Size legacyNumLinkers   = countLinkersByType("linker");
    Size legacyNumMotors    = countLinkersByType("motor");
    Size legacyNumBranchers = countLinkersByType("brancher");
    os << snapshot << " "
        << simulationTime << " "
        << filamentStruct.size() << " "
        << legacyNumLinkers << " "
        << legacyNumMotors << " "
        << legacyNumBranchers << " "
        << bubbleStruct.size() << " "
        << membraneStruct.size() << std::endl;
}

void OutputStructSnapshot::getFromOutput(std::istream& is, std::istringstream& iss) {
    Size numFilaments, numLinkers, numMotorGhosts, numBranchingPoints, numBubbles, numMembranes;
    iss >> snapshot
        >> simulationTime
        >> numFilaments
        >> numLinkers
        >> numMotorGhosts
        >> numBranchingPoints
        >> numBubbles
        >> numMembranes;

    std::string nextLine;
    do {
        std::getline(is, nextLine);
        if(!is || nextLine.empty()) break;

        std::istringstream newIss(nextLine);
        std::string name;
        newIss >> name;
        if(name == OutputStructFilament::name) {
            filamentStruct.emplace_back();
            filamentStruct.back().getFromOutput(is, newIss);
        } else if(name == "LINKER") {
            linkerStruct.emplace_back();
            linkerStruct.back().getFromOutput("linker", is, newIss);
        } else if(name == "MOTOR") {
            linkerStruct.emplace_back();
            linkerStruct.back().getFromOutput("motor", is, newIss);
        } else if(name == "BRANCHER") {
            linkerStruct.emplace_back();
            linkerStruct.back().getFromOutput("brancher", is, newIss);
        } else if(name == OutputStructBubble::name) {
            bubbleStruct.emplace_back();
            bubbleStruct.back().getFromOutput(is, newIss);
        } else if(name == OutputStructMembrane::name) {
            membraneStruct.emplace_back();
            membraneStruct.back().getFromOutput(is, newIss);
        } else {
            log::error("Unrecognized output structure name: {}", name);
            throw std::runtime_error("Unrecognized output structure name.");
        }

    } while(true);
}
//@}

} // namespace medyan

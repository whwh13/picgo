#include "ReadSnapshot.hpp"

#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <vector>

#include "utility.h"


#include "Analysis/Io/pdb.h"
#include "Util/Io/Log.hpp"
#include "Structure/OutputStruct.hpp"
#include "Structure/SubSystem.h"
#include "SysParams.h"

namespace medyan {
namespace analysis {

namespace {
    /// Helper struct to find number of beads (atoms) required for pdb-style output
    struct PdbMaxBead {
        std::vector<size_t> filament; // Maximum beads in each filament
        size_t linker = 0;
        size_t motor = 0;
        size_t membrane = 0;

        size_t maxBead = 0;

        void renew(const OutputStructSnapshot& snapshot) {
            size_t curMaxBead = 0;

            // Beads in filaments
            for(auto& eachFilament: snapshot.filamentStruct) {
                int curId = eachFilament.getId();
                if(curId >= filament.size()) {
                    filament.resize(curId + 1, 0);
                }
                if(filament[curId] < eachFilament.getNumBeads()) filament[curId] = eachFilament.getNumBeads();
                curMaxBead += filament[curId];
            }

            // Beads in linkers
            size_t newLinker = snapshot.linkerStruct.size() * 2;
            if(newLinker > linker) linker = newLinker;
            curMaxBead += linker;

            // Beads in membranes
            size_t newMembrane = 0;
            for(auto& eachMembrane: snapshot.membraneStruct) {
                /*
                newMembrane += eachMembrane.getNumEdges() * 2;
                */
                newMembrane += eachMembrane.getNumVertices();
            }
            if(newMembrane > membrane) membrane = newMembrane;
            curMaxBead += membrane;

            // Total beads
            if(curMaxBead > maxBead) maxBead = curMaxBead;
        }
    };

}

void SnapshotReader::readAndConvertToVmd(const size_t maxFrames, const RunAnalyzeParams& params) {
    using namespace std;

    vector<OutputStructSnapshot> snapshots;

    // Read snapshot
    ifstream is(_snapshotFilepath);

    PdbMaxBead maxBead;

    size_t curFrame = 0;
    const auto frameInterval = params.analyzeFrameInterval;

    LOG(STEP) << "Start reading " << _snapshotFilepath;
    if(frameInterval != 1) LOG(INFO) << "Frame interval is " << frameInterval;

    string line;
    while(maxFrames == 0 || curFrame < maxFrames) {
        getline(is, line);
        if(!is) break;
        if(line.empty()) continue;

        ++curFrame;
        if (curFrame % 20 == 0) LOG(INFO) << "Frame " << curFrame;

        istringstream iss(line);
        snapshots.emplace_back(0);
        snapshots.back().getFromOutput(is, iss);

        // Skip frames
        if(curFrame % frameInterval) {
            snapshots.pop_back();
            continue;
        }

        maxBead.renew(snapshots.back());
    }
    LOG(STEP) << "Reading complete. " << snapshots.size() << " frames to be processed.";

    // close input stream
    is.close();

    // Write to pdb
    PdbGenerator pg(_pdbFilepath);

    const size_t bondFrame = params.analyzeMembraneBondFrame;
    const bool   allBonds  = params.analyzeMembraneBondAllFrames;

    size_t numSnapshots = snapshots.size();
    for(size_t idx = 0; idx < numSnapshots; ++idx) {
        if (idx % 20 == 0) LOG(INFO) << "Generating model " << idx;
        pg.genModel(idx + 1);

        // Psf generation for fixed topology
        const bool shouldMakePsf = (allBonds || bondFrame == idx);
        std::unique_ptr< PsfGenerator > psfGen;
        if(shouldMakePsf) {
            const auto psfFilename = [&]() {
                std::ostringstream oss;
                oss << psfFileDir_ << '/' << psfFilenameMain_ << '-'
                    << std::setfill('0') << std::setw(6) << idx
                    << ".psf";
                return oss.str();
            }();
            psfGen = std::make_unique< PsfGenerator >(psfFilename);
            psfGen->genHeader();
            psfGen->genNatom(maxBead.maxBead);
        }

        char chain;
        size_t atomSerial; // Pdb file atom number
        size_t atomId = 0; // Psf file atom number. TER does not increase this variable.
        size_t atomCount; // Temporary counter
        size_t resSeq;

        SubSystem s; // Dummy subsystem

        // Filaments
        chain = 'F';
        atomSerial = 0;
        resSeq = 0;
        auto filamentPtr = snapshots[idx].filamentStruct.begin();
        size_t filamentAlloc = maxBead.filament.size();

        std::vector< std::pair< size_t, size_t >> bondList;

        for(size_t i = 0; i < filamentAlloc; ++i) {
            atomCount = 0;

            if(filamentPtr != snapshots[idx].filamentStruct.end() && filamentPtr->getId() == i) {
                // Filament exists for this id.
                const auto numBeads = filamentPtr->getNumBeads();
                auto& allCoords = filamentPtr->rawCoords;
                for(int beadIndex = 0; beadIndex < numBeads; ++beadIndex) {
                    auto eachCoord = allCoords.col(beadIndex);
                    ++atomSerial;
                    ++atomId;
                    ++atomCount;
                    ++resSeq;
                    pg.genAtom(
                        atomSerial, " CA ", ' ', "ARG", chain, resSeq, ' ',
                        eachCoord[0], eachCoord[1], eachCoord[2]
                    );
                    if(shouldMakePsf) {
                        psfGen->genAtom(atomId, "", resSeq, "ARG", "CA", "CA");
                    }
                }
                while(atomCount < maxBead.filament[i]) {
                    ++atomSerial;
                    ++atomId;
                    ++atomCount;
                    ++resSeq;
                    pg.genAtom(
                        atomSerial, " CA ", ' ', "ARG", chain, resSeq, ' ',
                        allCoords(0, numBeads - 1), allCoords(1, numBeads - 1), allCoords(2, numBeads - 1)
                    );
                    if(shouldMakePsf) {
                        psfGen->genAtom(atomId, "", resSeq, "ARG", "CA", "CA");
                    }
                }

                ++filamentPtr;
            } else {
                // Filament does not exist for this id.
                while(atomCount < maxBead.filament[i]) {
                    ++atomSerial;
                    ++atomId;
                    ++atomCount;
                    ++resSeq;
                    pg.genAtom(
                        atomSerial, " CA ", ' ', "ARG", chain, resSeq
                    );
                    if(shouldMakePsf) {
                        psfGen->genAtom(atomId, "", resSeq, "ARG", "CA", "CA");
                    }
                }
            }

            ++resSeq; // Separate trace.
        }
        pg.genTer(++atomSerial, "ARG", chain, resSeq);

        // Linkers
        atomCount = 0;
        resSeq = 0;
        for(auto& eachLinker: snapshots[idx].linkerStruct) {
            chain = eachLinker.type == "motor" ? 'M' :
                    eachLinker.type == "brancher" ? 'B' :
                    'L';
            for(int n = 0; n < 2; ++n) {
                auto eachCoord = eachLinker.rawCoords.col(n);
                ++atomSerial;
                ++atomId;
                ++atomCount;
                ++resSeq;
                pg.genAtom(
                    atomSerial, " CA ", ' ', "ARG", chain, resSeq, ' ',
                    eachCoord[0], eachCoord[1], eachCoord[2]
                );
                if(shouldMakePsf) {
                    psfGen->genAtom(atomId, "", resSeq, "ARG", "CA", "CA");
                }
            }
            ++resSeq; // Separate bonds
        }
        while(atomCount < maxBead.linker) {
            for(size_t b = 0; b < 2; ++b) {
                ++atomSerial;
                ++atomId;
                ++atomCount;
                ++resSeq;
                pg.genAtom(
                    atomSerial, " CA ", ' ', "ARG", chain, resSeq
                );
                if(shouldMakePsf) {
                    psfGen->genAtom(atomId, "", resSeq, "ARG", "CA", "CA");
                }
            }
            ++resSeq; // Separate bonds
        }
        pg.genTer(++atomSerial, "ARG", chain, resSeq);

        // Membranes
        chain = 'E';
        atomCount = 0;
        resSeq = 0;

        for(auto& eachMembrane: snapshots[idx].membraneStruct) {

            size_t atomIdSoFar = atomId;

            for(Index vi = 0; vi < eachMembrane.numVertices; ++vi) {
                auto v = eachMembrane.vertexDataFloat64.col(vi);
                ++atomSerial;
                ++atomId;
                ++atomCount;
                ++resSeq;
                pg.genAtom(
                    atomSerial, " CA ", ' ', "ARG", chain, 1/* resSeq */, ' ',
                    v[0], v[1], v[2]
                );
                if(shouldMakePsf) {
                    psfGen->genAtom(atomId, "", 1/* resSeq */, "ARG", "CA", "CA");
                }
            }
            ++resSeq;

            if(shouldMakePsf) {
                for(Index ti = 0; ti < eachMembrane.numTriangles; ++ti) {
                    auto t = eachMembrane.triangleDataInt64.col(ti);
                    for(size_t i = 0; i < 3; ++i) {
                        size_t i_next = (i+1) % 3;
                        if(t[i] < t[i_next]) {
                            bondList.emplace_back(
                                atomIdSoFar + t[i] + 1,
                                atomIdSoFar + t[i_next] + 1
                            );
                        }
                    }
                }

                Index curBorderCountIndex = 0;
                for(Index bi = 0; bi < eachMembrane.numBorders; ++bi) {
                    Size nv = eachMembrane.packedBorderVertices[curBorderCountIndex];
                    for(int i = 0; i < nv; ++i) {
                        auto i_next = (i+1) % nv;
                        auto v1 = eachMembrane.packedBorderVertices[curBorderCountIndex + 1 + i];
                        auto v2 = eachMembrane.packedBorderVertices[curBorderCountIndex + 1 + i_next];
                        if(v1 < v2) {
                            bondList.emplace_back(
                                atomIdSoFar + v1 + 1,
                                atomIdSoFar + v2 + 1
                            );
                        }
                    }
                    curBorderCountIndex += 1 + nv;
                }
            }

        }

        while(atomCount < maxBead.membrane) {
            ++atomSerial;
            ++atomId;
            ++atomCount;
            ++resSeq;
            pg.genAtom(
                atomSerial, " CA ", ' ', "ARG", chain, resSeq
            );
            if(shouldMakePsf) {
                psfGen->genAtom(atomId, "", resSeq, "ARG", "CA", "CA");
            }
        }

        // Generating bond
        if(shouldMakePsf) {
            size_t numBonds = bondList.size();
            psfGen->genNbond(numBonds);
            psfGen->genBondStart();

            for(const auto& bond : bondList) {
                psfGen->genBond(bond.first, bond.second);
            }

            psfGen->genBondEnd();
            LOG(INFO) << "Bond info generated.";
        }

        pg.genTer(++atomSerial, "ARG", chain, resSeq);

        // End of model
        pg.genEndmdl();

    } // End of looping through snapshots
    LOG(STEP) << "Writing complete. " << numSnapshots << " models created.";

}

} // namespace analysis
} // namespace medyan

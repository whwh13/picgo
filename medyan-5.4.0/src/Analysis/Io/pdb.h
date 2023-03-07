#ifndef MEDYAN_ANALYSIS_IO_PDB_H
#define MEDYAN_ANALYSIS_IO_PDB_H

#include <iostream>
#include <filesystem>
#include <fstream>
#include <string>

namespace medyan {
namespace analysis {

// Checkout www.wwpdb.org/

class PdbGenerator {
    std::filesystem::path _pdbFilepath;
    std::ofstream _ofs; ///< Handles the file output

public:
    PdbGenerator(const std::filesystem::path& filepath): _pdbFilepath(filepath), _ofs(filepath) {}

    void genModel(int serial);
    void genEndmdl();
    void genAtom(int serial, std::string name, char altLoc, std::string resName, char chainID,
        int resSeq, char iCode=' ', double x=0.0, double y=0.0, double z=0.0, double occupancy=1.0, double tempFactor=0.0,
        std::string element="  ", std::string charge="  ");
    void genTer(int serial, std::string resName, char chainID,
        int resSeq, char iCode=' ');

};

// Psf contains the bond information
// Uses implementation from www.ks.uiuc.edu/Research/vmd/plugins/doxygen/psfplugin_8c-source.html
class PsfGenerator {
    std::filesystem::path _psfFilepath;
    std::ofstream _ofs; ///< Handles the file output

    // State machine for genBond (0-4)
    int _bondOutputPos = 0;

public:
    PsfGenerator(const std::filesystem::path& filepath): _psfFilepath(filepath), _ofs(filepath) {}

    void genHeader();
    void genNatom(int num);
    void genAtom(int id, std::string segName, int resId, std::string resName, std::string name,
        std::string type, double charge=0.0, double mass=12.0110); ///< Starting from id 1. TER ignored.
    void genNbond(int num);
    void genBond(int id1, int id2);

    // Helper functions for genBond
    void genBondStart() { _bondOutputPos = 0; }
    void genBondEnd();
};


} // namespace analysis
} // namespace medyan

#endif
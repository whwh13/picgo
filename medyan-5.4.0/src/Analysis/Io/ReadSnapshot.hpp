#ifndef MEDYAN_ANALYSIS_IO_ReadSnapshot_hpp
#define MEDYAN_ANALYSIS_IO_ReadSnapshot_hpp

#include <filesystem>

#include "SysParams.h"

namespace medyan {
namespace analysis {

class SnapshotReader {
private:
    std::filesystem::path _snapshotFilepath;
    std::filesystem::path _pdbFilepath;
    std::filesystem::path psfFileDir_;
    std::filesystem::path psfFilenameMain_;

public:
    SnapshotReader(const std::filesystem::path& snapshotFilepath, const std::filesystem::path& pdbFilepath, const std::filesystem::path& psfFileDir, const std::filesystem::path& psfFilenameMain):
        _snapshotFilepath(snapshotFilepath), _pdbFilepath(pdbFilepath), psfFileDir_(psfFileDir), psfFilenameMain_(psfFilenameMain) {}

    void readAndConvertToVmd(const size_t maxFrames, const RunAnalyzeParams& params); // 0 means no limit on max frames
};

} // namespace analysis
} // namespace medyan

#endif
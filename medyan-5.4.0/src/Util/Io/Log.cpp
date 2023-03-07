#include "Util/Io/Log.hpp"

#include <chrono>
#include <ctime>
#include <iomanip>
#include <unordered_map>

#ifdef PLATFORM_UNIX_LIKE
    #include <unistd.h>
#elif defined(PLATFORM_WINDOWS)
    #include <Windows.h>
#endif

namespace medyan {
namespace logger {

namespace {

// Helper functions

std::string timeLiteralGeneration() {
    using namespace std::chrono;

    system_clock::time_point p = system_clock::now();
    milliseconds ms = duration_cast<milliseconds>(p.time_since_epoch());
    seconds s = duration_cast<seconds>(ms);

    std::time_t timeToSec = s.count();
    tm timeinfoToSec;
#ifdef COMPILER_MSVC
    localtime_s(&timeinfoToSec, &timeToSec);
#elif defined(COMPILER_GCC) || defined(COMPILER_CLANG)
    localtime_r(&timeToSec, &timeinfoToSec);
#else
    // Not thread safe
    timeinfoToSec = *localtime(&timeToSec);
#endif

    std::size_t msRemain = ms.count() % 1000;

    std::stringstream ss;
    ss << timeinfoToSec.tm_year + 1900 << '-'
        << std::setfill('0') << std::setw(2) << timeinfoToSec.tm_mon + 1 << '-'
        << std::setw(2) << timeinfoToSec.tm_mday << ' '
        << std::setw(2) << timeinfoToSec.tm_hour << ':'
        << std::setw(2) << timeinfoToSec.tm_min << ':'
        << std::setw(2) << timeinfoToSec.tm_sec << '.'
        << std::setw(3) << msRemain;

    return ss.str();
}

} // namespace


// Level color codes.
const std::unordered_map<LogLevel, const char*> logLevelColorAnsi {
    {LogLevel::Debug,   "\033[37m"}, // White
    {LogLevel::Info,    "\033[97m"}, // Bright white
    {LogLevel::Step,    "\033[96m"}, // Bright cyan
    {LogLevel::Note,    "\033[92m"}, // Bright green
    {LogLevel::Warning, "\033[93m"}, // Bright yellow
    {LogLevel::Error,   "\033[91m"}, // Bright red
    {LogLevel::Fatal,   "\033[91m"}  // Bright red
};
constexpr const char * resetAnsi = "\033[0m";


namespace internal {

void LogWriter::logDispatch() {
    bool genTime = true; // Generate time literal only once
    std::string strTime;

    const auto& settings = _logger.settings;

    for(auto& eachOs: _logger.getOsContainers()) {
        if(eachOs.disp.isOnWith(_lv)) {
            // The oss for final output
            std::ostringstream finalOss;

            // Prefix generation
            if(eachOs.dispColor.isOnWith(_lv)) {
                finalOss << logLevelColorAnsi.find(_lv)->second;
            }
            if(eachOs.dispTime.isOnWith(_lv)) {
                if(genTime) {
                    strTime = timeLiteralGeneration();
                    genTime = false;
                }
                finalOss << settings.delimiterBefore << strTime << settings.delimiterAfter
                    << (settings.spaceAfterDelimiter? " ": "");
            }
            if(eachOs.dispLevel.isOnWith(_lv)) {
                finalOss << settings.delimiterBefore << literal(_lv) << settings.delimiterAfter
                    << (settings.spaceAfterDelimiter? " ": "");
            }

            // Attach log content
            finalOss << oss_.str();

            // Suffix generation
            if(eachOs.dispFile.isOnWith(_lv)) {
                finalOss << (settings.spaceAfterDelimiter? " ": "")
                    << settings.delimiterBefore << "File " << _curFile << settings.delimiterAfter;;
            }
            if(eachOs.dispLine.isOnWith(_lv)) {
                finalOss << (settings.spaceAfterDelimiter? " ": "")
                    << settings.delimiterBefore << "Line " << _curLine << settings.delimiterAfter;
            }
            if(eachOs.dispFunc.isOnWith(_lv)) {
                finalOss << (settings.spaceAfterDelimiter? " ": "")
                    << settings.delimiterBefore << "Function " << _curFunc << settings.delimiterAfter;
            }

            if(eachOs.dispColor.isOnWith(_lv)) {
                finalOss << resetAnsi;
            }
            if(settings.newLineAfterLog) finalOss << '\n';

            // This output should not cause data races
            (*eachOs.os) << finalOss.str();

            // Flush
            if(eachOs.flushLevel.isOnWith(_lv)) (*eachOs.os) << std::flush;
        }
    }

}
} // namespace internal


} // namespace logger
} // namespace medyan

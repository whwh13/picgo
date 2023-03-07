#ifndef MEDYAN_Util_Io_Log_Hpp
#define MEDYAN_Util_Io_Log_Hpp

#include <algorithm> // find_if
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <utility>

#include <spdlog/fmt/ostr.h>    // For formatting types with operator<< defined.
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/spdlog.h>

#include "Util/Environment.hpp"

namespace medyan {
namespace log {

using spdlog::debug;
using spdlog::info;
using spdlog::warn;
using spdlog::error;
using spdlog::critical;

inline void configureLoggers(const std::filesystem::path& outDir) {
    // Console output.
    auto sinkConsole = std::make_shared< spdlog::sinks::stdout_color_sink_mt >();
    sinkConsole->set_level(spdlog::level::info);
    sinkConsole->set_pattern("[%^%l%$] %v");

    // File output.
    const auto logFile = outDir / "medyan.log";
    auto sinkFile = std::make_shared< spdlog::sinks::basic_file_sink_mt >(logFile.string(), true);
    sinkFile->set_level(spdlog::level::trace);

    // Create logger.
    auto logger = std::make_shared< spdlog::logger >("medyan", spdlog::sinks_init_list{sinkConsole, sinkFile});
    logger->set_level(spdlog::level::trace);
    logger->flush_on(spdlog::level::trace);

    // Set as default logger.
    spdlog::set_default_logger(logger);
}

} // namespace log

namespace logger {

/**
 * This file implements a versatile yet easy to use logger.
 * 
 * To use the logger, first configure the output file path via macro
 *     MEDYAN_LOG_DEFAULT_CONFIGURATION(filepath)
 * Then use
 *     LOG(level) << ... << ... ;
 * just as `cout` to invoke logging. Note that `endl` is not needed for default
 * configuration.
 * 
 * The LOG macro actually creates a temporary LogWriter object which will be
 * destroyed at the end of the expression when the log will be printed. The
 * logging is supposed to be thread-safe as long as the underlying `ostream`s
 * are thread-safe.
 * 
 * Note: the logger is designed for easy logging of run-time information with
 * different severity levels, and it should NOT handle the heavy duty output,
 * since it is less efficient than the standard output methods.
 */


// Log level definition
//-----------------------------------------------------------------------------
enum class LogLevel: int {
    Debug   = 0,
    Info    = 1,
    Step    = 2,
    Note    = 3,
    Warning = 4,
    Error   = 5,
    Fatal   = 6
};

// Converts log level to string
constexpr const char* literal(LogLevel lv) {
    switch(lv) {

        case LogLevel::Debug:   return "Debug";
        case LogLevel::Info:    return "Info";
        case LogLevel::Step:    return "Step";
        case LogLevel::Note:    return "Note";
        case LogLevel::Warning: return "Warning";
        case LogLevel::Error:   return "Error";
        case LogLevel::Fatal:   return "Fatal";
        default:                return "Unknown";
    }
}


class LoggerLevelFlag {
    int _flag = 0;

    /// Helper function to convert loglevel to flag bit.
    static constexpr int _convert(LogLevel lv) { return 1 << static_cast<int>(lv); }
public:
    int value()const { return _flag; }

    /// Flag manipulations
    bool isOnWith(LogLevel lv)const { return _flag & _convert(lv); }

    void turnOn(LogLevel lv) { _flag |= _convert(lv); }
    void turnOff(LogLevel lv) { _flag &= ~_convert(lv); }
    void turnOnAtLeast(LogLevel lv) {
        switch(lv) {
            case LogLevel::Debug:   turnOn(LogLevel::Debug);   [[fallthrough]];
            case LogLevel::Info:    turnOn(LogLevel::Info);    [[fallthrough]];
            case LogLevel::Step:    turnOn(LogLevel::Step);    [[fallthrough]];
            case LogLevel::Note:    turnOn(LogLevel::Note);    [[fallthrough]];
            case LogLevel::Warning: turnOn(LogLevel::Warning); [[fallthrough]];
            case LogLevel::Error:   turnOn(LogLevel::Error);   [[fallthrough]];
            case LogLevel::Fatal:   turnOn(LogLevel::Fatal);
        }
    }
};
struct LoggerOstreamContainer {
    std::ostream* os;

    /// For ofstream only
    bool isOfstream;
    std::string filepath;

    /// Display flags
    LoggerLevelFlag disp; ///< Main switch
    LoggerLevelFlag dispColor;
    LoggerLevelFlag dispTime;
    LoggerLevelFlag dispFile;
    LoggerLevelFlag dispLine;
    LoggerLevelFlag dispFunc;
    LoggerLevelFlag dispLevel;

    /// Other settings
    LoggerLevelFlag flushLevel; ///< Whether to flush after each log.

    LoggerOstreamContainer(std::ostream* os, bool isOfstream): os(os), isOfstream(isOfstream) {}
};

struct LoggerSettings {
    bool supressColorIfRedirected = true;

    std::string delimiterBefore = "[";
    std::string delimiterAfter = "]";
    bool spaceAfterDelimiter = true;

    bool newLineAfterLog = true;
};

/// Stores configurations of the logger
class Logger {
public:
    /// Logger settings
    LoggerSettings settings;

    /// Getters and setters
    const std::vector<LoggerOstreamContainer>& getOsContainers()const { return _osContainers; }

    /// New ostreams
    LoggerOstreamContainer& attachOstream(std::ostream* os, bool isOfstream) {
        _osContainers.emplace_back(os, isOfstream);
        return _osContainers.back();
    }
    LoggerOstreamContainer& addOfstream(const std::string& filepath) {
        _osManaged.emplace_back(new std::ofstream(filepath));
        _osContainers.emplace_back(_osManaged.back().get(), true);
        _osContainers.back().filepath = filepath;
        return _osContainers.back();
    }
    void removeOstream(std::ostream* os) {
        auto itContainer = std::find_if(
            _osContainers.begin(), _osContainers.end(),
            [&](LoggerOstreamContainer& loc) { return loc.os == os; }
        );
        if(itContainer != _osContainers.end()) _osContainers.erase(itContainer);

        auto itManaged = std::find_if(
            _osManaged.begin(), _osManaged.end(),
            [&](std::unique_ptr<std::ostream>& p) { return p.get() == os; }
        );
        if(itManaged != _osManaged.end()) _osManaged.erase(itManaged);
    }

    // Default settings
    void defaultInitialization() {

        auto& scrn = attachOstream(&std::cout, false);
        scrn.disp.turnOnAtLeast(LogLevel::Debug);
        if(!(ioEnv().stdoutRedirected && settings.supressColorIfRedirected)) scrn.dispColor.turnOnAtLeast(LogLevel::Debug);
        scrn.dispTime.turnOnAtLeast(LogLevel::Note);
        scrn.dispFile.turnOnAtLeast(LogLevel::Error);
        scrn.dispLine.turnOnAtLeast(LogLevel::Error);
        scrn.dispFunc.turnOnAtLeast(LogLevel::Error);
        scrn.dispLevel.turnOnAtLeast(LogLevel::Debug);
        scrn.flushLevel.turnOnAtLeast(LogLevel::Debug);
    }

private:
    /// The actual stringstream
    std::ostringstream _oss;

    /// Ostream containers
    std::vector<LoggerOstreamContainer> _osContainers;

    /// Managed ostreams. will be destroyed when instance of this class is going out of scope
    std::vector<std::unique_ptr<std::ostream>> _osManaged;
};

// Default logger
inline Logger& defaultLogger() {
    struct DefaultLogger {
        Logger l;
        DefaultLogger() {
            l.defaultInitialization();
        }
    };
    static DefaultLogger dl;
    return dl.l;
}

namespace internal {
/// This is the class that prints the log when destructed.
class LogWriter {
public:
    /// Constructor accepts environmental information.
    LogWriter(const char* curFile, int curLine, const char* curFunc, LogLevel lv, const Logger& logger)
        : _logger(logger), _lv(lv), _curFile(curFile), _curLine(curLine), _curFunc(curFunc)
    {
        oss_.precision(17);
    }

    /// Destructor dispatches the log.
    ~LogWriter() { logDispatch(); }

    /// Copy is not allowed
    LogWriter(const LogWriter&) = delete;
    LogWriter& operator=(const LogWriter&) = delete;

    /// Log generation and dispatch
    void logDispatch();

    /// Overloads operator<< for normal type and user defined class type
    template<typename MsgType>
    LogWriter& operator<<(MsgType&& msg) {
        oss_ << std::forward<MsgType>(msg);
        return *this;
    }
    LogWriter& operator<<(std::ostream& (*manip)(std::ostream&)) {
        oss_ << manip;
        return *this;
    }

private:
    /// The ref to the logger information
    const Logger& _logger;

    /// The content of the log
    std::ostringstream oss_;
    /// Level of current log
    const LogLevel _lv;
    /// Other information of the log
    const char* _curFile;
    const int _curLine;
    const char* _curFunc;
};
} // namespace internal

} // namespace logger
} // namespace medyan

/// Exposed macro
#undef ERROR // Fuck <windows.h>

#if defined(COMPILER_CLANG) || defined(COMPILER_GCC)
    #define MEDYAN_LOG_FUNC __PRETTY_FUNCTION__
#elif defined(COMPILER_MSVC)
    #define MEDYAN_LOG_FUNC __FUNCSIG__
#else
    #define MEDYAN_LOG_FUNC __func__
#endif
#define MEDYAN_WRITE_LOG(whichLogger, logLevel) ::medyan::logger::internal::LogWriter(__FILE__, __LINE__, MEDYAN_LOG_FUNC, logLevel, whichLogger)

#define MEDYAN_LOG_GEN_DEBUG(whichLogger)   MEDYAN_WRITE_LOG(whichLogger, ::medyan::logger::LogLevel::Debug)
#define MEDYAN_LOG_GEN_INFO(whichLogger)    MEDYAN_WRITE_LOG(whichLogger, ::medyan::logger::LogLevel::Info)
#define MEDYAN_LOG_GEN_STEP(whichLogger)    MEDYAN_WRITE_LOG(whichLogger, ::medyan::logger::LogLevel::Step)
#define MEDYAN_LOG_GEN_NOTE(whichLogger)    MEDYAN_WRITE_LOG(whichLogger, ::medyan::logger::LogLevel::Note)
#define MEDYAN_LOG_GEN_WARNING(whichLogger) MEDYAN_WRITE_LOG(whichLogger, ::medyan::logger::LogLevel::Warning)
#define MEDYAN_LOG_GEN_ERROR(whichLogger)   MEDYAN_WRITE_LOG(whichLogger, ::medyan::logger::LogLevel::Error)
#define MEDYAN_LOG_GEN_FATAL(whichLogger)   MEDYAN_WRITE_LOG(whichLogger, ::medyan::logger::LogLevel::Fatal)

#define MEDYAN_LOG_GEN(logLevel) MEDYAN_LOG_GEN_##logLevel(::medyan::logger::defaultLogger())

/// User interface
#define LOG(logLevel) MEDYAN_LOG_GEN(logLevel)

#endif

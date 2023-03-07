#ifndef MEDYAN_Util_Environment_hpp
#define MEDYAN_Util_Environment_hpp

// Target platforms
//-----------------------------------------------------------------------------
#if defined(__APPLE__) && defined(__MACH__)
    #define PLATFORM_MACOS
#endif
#if defined(__linux__)
    #define PLATFORM_LINUX
#endif
#if defined(__unix__) || defined(PLATFORM_LINUX) || defined(PLATFORM_MACOS)
    #define PLATFORM_UNIX_LIKE
#endif
#if defined(_WIN32)
    #define PLATFORM_WINDOWS
#endif

// Compiler
//-----------------------------------------------------------------------------
#if defined(__clang__)
    #define COMPILER_CLANG
#elif defined(__GNUC__)
    #define COMPILER_GCC
#elif defined(_MSC_VER)
    #define COMPILER_MSVC
#endif

// IO environments
//-----------------------------------------------------------------------------

// Struct for getting and initializing io environment
struct IoEnv {
    bool stdoutRedirected = false;
    bool stderrRedirected = false;
    IoEnv();
};

inline const IoEnv& ioEnv() {
    static IoEnv ie;
    return ie;
}


namespace medyan {

// Floating point environments
//-----------------------------------------------------------------------------
void trapInvalidFP(bool trap = true);


// GUI built
//-----------------------------------------------------------------------------
#ifdef NO_GUI
    constexpr bool builtGui = false;
#else
    constexpr bool builtGui = true;
#endif


} // namespace medyan

#endif

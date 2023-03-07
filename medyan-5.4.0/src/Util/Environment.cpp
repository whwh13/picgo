#include "Environment.hpp"

#include <iostream>

#ifdef PLATFORM_UNIX_LIKE
    #include <fenv.h>
    #include <unistd.h>
#elif defined(PLATFORM_WINDOWS)
    #include <float.h>
    #include <Windows.h>
#endif

IoEnv::IoEnv() {

#ifdef PLATFORM_UNIX_LIKE
    // Detect redirection
    stdoutRedirected = !isatty(STDOUT_FILENO);
    stderrRedirected = !isatty(STDERR_FILENO);

#elif defined(PLATFORM_WINDOWS)
    // Detect redirection and set VT mode
    HANDLE hOut = GetStdHandle(STD_OUTPUT_HANDLE);
    if(hOut == INVALID_HANDLE_VALUE) {
        stdoutRedirected = true;
    } else {
        DWORD mode = 0;
        stdoutRedirected = !GetConsoleMode(hOut, &mode);
        if(!stdoutRedirected) {
            mode |= ENABLE_VIRTUAL_TERMINAL_PROCESSING;
            SetConsoleMode(hOut, mode);
        }
    }
    HANDLE hErr = GetStdHandle(STD_ERROR_HANDLE);
    if(hErr == INVALID_HANDLE_VALUE) {
        stderrRedirected = true;
    } else {
        DWORD mode = 0;
        stderrRedirected = !GetConsoleMode(hErr, &mode);
        if(!stderrRedirected) {
            mode |= ENABLE_VIRTUAL_TERMINAL_PROCESSING;
            SetConsoleMode(hErr, mode);
        }
    }

#endif

} // IoEnv::IoEnv()

namespace medyan {

// Floating point environments
//-----------------------------------------------------------------------------
void trapInvalidFP(bool trap) {
#ifdef PLATFORM_LINUX
    if(trap) {
        feenableexcept(FE_INVALID);
    } else {
        fedisableexcept(FE_INVALID);
    }
#elif defined(PLATFORM_WINDOWS)
    unsigned int currentFpControl;
    _controlfp_s(&currentFpControl, 0, 0);
    if(trap) {
        _controlfp_s(&currentFpControl, currentFpControl & ~_EM_INVALID, _MCW_EM);
    } else {
        _controlfp_s(&currentFpControl, currentFpControl | _EM_INVALID, _MCW_EM);
    }
#else
    std::cout << "Invalid floating point operation trapping is not supported on this platform yet." << std::endl;
#endif
}


} // namespace medyan

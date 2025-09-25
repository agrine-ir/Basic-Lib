#pragma once
#include <string>

#ifdef AGRINECMD_EXPORTS
#define AGRINECMD_API __declspec(dllexport)
#else
#define AGRINECMD_API __declspec(dllimport)
#endif

extern "C" {
    AGRINECMD_API const char* RunCmd(const char* command);
    AGRINECMD_API const char* RunPowerShell(const char* command);
}

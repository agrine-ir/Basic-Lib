#include "pch.h"
#include "AgrineCmdAPI.h"
#include "CmdExecutor.h"
#include "PowerShellExecutor.h"

// Static buffer to hold output
static std::string resultBuffer;

const char* RunCmd(const char* command) {
    resultBuffer = CmdExecutor::RunCommand(command);
    return resultBuffer.c_str();
}

const char* RunPowerShell(const char* command) {
    resultBuffer = PowerShellExecutor::RunCommand(command);
    return resultBuffer.c_str();
}

#include "pch.h"
#include "PowerShellExecutor.h"
#include "CmdExecutor.h"

std::string PowerShellExecutor::RunCommand(const std::string& command) {
    std::string psCommand = "powershell -Command \"" + command + "\"";
    return CmdExecutor::RunCommand(psCommand);
}

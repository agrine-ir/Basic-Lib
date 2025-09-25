#include "pch.h"
#include "PowerShellExecutor.h"
#include "CmdExecutor.h"

/**
 * @brief Executes a PowerShell command by wrapping it for CMD execution
 * @param command Command string
 * @return Output from PowerShell execution
 */
std::string PowerShellExecutor::Execute(const std::string& command) {
    CmdExecutor executor;
    std::string psCommand = "powershell -Command \"" + command + "\"";
    return executor.Execute(psCommand);
}

#include "pch.h"
#include "AgrineCmdAPI.h"
#include "CmdExecutor.h"
#include "PowerShellExecutor.h"

#include <mutex>
#include <unordered_map>

/**
 * @brief Thread-safe storage for result strings
 */
static std::mutex resultMutex;
static std::unordered_map<std::thread::id, std::string> resultStorage;

const char* RunCmd(const char* command) {
    CmdExecutor executor;
    std::string result = executor.Execute(command);

    std::lock_guard<std::mutex> lock(resultMutex);
    resultStorage[std::this_thread::get_id()] = result;

    return resultStorage[std::this_thread::get_id()].c_str();
}

const char* RunPowerShell(const char* command) {
    PowerShellExecutor executor;
    std::string result = executor.Execute(command);

    std::lock_guard<std::mutex> lock(resultMutex);
    resultStorage[std::this_thread::get_id()] = result;

    return resultStorage[std::this_thread::get_id()].c_str();
}

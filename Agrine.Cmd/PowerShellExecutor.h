#pragma once
#include <string>

class PowerShellExecutor {
public:
    static std::string RunCommand(const std::string& command);
};

#pragma once
#include <string>

class CmdExecutor {
public:
    static std::string RunCommand(const std::string& command);
};

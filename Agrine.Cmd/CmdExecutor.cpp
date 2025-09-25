#include "pch.h"
#include "CmdExecutor.h"

#include <cstdio>

std::string CmdExecutor::RunCommand(const std::string& command) {
    char buffer[256];
    std::string result;

    FILE* pipe = _popen(command.c_str(), "r");
    if (!pipe) return "Error: cannot open pipe";

    while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
        result += buffer;
    }

    _pclose(pipe);
    return result;
}

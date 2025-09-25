#pragma once
#include "ICommandExecutor.h"

/**
 * @brief Executes commands using the standard Windows Command Prompt (cmd.exe)
 */
class CmdExecutor : public ICommandExecutor {
public:
    std::string Execute(const std::string& command) override;
};

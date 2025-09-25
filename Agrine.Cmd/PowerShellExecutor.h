#pragma once
#include "ICommandExecutor.h"

/**
 * @brief Executes commands using Windows PowerShell
 */
class PowerShellExecutor : public ICommandExecutor {
public:
    std::string Execute(const std::string& command) override;
};

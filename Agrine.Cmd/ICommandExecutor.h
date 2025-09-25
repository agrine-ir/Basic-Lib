#pragma once
#include <string>

/**
 * @brief Interface for command executors
 */
class ICommandExecutor {
public:
    virtual ~ICommandExecutor() = default;

    /**
     * @brief Executes a given command and returns the result.
     * @param command Command string to execute
     * @return Output from the command execution
     */
    virtual std::string Execute(const std::string& command) = 0;
};

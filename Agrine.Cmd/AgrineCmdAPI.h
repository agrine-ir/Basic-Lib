#pragma once
#include <string>

#ifdef AGRINECMD_EXPORTS
#define AGRINECMD_API __declspec(dllexport)
#else
#define AGRINECMD_API __declspec(dllimport)
#endif

extern "C" {

	/**
	 * @brief Runs a CMD command.
	 * @param command Command string
	 * @return Pointer to result string (needs copying in managed code)
	 */
	AGRINECMD_API const char* RunCmd(const char* command);

	/**
	 * @brief Runs a PowerShell command.
	 * @param command Command string
	 * @return Pointer to result string (needs copying in managed code)
	 */
	AGRINECMD_API const char* RunPowerShell(const char* command);
}

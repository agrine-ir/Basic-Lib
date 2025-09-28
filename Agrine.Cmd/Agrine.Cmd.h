#pragma once

// Main header file for the Agrine.Cmd library
// This header exposes the public API of the DLL.

using namespace System;

namespace Agrine {
	namespace Cmd {

		public ref class CmdExecutor
		{
		public:
			// Executes a command using CMD.exe and returns the output.
			static String^ Execute(String^ command);
		};

		public ref class PowerShellExecutor
		{
		public:
			// Executes a command using PowerShell and returns the output.
			static String^ Execute(String^ command);
		};
	}
}

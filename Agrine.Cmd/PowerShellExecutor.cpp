#include "pch.h"
#include "PowerShellExecutor.h"

using namespace System;
using namespace System::Diagnostics;

namespace Agrine {
    namespace Cmd {

        String^ PowerShellExecutor::Execute(String^ command)
        {
            ProcessStartInfo^ psi = gcnew ProcessStartInfo();
            psi->FileName = "powershell.exe";
            psi->Arguments = "-NoProfile -Command " + command;
            psi->RedirectStandardOutput = true;
            psi->RedirectStandardError = true;
            psi->UseShellExecute = false;
            psi->CreateNoWindow = true;

            Process^ process = gcnew Process();
            process->StartInfo = psi;

            process->Start();

            String^ output = process->StandardOutput->ReadToEnd();
            String^ error = process->StandardError->ReadToEnd();

            process->WaitForExit();

            if (!String::IsNullOrEmpty(error))
                return output + "\nERROR:\n" + error;

            return output;
        }
    }
}

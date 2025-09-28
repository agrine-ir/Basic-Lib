#include "pch.h"
#include "CmdExecutor.h"

using namespace System;
using namespace System::Diagnostics;
using namespace System::Text;

namespace Agrine {
    namespace Cmd {

        String^ CmdExecutor::Execute(String^ command)
        {
            ProcessStartInfo^ psi = gcnew ProcessStartInfo();
            psi->FileName = "cmd.exe";
            psi->Arguments = "/c " + command; // "/c" tells cmd to execute and exit
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

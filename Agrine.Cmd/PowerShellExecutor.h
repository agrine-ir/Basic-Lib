#pragma once

using namespace System;

namespace Agrine {
    namespace Cmd {

        public ref class PowerShellExecutor
        {
        public:
            static String^ Execute(String^ command);
        };
    }
}


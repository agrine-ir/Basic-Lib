#pragma once

using namespace System;

namespace Agrine {
    namespace Cmd {

        public ref class CmdExecutor
        {
        public:
            static String^ Execute(String^ command);
        };

    }
}

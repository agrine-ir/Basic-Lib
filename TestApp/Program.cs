using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using Agrine.Cmd.Wrapper;
using Agrine.Cmd.Wrapper.Core;

namespace TestApp
{
    internal class Program
    {

        static void Main(string[] args)
        {

            try
            {
                Console.WriteLine("Testing Agrine.Cmd.Wrapper...");

                var cmdResult = AgrineCmdFacade.RunCmd("dir");
                Console.WriteLine("CMD Output:\n" + cmdResult);

                var psResult = AgrineCmdFacade.RunPowerShell("Get-Process");
                Console.WriteLine("PowerShell Output:\n" + psResult);
            }
            catch (Exception ex)
            {
                Console.WriteLine("Error: " + ex.Message);
            }

            Console.WriteLine("Press any key to exit...");
            Console.ReadKey();
        }
    }
}


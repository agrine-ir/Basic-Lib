using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Agrine.Cmd;

namespace TestPro
{
    internal class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("Microsoft");
            string result = Agrine.Cmd.PowerShellExecutor.Execute("Get-Location");
            Console.WriteLine(result);
            Console.ReadKey();
        }
    }
}

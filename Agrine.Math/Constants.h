#pragma once

// Mathematical constants for Agrine.Math
// Project: Agrine.Math (C++/CLI, .NET Framework 4.6.2)

using namespace System;

namespace Agrine {
    namespace Math {
        namespace Core {

            /// <summary>
            /// Provides common mathematical constants.
            /// These are declared as literal constants for direct use.
            /// </summary>
            public ref class Constants abstract sealed
            {
            public:
                // Fundamental constants
                literal double Pi = 3.14159265358979323846;  // π
                literal double TwoPi = 6.28318530717958647692;  // 2π
                literal double HalfPi = 1.57079632679489661923;  // π/2
                literal double E = 2.71828182845904523536;  // Euler's number
                literal double Sqrt2 = 1.41421356237309504880;  // √2
                literal double Sqrt3 = 1.73205080756887729353;  // √3
                literal double GoldenR = 1.61803398874989484820;  // Golden ratio φ

                // Precision-related
                literal double Epsilon = 1e-12;   // Small threshold for floating-point comparisons
                literal double Tolerance = 1e-9;    // Default numerical tolerance

                // Conversion factors
                literal double Deg2Rad = Pi / 180.0;  // Degrees → Radians
                literal double Rad2Deg = 180.0 / Pi;  // Radians → Degrees
            };
        }
    }
}

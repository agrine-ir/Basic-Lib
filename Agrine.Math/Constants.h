#pragma once

// Core/Constants.h
// Agrine.Math (C++/CLI, .NET Framework 4.6.2)

using namespace System;

namespace Agrine {
    namespace Math {
        namespace Core {

            /// <summary>
            /// Provides common mathematical constants.
            /// Defined explicitly to avoid initializer ordering / compile-time issues.
            /// </summary>
            public ref class Constants abstract sealed
            {
            public:
                // Fundamental constants (explicit numeric literals)
                literal double Pi = 3.14159265358979323846;      // π
                literal double TwoPi = 6.28318530717958647692;      // 2π
                literal double HalfPi = 1.57079632679489661923;      // π/2
                literal double E = 2.71828182845904523536;      // Euler's number
                literal double Sqrt2 = 1.41421356237309504880;      // √2
                literal double Sqrt3 = 1.73205080756887729353;      // √3
                literal double GoldenR = 1.61803398874989484820;      // φ

                // Precision-related
                literal double Epsilon = 1e-12;                      // Small threshold for comparisons
                literal double Tolerance = 1e-9;                       // Default numerical tolerance

                // Conversion factors (explicit)
                literal double Deg2Rad = 0.01745329251994329576923690768489; // Pi / 180
                literal double Rad2Deg = 57.295779513082320876798154814105;  // 180 / Pi
            };
        }
    }
}

#pragma once

// ============================================================
// Constants.h - Mathematical Constants for Agrine.Math Library
// ============================================================
//
// This header defines commonly used mathematical constants
// with high precision for numerical calculations.
//
// These constants are available to all modules in the library.
//
// ============================================================

namespace Agrine
{
	namespace Math
	{
		namespace Constants
		{
			// Mathematical constants
			constexpr double PI = 3.14159265358979323846; // π: Circle constant
			constexpr double TWO_PI = 6.28318530717958647692; // 2π
			constexpr double HALF_PI = 1.57079632679489661923; // π/2
			constexpr double E = 2.71828182845904523536; // Euler's number
			constexpr double LN2 = 0.69314718055994530942; // Natural log of 2
			constexpr double LN10 = 2.30258509299404568402; // Natural log of 10
			constexpr double SQRT2 = 1.41421356237309504880; // √2
			constexpr double SQRT1_2 = 0.70710678118654752440; // 1/√2
			constexpr double GAMMA = 0.57721566490153286060; // Euler–Mascheroni constant

			// Conversion constants
			constexpr double DEG_TO_RAD = PI / 180.0; // Degrees to radians
			constexpr double RAD_TO_DEG = 180.0 / PI; // Radians to degrees

			// Infinity and epsilon
			constexpr double INF = std::numeric_limits<double>::infinity(); // Positive infinity
			constexpr double NAN_VAL = std::numeric_limits<double>::quiet_NaN(); // Not-a-Number
			constexpr double EPSILON = std::numeric_limits<double>::epsilon();  // Machine epsilon
		}

	}
}

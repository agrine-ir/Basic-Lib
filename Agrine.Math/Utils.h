#pragma once

// Utility functions for Agrine.Math
// Project: Agrine.Math (C++/CLI, .NET Framework 4.6.2)

using namespace System;

#include "Agrine.Math.h"

namespace Agrine {
    namespace Math {
        namespace Core {

            /// <summary>
            /// Provides common mathematical utility functions.
            /// </summary>
            public ref class Utils abstract sealed
            {
            public:
                /// <summary>
                /// Checks if a double is approximately zero within tolerance.
                /// </summary>
                static bool IsZero(double value, double eps = Constants::Epsilon) {
                    return System::Math::Abs(value) <= eps;
                }

                /// <summary>
                /// Compares two doubles for approximate equality.
                /// </summary>
                static bool AreEqual(double a, double b, double tol = Constants::Tolerance) {
                    return System::Math::Abs(a - b) <= tol;
                }

                /// <summary>
                /// Clamps a value between min and max.
                /// </summary>
                static double Clamp(double value, double min, double max) {
                    if (value < min) return min;
                    if (value > max) return max;
                    return value;
                }

                /// <summary>
                /// Swaps two values by reference.
                /// </summary>
                generic <typename T>
                static void Swap(T% a, T% b) {
                    T temp = a;
                    a = b;
                    b = temp;
                }

                /// <summary>
                /// Returns the factorial of n (n!).
                /// Uses iterative approach to avoid recursion overhead.
                /// </summary>
                static UInt64 Factorial(int n) {
                    if (n < 0)
                        throw gcnew ArgumentException("Factorial is not defined for negative numbers.");

                    UInt64 result = 1;
                    for (int i = 2; i <= n; i++) {
                        result *= (UInt64)i;
                    }
                    return result;
                }

                /// <summary>
                /// Computes the greatest common divisor (GCD) using Euclidean algorithm.
                /// </summary>
                static int GCD(int a, int b) {
                    a = System::Math::Abs(a);
                    b = System::Math::Abs(b);
                    while (b != 0) {
                        int temp = b;
                        b = a % b;
                        a = temp;
                    }
                    return a;
                }

                /// <summary>
                /// Computes the least common multiple (LCM).
                /// </summary>
                static int LCM(int a, int b) {
                    if (a == 0 || b == 0) return 0;
                    return (System::Math::Abs(a) / GCD(a, b)) * System::Math::Abs(b);
                }

                /// <summary>
                /// Normalizes an angle (in radians) to [0, 2π).
                /// </summary>
                static double NormalizeAngle(double rad) {
                    double result = fmod(rad, Constants::TwoPi);
                    if (result < 0) result += Constants::TwoPi;
                    return result;
                }

                /// <summary>
                /// Converts degrees to radians.
                /// </summary>
                static double ToRadians(double deg) {
                    return deg * Constants::Deg2Rad;
                }

                /// <summary>
                /// Converts radians to degrees.
                /// </summary>
                static double ToDegrees(double rad) {
                    return rad * Constants::Rad2Deg;
                }
            };
        }
    }
}

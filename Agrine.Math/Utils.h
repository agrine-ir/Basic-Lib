#pragma once
#include "MathCore.h"
#include <vector>
#include <utility>

namespace Agrine
{
    namespace Math
    {
        namespace Utils
        {
            // ============================================================
            // Utils.h - Utility Functions for Agrine.Math
            // ============================================================
            // This module provides helper and utility functions that can be
            // reused across other modules (Algebra, Geometry, Statistics, ...).
            // ============================================================

            // ----------------------------
            // Basic Math Helpers
            // ----------------------------
            AGRINE_API int GCD(int a, int b);
            AGRINE_API int LCM(int a, int b);
            AGRINE_API long long Factorial(int n);
            AGRINE_API double Power(double base, int exponent);

            // ----------------------------
            // Number Utilities
            // ----------------------------
            AGRINE_API bool IsPrime(int n);
            AGRINE_API std::vector<int> PrimeFactors(int n);
            AGRINE_API double Clamp(double value, double minVal, double maxVal);
            AGRINE_API bool AlmostEqual(double a, double b, double epsilon = 1e-9);

            // ----------------------------
            // Random Utilities
            // ----------------------------
            AGRINE_API double RandomDouble(double minVal = 0.0, double maxVal = 1.0);
            AGRINE_API int RandomInt(int minVal, int maxVal);
            AGRINE_API double RandomChoice(const std::vector<double>& data);

            // ----------------------------
            // Vector Utilities
            // ----------------------------
            AGRINE_API std::vector<double> Normalize(const std::vector<double>& data);
            AGRINE_API double DotProduct(const std::vector<double>& a, const std::vector<double>& b);
            AGRINE_API std::pair<double, double> MinMax(const std::vector<double>& data);

            // ----------------------------
            // Matrix Utilities
            // ----------------------------
            AGRINE_API std::vector<std::vector<double>> IdentityMatrix(size_t n);
            AGRINE_API std::vector<std::vector<double>> ZeroMatrix(size_t rows, size_t cols);
            AGRINE_API std::vector<std::vector<double>> RandomMatrix(size_t rows, size_t cols,
                double minVal = 0.0, double maxVal = 1.0);
            AGRINE_API std::vector<std::vector<double>> Transpose(const std::vector<std::vector<double>>& mat);
            AGRINE_API bool SameDimensions(const std::vector<std::vector<double>>& a,
                const std::vector<std::vector<double>>& b);
        }
    }
}

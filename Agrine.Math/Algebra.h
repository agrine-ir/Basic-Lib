#pragma once
#include "MathCore.h"
#include <cstdint> // For integer types like int64_t

// ============================================================
// Algebra Module - Basic Arithmetic and Elementary Operations
// ============================================================
// This module provides fundamental arithmetic operations
// such as addition, subtraction, multiplication, division,
// exponentiation, and factorial.
// All functions are exported for external use in C# and other projects.
// ============================================================

extern "C"
{
    // Addition of two real numbers
    AGRINE_API double Add(double a, double b);

    // Subtraction of two real numbers
    AGRINE_API double Subtract(double a, double b);

    // Multiplication of two real numbers
    AGRINE_API double Multiply(double a, double b);

    // Division of two real numbers (returns 0 if division by zero occurs)
    AGRINE_API double Divide(double a, double b);

    // Exponentiation (base^exponent)
    AGRINE_API double Power(double base, double exponent);

    // Factorial of a non-negative integer (n!)
    AGRINE_API int64_t Factorial(int n);
}

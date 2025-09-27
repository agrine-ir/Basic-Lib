#pragma once
#include "MathCore.h"

/**
 * ============================================================
 * Algebra Module - Agrine.Math
 * ============================================================
 * This module provides both basic and advanced algebraic operations:
 * - Basic arithmetic
 * - Power and factorial functions
 * - Logarithmic and exponential functions
 * - Numerical derivative and limit calculations
 * - Solving linear equations (simple cases)
 * - Polynomial evaluation
 * ============================================================
 */

extern "C"
{
    // Basic Arithmetic
    AGRINE_API double Add(double a, double b);
    AGRINE_API double Subtract(double a, double b);
    AGRINE_API double Multiply(double a, double b);
    AGRINE_API double Divide(double a, double b);

    // Advanced Arithmetic
    AGRINE_API double Power(double base, double exponent);
    AGRINE_API int64_t Factorial(int n);
    AGRINE_API double Logarithm(double value, double base);
    AGRINE_API double Exponential(double exponent);

    // Numerical Calculus
    AGRINE_API double Derivative(double (*func)(double), double x, double h);
    AGRINE_API double Limit(double (*func)(double), double x, double h);

    // Polynomial Functions
    AGRINE_API double EvaluatePolynomial(const double* coefficients, size_t degree, double x);

    // Linear Equation Solver (simple Ax = b for 2x2 and 3x3 cases)
    AGRINE_API bool SolveLinear2x2(const double* A, const double* b, double* x);
    AGRINE_API bool SolveLinear3x3(const double* A, const double* b, double* x);
}

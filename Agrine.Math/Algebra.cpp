#include "pch.h"
#include "Algebra.h"
#include <cmath>     // For pow()
#include <stdexcept> // For exceptions

// ============================================================
// Implementation of Algebra Module
// ============================================================

double Add(double a, double b)
{
    return a + b;
}

double Subtract(double a, double b)
{
    return a - b;
}

double Multiply(double a, double b)
{
    return a * b;
}

double Divide(double a, double b)
{
    if (b == 0.0)
    {
        // Division by zero handling: return 0 (could also throw an error)
        return 0.0;
    }
    return a / b;
}

double Power(double base, double exponent)
{
    return std::pow(base, exponent);
}

int64_t Factorial(int n)
{
    if (n < 0)
    {
        // Negative factorial is undefined
        return -1;
    }

    int64_t result = 1;
    for (int i = 2; i <= n; ++i)
    {
        result *= i;
    }
    return result;
}

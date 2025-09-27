#include "pch.h"
#include "Algebra.h"
#include <cmath>
#include <stdexcept>

/** ---------------- Basic Arithmetic ---------------- */

double Add(double a, double b) { return a + b; }
double Subtract(double a, double b) { return a - b; }
double Multiply(double a, double b) { return a * b; }
double Divide(double a, double b) { return (b != 0.0) ? (a / b) : 0.0; }

/** ---------------- Advanced Arithmetic ---------------- */

double Power(double base, double exponent) { return std::pow(base, exponent); }

int64_t Factorial(int n)
{
    if (n < 0) return -1; // Undefined for negative integers
    int64_t result = 1;
    for (int i = 2; i <= n; ++i) result *= i;
    return result;
}

double Logarithm(double value, double base)
{
    if (value <= 0.0 || base <= 0.0) return NAN;
    return std::log(value) / std::log(base);
}

double Exponential(double exponent) { return std::exp(exponent); }

/** ---------------- Numerical Calculus ---------------- */

// Numerical derivative using central difference method
double Derivative(double (*func)(double), double x, double h)
{
    return (func(x + h) - func(x - h)) / (2 * h);
}

// Numerical limit estimation
double Limit(double (*func)(double), double x, double h)
{
    return func(x + h);
}

/** ---------------- Polynomial Functions ---------------- */

// Evaluate polynomial given coefficients array
double EvaluatePolynomial(const double* coefficients, size_t degree, double x)
{
    double result = 0.0;
    for (size_t i = 0; i <= degree; ++i)
    {
        result += coefficients[i] * std::pow(x, static_cast<double>(i));
    }
    return result;
}

/** ---------------- Linear Equation Solver ---------------- */

// Solve 2x2 system: A * x = b
bool SolveLinear2x2(const double* A, const double* b, double* x)
{
    double det = A[0] * A[3] - A[1] * A[2];
    if (det == 0.0) return false;
    x[0] = (b[0] * A[3] - b[1] * A[1]) / det;
    x[1] = (A[0] * b[1] - A[2] * b[0]) / det;
    return true;
}

// Solve 3x3 system: A * x = b (Cramer's rule)
bool SolveLinear3x3(const double* A, const double* b, double* x)
{
    double det = A[0] * (A[4] * A[8] - A[5] * A[7]) - A[1] * (A[3] * A[8] - A[5] * A[6]) + A[2] * (A[3] * A[7] - A[4] * A[6]);
    if (det == 0.0) return false;

    x[0] = (b[0] * (A[4] * A[8] - A[5] * A[7]) - A[1] * (b[1] * A[8] - A[5] * b[2]) + A[2] * (b[1] * A[7] - A[4] * b[2])) / det;
    x[1] = (A[0] * (b[1] * A[8] - A[5] * b[2]) - b[0] * (A[3] * A[8] - A[5] * A[6]) + A[2] * (A[3] * b[2] - b[1] * A[6])) / det;
    x[2] = (A[0] * (A[4] * b[2] - b[1] * A[7]) - A[1] * (A[3] * b[2] - b[1] * A[6]) + b[0] * (A[3] * A[7] - A[4] * A[6])) / det;

    return true;
}

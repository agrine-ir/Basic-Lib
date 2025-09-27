#include "pch.h"
#include "Complex.h"
#include <cmath>
#include <stdexcept>

/** ---------------- Complex Creation ---------------- */
ComplexNumber CreateComplex(double real, double imag)
{
    return ComplexNumber(real, imag);
}

/** ---------------- Basic Arithmetic ---------------- */
ComplexNumber AddComplex(ComplexNumber a, ComplexNumber b)
{
    return ComplexNumber(a.real + b.real, a.imag + b.imag);
}

ComplexNumber SubtractComplex(ComplexNumber a, ComplexNumber b)
{
    return ComplexNumber(a.real - b.real, a.imag - b.imag);
}

ComplexNumber MultiplyComplex(ComplexNumber a, ComplexNumber b)
{
    return ComplexNumber(
        a.real * b.real - a.imag * b.imag,
        a.real * b.imag + a.imag * b.real
    );
}

ComplexNumber DivideComplex(ComplexNumber a, ComplexNumber b)
{
    double denominator = b.real * b.real + b.imag * b.imag;
    if (denominator == 0.0) throw std::runtime_error("Division by zero in complex division");
    return ComplexNumber(
        (a.real * b.real + a.imag * b.imag) / denominator,
        (a.imag * b.real - a.real * b.imag) / denominator
    );
}

/** ---------------- Advanced Operations ---------------- */
double ComplexMagnitude(ComplexNumber z)
{
    return std::sqrt(z.real * z.real + z.imag * z.imag);
}

double ComplexPhase(ComplexNumber z)
{
    return std::atan2(z.imag, z.real); // returns phase in radians
}

ComplexNumber ComplexConjugate(ComplexNumber z)
{
    return ComplexNumber(z.real, -z.imag);
}

ComplexNumber ComplexPower(ComplexNumber z, double exponent)
{
    double r = ComplexMagnitude(z);
    double theta = ComplexPhase(z);
    double r_pow = std::pow(r, exponent);
    double new_theta = theta * exponent;
    return ComplexNumber(r_pow * std::cos(new_theta), r_pow * std::sin(new_theta));
}

ComplexNumber ComplexRoot(ComplexNumber z, int n)
{
    if (n <= 0) throw std::runtime_error("Root degree must be positive");
    double r = ComplexMagnitude(z);
    double theta = ComplexPhase(z);
    double r_root = std::pow(r, 1.0 / n);
    double new_theta = theta / n;
    return ComplexNumber(r_root * std::cos(new_theta), r_root * std::sin(new_theta));
}

/** ---------------- Complex Exponential and Logarithm ---------------- */
ComplexNumber ComplexExp(ComplexNumber z)
{
    double exp_r = std::exp(z.real);
    return ComplexNumber(exp_r * std::cos(z.imag), exp_r * std::sin(z.imag));
}

ComplexNumber ComplexLog(ComplexNumber z)
{
    double r = ComplexMagnitude(z);
    double theta = ComplexPhase(z);
    return ComplexNumber(std::log(r), theta);
}

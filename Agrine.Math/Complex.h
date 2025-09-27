#pragma once
#include "MathCore.h"

/**
 * ============================================================
 * Complex Module - Agrine.Math
 * ============================================================
 * This module provides full support for complex numbers:
 * - Basic arithmetic (addition, subtraction, multiplication, division)
 * - Magnitude and phase calculation
 * - Complex conjugate
 * - Power and root calculations
 * - Complex exponential and logarithm
 * ============================================================
 */

struct ComplexNumber
{
    double real;    // Real part
    double imag;    // Imaginary part

    ComplexNumber(double r = 0.0, double i = 0.0) : real(r), imag(i) {}
};

extern "C"
{
    // Creation and basic operations
    AGRINE_API ComplexNumber CreateComplex(double real, double imag);
    AGRINE_API ComplexNumber AddComplex(ComplexNumber a, ComplexNumber b);
    AGRINE_API ComplexNumber SubtractComplex(ComplexNumber a, ComplexNumber b);
    AGRINE_API ComplexNumber MultiplyComplex(ComplexNumber a, ComplexNumber b);
    AGRINE_API ComplexNumber DivideComplex(ComplexNumber a, ComplexNumber b);

    // Advanced operations
    AGRINE_API double ComplexMagnitude(ComplexNumber z);
    AGRINE_API double ComplexPhase(ComplexNumber z);
    AGRINE_API ComplexNumber ComplexConjugate(ComplexNumber z);
    AGRINE_API ComplexNumber ComplexPower(ComplexNumber z, double exponent);
    AGRINE_API ComplexNumber ComplexRoot(ComplexNumber z, int n);

    // Complex exponential and logarithm
    AGRINE_API ComplexNumber ComplexExp(ComplexNumber z);
    AGRINE_API ComplexNumber ComplexLog(ComplexNumber z);
}

#pragma once

// Complex.h - Advanced complex number utilities (C++/CLI)
// Part of Agrine.Math (C++/CLI, .NET Framework 4.6.2)

#include "Constants.h"
#include "Utils.h"
#include "Exceptions.h"
#include <cmath>

using namespace System;
using namespace Agrine::Math::Core;

namespace Agrine {
    namespace Math {
        namespace Advanced {

            /// <summary>
            /// Complex - value type representing a complex number (double precision).
            /// Designed for numeric computations (fast, value-type semantics).
            /// </summary>
            public value struct Complex
            {
            public:
                // Fields
                double Re;
                double Im;

                // Constructors
                Complex(double real, double imag);
                Complex(double real); // imag = 0

                // Properties
                property double Real { double get(); void set(double); }
                property double Imag { double get(); void set(double); }

                // Basic constants
                static Complex Zero();
                static Complex One();
                static Complex I();

                // Basic operations
                double Abs();               // magnitude
                double AbsSquared();        // squared magnitude
                double Arg();               // principal argument in (-pi, pi]
                Complex Conjugate();
                bool IsZero();              // uses Core::Constants::Tolerance

                // Operators
                static Complex operator+(Complex a, Complex b);
                static Complex operator-(Complex a, Complex b);
                static Complex operator*(Complex a, Complex b);
                static Complex operator/(Complex a, Complex b);
                static Complex operator-(Complex a); // unary neg

                // scalar ops
                static Complex operator*(Complex a, double k);
                static Complex operator*(double k, Complex a);
                static Complex operator/(Complex a, double k);

                // Advanced functions (static)
                static Complex FromPolar(double r, double theta);
                static Complex Exp(Complex z);
                static Complex Log(Complex z);                  // principal log
                static Complex Pow(Complex z, Complex w);       // z^w = exp(w*log(z))
                static Complex PowReal(Complex z, double a);    // z^a
                static Complex Sqrt(Complex z);                 // principal sqrt
                static Complex Sin(Complex z);
                static Complex Cos(Complex z);
                static Complex Tan(Complex z);

                // roots of unity: n >= 1
                static array<Complex>^ RootsOfUnity(int n);

                // overrides
                virtual String^ ToString() override;
                virtual bool Equals(Object^ obj) override;
                virtual int GetHashCode() override;
            };

        } // namespace Advanced
    } // namespace Math
} // namespace Agrine

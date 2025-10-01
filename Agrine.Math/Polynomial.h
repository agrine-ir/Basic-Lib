#pragma once

// Polynomial.h - Algebra module
// Part of Agrine.Math (C++/CLI, .NET Framework 4.6.2)

#include "Types.h"
#include "Constants.h"
#include "Utils.h"
#include "Exceptions.h"

using namespace System;
using namespace Agrine::Math::Core;

namespace Agrine {
    namespace Math {
        namespace Algebra {

            /// <summary>
            /// Represents a polynomial with real coefficients.
            /// Coefficients are stored in increasing order:
            /// P(x) = coeffs[0] + coeffs[1]*x + ... + coeffs[n]*x^n
            /// </summary>
            public ref class Polynomial
            {
            private:
                array<double>^ coeffs;

                void Trim(); // remove leading zeros

            public:
                // ===== Properties =====
                property int Degree {
                    int get();
                }

                property double default[int] {
                    double get(int i);
                    void set(int i, double value);
                    }

                    // ===== Constructors =====
                Polynomial(... array<double>^ coefficients);
                Polynomial(array<double>^ coefficients);

                // ===== Utility =====
                Polynomial^ Clone();

                // ===== Basic operations =====
                double Evaluate(double x);
                Polynomial^ Add(Polynomial^ other);
                Polynomial^ Subtract(Polynomial^ other);
                Polynomial^ Multiply(Polynomial^ other);
                void Divide(Polynomial^ divisor, Polynomial^% quotient, Polynomial^% remainder);

                // ===== Advanced operations =====
                Polynomial^ Derivative();
                // Remove default argument and add overload
                Polynomial^ Integral(double constantTerm);
                Polynomial^ Integral(); // overload that uses 0.0 as default
                array<double>^ Roots(); // quadratic closed form, Newton otherwise

                static bool AreEqual(Polynomial^ a, Polynomial^ b, double tol = Constants::Tolerance);

                // ===== ToString override =====
                virtual String^ ToString() override;
            };
        }
    }
}

#pragma once

// Numerics.h - Numerical analysis algorithms
// Part of Agrine.Math (C++/CLI, .NET Framework 4.6.2)

#include "Constants.h"
#include "Utils.h"
#include "Exceptions.h"
#include "Polynomial.h"
#include "Matrix.h"

using namespace System;
using namespace Agrine::Math::Core;
using namespace Agrine::Math::Algebra;

namespace Agrine {
    namespace Math {
        namespace Advanced {

            /// <summary>
            /// Numerics - Numerical analysis utilities (root finding, linear systems, differentiation, integration).
            /// </summary>
            public ref class Numerics abstract sealed
            {
            public:
                // ---------- Root finding ----------
                static double Bisection(Func<double, double>^ f, double a, double b, double tol = Constants::Tolerance, int maxIter = 1000);
                static double NewtonRaphson(Func<double, double>^ f, Func<double, double>^ df, double x0, double tol = Constants::Tolerance, int maxIter = 100);
                static double Secant(Func<double, double>^ f, double x0, double x1, double tol = Constants::Tolerance, int maxIter = 100);

                // Polynomial roots (1 real root with Newton)
                static double PolynomialRoot(Polynomial^ p, double x0, double tol = Constants::Tolerance, int maxIter = 100);

                // ---------- Linear systems ----------
                static array<double>^ SolveGaussian(Matrix^ A, array<double>^ b);
                static Tuple<Matrix^, Matrix^>^ LUDecomposition(Matrix^ A); // returns (L,U)
                static array<double>^ SolveLU(Matrix^ L, Matrix^ U, array<double>^ b);

                // ---------- Numerical differentiation ----------
                static double DerivativeForward(Func<double, double>^ f, double x, double h = 1e-5);
                static double DerivativeCentral(Func<double, double>^ f, double x, double h = 1e-5);

                // ---------- Numerical integration ----------
                static double Trapezoidal(Func<double, double>^ f, double a, double b, int n = 1000);
                static double Simpson(Func<double, double>^ f, double a, double b, int n = 1000);
            };

        } // namespace Advanced
    } // namespace Math
} // namespace Agrine

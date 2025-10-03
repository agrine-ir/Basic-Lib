#pragma once

// Optimization.h - Optimization algorithms (1D & multivariate)
// Part of Agrine.Math (C++/CLI, .NET Framework 4.6.2)

#include "Constants.h"
#include "Utils.h"
#include "Exceptions.h"
#include "Matrix.h"
#include "Numerics.h" // optional helpers

using namespace System;

namespace Agrine {
    namespace Math {
        namespace Advanced {

            /// <summary>
            /// OptimizationTools - contains common optimization routines:
            ///  - Golden section search (1D)
            ///  - Newton for extrema (1D)
            ///  - Gradient descent (multidimensional)
            ///  - Conjugate gradient (linear system and quadratic minimization)
            /// </summary>
            public ref class OptimizationTools abstract sealed
            {
            public:
                // ---------- 1D optimization ----------
                // Golden-section search for unimodal function on [a,b]
                static double GoldenSection(Func<double, double>^ f, double a, double b, double tol = Core::Constants::Tolerance, int maxIter = 1000);

                // Newton method for extremum: requires f and f' (df). If second derivative df2==nullptr,
                // numerical second derivative is used.
                static double NewtonForExtremum(Func<double, double>^ f, Func<double, double>^ df, Func<double, double>^ df2, double x0, double tol = Core::Constants::Tolerance, int maxIter = 100);

                // ---------- Multivariate optimization ----------
                // GradientDescent: minimize f(x) using gradient g(x).
                // - x0: initial guess
                // - alpha: initial step size (if useLineSearch==true, acts as initial step)
                // - useLineSearch: if true, perform simple backtracking Armijo line search
                static array<double>^ GradientDescent(
                    Func<array<double>^, double>^ f,
                    Func<array<double>^, array<double>^>^ grad,
                    array<double>^ x0,
                    double alpha = 1e-2,
                    bool useLineSearch = true,
                    double tol = Core::Constants::Tolerance,
                    int maxIter = 1000);

                // ---------- Conjugate Gradient ----------
                // Solve A x = b for symmetric positive-definite A using Conjugate Gradient.
                // Matrix class is used (A[i,j] indexable).
                static array<double>^ ConjugateGradientLinear(Agrine::Math::Algebra::Matrix^ A, array<double>^ b, array<double>^ x0, double tol = 1e-10, int maxIter = 1000);

                // Minimize quadratic f(x) = 0.5 x^T A x - b^T x   (A symmetric positive-definite)
                // returns minimizer x (wrapper around ConjugateGradientLinear)
                static array<double>^ ConjugateGradientMinimize(Agrine::Math::Algebra::Matrix^ A, array<double>^ b, array<double>^ x0, double tol = 1e-10, int maxIter = 1000);

            private:
                // helpers
                static double Dot(array<double>^ a, array<double>^ b);
                static void Axpy(double alpha, array<double>^ x, array<double>^ y); // y += alpha * x
                static double Norm(array<double>^ a);
                static array<double>^ CloneVec(array<double>^ v);
                static void CheckVectorSameLength(array<double>^ a, array<double>^ b, String^ nameA, String^ nameB);
            };

        } // namespace Advanced
    } // namespace Math
} // namespace Agrine

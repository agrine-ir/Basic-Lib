#pragma once

// Integration.h - Calculus module (numerical and analytic integration)
// Part of Agrine.Math (C++/CLI, .NET Framework 4.6.2)

#include "Types.h"
#include "Constants.h"
#include "Utils.h"
#include "Exceptions.h"
#include "Polynomial.h"

using namespace System;
using namespace Agrine::Math::Core;
using namespace Agrine::Math::Algebra;

namespace Agrine {
    namespace Math {
        namespace Calculus {

            /// <summary>
            /// Integration utilities: analytic (Polynomial) and numerical methods.
            /// </summary>
            public ref class Integration abstract sealed
            {
            public:
                // ===== Analytic =====
                /// <summary>Return indefinite integral of polynomial with specified constant term.</summary>
                static Polynomial^ AnalyticIntegral(Polynomial^ poly, double constantTerm);
                /// <summary>Return indefinite integral of polynomial with zero constant.</summary>
                static Polynomial^ AnalyticIntegral(Polynomial^ poly);

                // ===== Composite (fixed-step) numerical integration =====
                /// <summary>Composite rectangle rule (left endpoints).</summary>
                static double Rectangle(Func<double, double>^ f, double a, double b, int n);
                /// <summary>Composite trapezoid rule.</summary>
                static double Trapezoid(Func<double, double>^ f, double a, double b, int n);
                /// <summary>Composite Simpson's rule. n must be even.</summary>
                static double Simpson(Func<double, double>^ f, double a, double b, int n);

                // ===== Adaptive numerical integration =====
                /// <summary>Adaptive Simpson's integration (tolerance-controlled).</summary>
                static double AdaptiveSimpson(Func<double, double>^ f, double a, double b, double tol);

                // ===== Multiple integrals (basic iterative) =====
                /// <summary>Simple iterated 2D integral using composite trapezoid rule.</summary>
                static double Integrate2D(Func<double, double, double>^ f, double ax, double bx, int nx, double ay, double by, int ny);

            private:
                // Private helpers for adaptive Simpson
                static double SimpsonEstimate(Func<double, double>^ f, double a, double b);
                static double AdaptiveSimpsonRecursive(Func<double, double>^ f, double a, double b, double eps, double whole, int depth);
                literal static int MAX_ADAPTIVE_RECURSION = 20;
            };
        }
    }
}

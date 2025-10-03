#pragma once

// Laplace.h - Laplace Transform utilities
// Part of Agrine.Math (C++/CLI, .NET Framework 4.6.2)

#include "Constants.h"
#include "Utils.h"
#include "Exceptions.h"
#include "Complex.h"
#include "Integration.h"

using namespace System;
using namespace Agrine::Math::Core;
using namespace Agrine::Math::Advanced;
using namespace Agrine::Math::Calculus;

namespace Agrine {
    namespace Math {
        namespace Advanced {

            /// <summary>
            /// Laplace Transform utilities (numeric & analytical approximations).
            /// </summary>
            public ref class Laplace abstract sealed
            {
            public:
                /// <summary>
                /// Numeric Laplace transform for sampled real function f(t) on [0,T].
                /// Returns F(s) ≈ ∫0^T f(t) e^{-st} dt
                /// </summary>
                static Complex LaplaceTransform(
                    Func<double, double>^ f,
                    Complex s,
                    double T,
                    int steps = 1000);

                /// <summary>
                /// Approximate inverse Laplace transform using Bromwich integral.
                /// </summary>
                static double InverseLaplaceTransform(
                    Func<Complex, Complex>^ F,
                    double t,
                    double sigma = 1.0,
                    int N = 50);

                /// <summary>
                /// Some standard Laplace transforms.
                /// </summary>
                static Complex TransformExpDecay(double a, Complex s);  // f(t) = e^{-a t}, t >= 0
                static Complex TransformSin(double w, Complex s);       // f(t) = sin(wt), t >= 0
                static Complex TransformCos(double w, Complex s);       // f(t) = cos(wt), t >= 0
                static Complex TransformUnitStep(Complex s);            // f(t) = 1, t >= 0
            };

        } // namespace Advanced
    } // namespace Math
} // namespace Agrine

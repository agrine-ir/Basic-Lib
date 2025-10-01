#pragma once

// DiffEq.h - Calculus module (ODE solvers)
// Part of Agrine.Math (C++/CLI, .NET Framework 4.6.2)

#include "Types.h"
#include "Constants.h"
#include "Utils.h"
#include "Exceptions.h"
#include "Vector.h"

using namespace System;
using namespace System::Collections::Generic;
using namespace Agrine::Math::Core;
using namespace Agrine::Math::Algebra;

namespace Agrine {
    namespace Math {
        namespace Calculus {

            /// <summary>
            /// Scalar ODE solution container: Times and Values arrays.
            /// </summary>
            public ref class ODESolutionScalar
            {
            public:
                array<double>^ Times;
                array<double>^ Values;

                ODESolutionScalar(array<double>^ t, array<double>^ y) {
                    Times = t; Values = y;
                }
            };

            /// <summary>
            /// System ODE solution container: Times and States (vector at each time).
            /// </summary>
            public ref class ODESolutionSystem
            {
            public:
                array<double>^ Times;
                array<Vector^>^ States;

                ODESolutionSystem(array<double>^ t, array<Vector^>^ s) {
                    Times = t; States = s;
                }
            };

            /// <summary>
            /// DiffEq - collection of ODE solvers (scalar and system).
            /// </summary>
            public ref class DiffEq abstract sealed
            {
            public:
                // ============================
                // Scalar ODE solvers
                // ============================

                /// <summary>
                /// Euler explicit solver for dy/dt = f(t,y)
                /// returns arrays of times and y-values.
                /// </summary>
                static ODESolutionScalar^ EulerScalar(Func<double, double, double>^ f,
                    double t0, double y0, double tf, double h);

                /// <summary>
                /// Improved Euler (Heun) for scalar ODE.
                /// </summary>
                static ODESolutionScalar^ HeunScalar(Func<double, double, double>^ f,
                    double t0, double y0, double tf, double h);

                /// <summary>
                /// Classical RK4 for scalar ODE.
                /// </summary>
                static ODESolutionScalar^ RK4Scalar(Func<double, double, double>^ f,
                    double t0, double y0, double tf, double h);

                /// <summary>
                /// Adaptive RKF45 (scalar) with tolerance controlling local error.
                /// </summary>
                static ODESolutionScalar^ AdaptiveRK45Scalar(Func<double, double, double>^ f,
                    double t0, double y0, double tf,
                    double h0, double tol);

                // ============================
                // System ODE solvers (Vector)
                // f: (t, Vector) -> Vector (d/dt state)
                // ============================

                /// <summary>
                /// Euler explicit for system.
                /// </summary>
                static ODESolutionSystem^ EulerSystem(Func<double, Vector^, Vector^>^ f,
                    double t0, Vector^ y0, double tf, double h);

                /// <summary>
                /// Heun (improved Euler) for system.
                /// </summary>
                static ODESolutionSystem^ HeunSystem(Func<double, Vector^, Vector^>^ f,
                    double t0, Vector^ y0, double tf, double h);

                /// <summary>
                /// RK4 for system.
                /// </summary>
                static ODESolutionSystem^ RK4System(Func<double, Vector^, Vector^>^ f,
                    double t0, Vector^ y0, double tf, double h);

                /// <summary>
                /// Adaptive RKF45 for system; uses vector norm to control error.
                /// </summary>
                static ODESolutionSystem^ AdaptiveRK45System(Func<double, Vector^, Vector^>^ f,
                    double t0, Vector^ y0, double tf,
                    double h0, double tol);

            private:
                // small helpers
                static double SafeStepCount(double a, double b, double h);
                static double Max(double a, double b) { return (a > b) ? a : b; }
            };
        }
    }
}

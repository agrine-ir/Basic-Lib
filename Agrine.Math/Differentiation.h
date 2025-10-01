#pragma once

// Differentiation.h - Calculus module
// Part of Agrine.Math (C++/CLI, .NET Framework 4.6.2)

#include "Types.h"
#include "Constants.h"
#include "Utils.h"
#include "Exceptions.h"
#include "Polynomial.h"
#include "Vector.h"

using namespace System;
using namespace Agrine::Math::Core;
using namespace Agrine::Math::Algebra;

namespace Agrine {
    namespace Math {
        namespace Calculus {

            /// <summary>
            /// Provides methods for analytic and numerical differentiation.
            /// </summary>
            public ref class Differentiation abstract sealed
            {
            public:
                // ===== Analytic differentiation =====
                static Polynomial^ Derivative(Polynomial^ poly, int order = 1);

                // ===== Numerical differentiation (single-variable) =====
                static double ForwardDiff(Func<double, double>^ f, double x, double h = Constants::Epsilon);
                static double BackwardDiff(Func<double, double>^ f, double x, double h = Constants::Epsilon);
                static double CentralDiff(Func<double, double>^ f, double x, double h = Constants::Epsilon);

                static double Derivative(Func<double, double>^ f, double x, int order, double h = Constants::Epsilon);

                // ===== Multi-variable differentiation =====
                static Vector^ Gradient(Func<Vector^, double>^ f, Vector^ point, double h = Constants::Epsilon);
                static double PartialDerivative(Func<Vector^, double>^ f, Vector^ point, int variableIndex, double h = Constants::Epsilon);
            };
        }
    }
}

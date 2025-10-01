#pragma once

// Series.h - Calculus module
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
            /// Provides utilities for infinite series, Taylor expansions, and convergence tests.
            /// </summary>
            public ref class Series abstract sealed
            {
            public:
                // ===== Taylor / Maclaurin expansions =====
                static Polynomial^ Maclaurin(Func<double, double>^ f, int order, double h = Constants::Epsilon);
                static Polynomial^ Taylor(Func<double, double>^ f, double a, int order, double h = Constants::Epsilon);

                // ===== Series summation =====
                static double PartialSum(Func<int, double>^ term, int n);
                static double InfiniteSum(Func<int, double>^ term, int maxTerms, double tol = Constants::Tolerance);
                static double EulerTransform(Func<int, double>^ term, int maxTerms);

                // ===== Convergence tests =====
                static bool RatioTest(Func<int, double>^ term, int n, double tol = Constants::Tolerance);
                static bool RootTest(Func<int, double>^ term, int n, double tol = Constants::Tolerance);
                static bool AlternatingSeriesTest(Func<int, double>^ term, int n, double tol = Constants::Tolerance);
            };
        }
    }
}

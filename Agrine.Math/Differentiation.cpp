#include "Differentiation.h"

using namespace System;
using namespace Agrine::Math::Core;
using namespace Agrine::Math::Algebra;

namespace Agrine {
    namespace Math {
        namespace Calculus {

            // ===== Analytic differentiation =====
            Polynomial^ Differentiation::Derivative(Polynomial^ poly, int order)
            {
                if (order < 1)
                    throw gcnew InvalidArgumentException("Order of derivative must be >= 1");

                Polynomial^ result = poly->Clone();
                for (int i = 0; i < order; i++)
                    result = result->Derivative();

                return result;
            }

            // ===== Numerical differentiation (single-variable) =====
            double Differentiation::ForwardDiff(Func<double, double>^ f, double x, double h)
            {
                return (f(x + h) - f(x)) / h;
            }

            double Differentiation::BackwardDiff(Func<double, double>^ f, double x, double h)
            {
                return (f(x) - f(x - h)) / h;
            }

            double Differentiation::CentralDiff(Func<double, double>^ f, double x, double h)
            {
                return (f(x + h) - f(x - h)) / (2.0 * h);
            }

            double Differentiation::Derivative(Func<double, double>^ f, double x, int order, double h)
            {
                if (order < 1)
                    throw gcnew InvalidArgumentException("Order of derivative must be >= 1");

                if (order == 1)
                    return CentralDiff(f, x, h);

                // Higher order derivatives using recursive central difference
                // (not very efficient, but works as a starting point)
                return (Derivative(f, x + h, order - 1, h) - Derivative(f, x - h, order - 1, h)) / (2.0 * h);
            }

            // ===== Multi-variable differentiation =====
            double Differentiation::PartialDerivative(Func<Vector^, double>^ f, Vector^ point, int variableIndex, double h)
            {
                if (variableIndex < 0 || variableIndex >= point->Dimension)
                    throw gcnew IndexOutOfRangeException("Invalid variable index for partial derivative.");

                Vector^ pPlus = point->Clone();
                Vector^ pMinus = point->Clone();

                pPlus[variableIndex] += h;
                pMinus[variableIndex] -= h;

                return (f(pPlus) - f(pMinus)) / (2.0 * h);
            }

            Vector^ Differentiation::Gradient(Func<Vector^, double>^ f, Vector^ point, double h)
            {
                array<double>^ grad = gcnew array<double>(point->Dimension);

                for (int i = 0; i < point->Dimension; i++)
                    grad[i] = PartialDerivative(f, point, i, h);

                return gcnew Vector(grad);
            }

        } // namespace Calculus
    } // namespace Math
} // namespace Agrine

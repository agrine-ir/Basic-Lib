#include "Integration.h"

using namespace System;
using namespace Agrine::Math::Core;
using namespace Agrine::Math::Algebra;

namespace Agrine {
    namespace Math {
        namespace Calculus {

            // ===== Analytic =====
            Polynomial^ Integration::AnalyticIntegral(Polynomial^ poly, double constantTerm)
            {
                if (poly == nullptr)
                    throw gcnew ArgumentNullException("poly");

                return poly->Integral(constantTerm);
            }

            Polynomial^ Integration::AnalyticIntegral(Polynomial^ poly)
            {
                return AnalyticIntegral(poly, 0.0);
            }

            // ===== Composite (fixed-step) numerical integration =====
            double Integration::Rectangle(Func<double, double>^ f, double a, double b, int n)
            {
                if (f == nullptr) throw gcnew ArgumentNullException("f");
                if (n <= 0) throw gcnew InvalidArgumentException("n must be positive.");
                if (Utils::IsZero(b - a)) return 0.0;

                double h = (b - a) / n;
                double sum = 0.0;
                double x = a;
                for (int i = 0; i < n; i++)
                {
                    sum += f(x); // left-endpoint rule
                    x += h;
                }
                return sum * h;
            }

            double Integration::Trapezoid(Func<double, double>^ f, double a, double b, int n)
            {
                if (f == nullptr) throw gcnew ArgumentNullException("f");
                if (n <= 0) throw gcnew InvalidArgumentException("n must be positive.");
                if (Utils::IsZero(b - a)) return 0.0;

                double h = (b - a) / n;
                double sum = 0.5 * (f(a) + f(b));
                for (int i = 1; i < n; i++)
                {
                    double x = a + i * h;
                    sum += f(x);
                }
                return sum * h;
            }

            double Integration::Simpson(Func<double, double>^ f, double a, double b, int n)
            {
                if (f == nullptr) throw gcnew ArgumentNullException("f");
                if (n <= 0) throw gcnew InvalidArgumentException("n must be positive.");
                if (n % 2 != 0) throw gcnew InvalidArgumentException("n must be even for Simpson's rule.");
                if (Utils::IsZero(b - a)) return 0.0;

                double h = (b - a) / n;
                double sum = f(a) + f(b);

                // odd terms
                for (int i = 1; i < n; i += 2) {
                    double x = a + i * h;
                    sum += 4.0 * f(x);
                }
                // even terms
                for (int i = 2; i < n; i += 2) {
                    double x = a + i * h;
                    sum += 2.0 * f(x);
                }
                return sum * (h / 3.0);
            }

            // ===== Adaptive Simpson =====
            double Integration::SimpsonEstimate(Func<double, double>^ f, double a, double b)
            {
                double c = 0.5 * (a + b);
                return (f(a) + 4.0 * f(c) + f(b)) * (b - a) / 6.0;
            }

            double Integration::AdaptiveSimpson(Func<double, double>^ f, double a, double b, double tol)
            {
                if (f == nullptr) throw gcnew ArgumentNullException("f");
                if (tol <= 0.0) throw gcnew InvalidArgumentException("tol must be positive.");
                if (Utils::IsZero(b - a)) return 0.0;

                double whole = SimpsonEstimate(f, a, b);
                return AdaptiveSimpsonRecursive(f, a, b, tol, whole, 0);
            }

            double Integration::AdaptiveSimpsonRecursive(Func<double, double>^ f, double a, double b, double eps, double whole, int depth)
            {
                // Prevent runaway recursion
                if (depth > MAX_ADAPTIVE_RECURSION)
                    throw gcnew ConvergenceException("Adaptive Simpson exceeded maximum recursion depth.");

                double c = 0.5 * (a + b);
                double left = SimpsonEstimate(f, a, c);
                double right = SimpsonEstimate(f, c, b);
                double delta = left + right - whole;

                if (System::Math::Abs(delta) <= 15.0 * eps)
                    return left + right + delta / 15.0; // Richardson extrapolation
                else
                {
                    double leftRes = AdaptiveSimpsonRecursive(f, a, c, eps / 2.0, left, depth + 1);
                    double rightRes = AdaptiveSimpsonRecursive(f, c, b, eps / 2.0, right, depth + 1);
                    return leftRes + rightRes;
                }
            }

            // ===== 2D iterated integral (simple composite trapezoid) =====
            double Integration::Integrate2D(Func<double, double, double>^ f,
                double ax, double bx, int nx,
                double ay, double by, int ny)
            {
                if (f == nullptr) throw gcnew ArgumentNullException("f");
                if (nx <= 0 || ny <= 0) throw gcnew InvalidArgumentException("nx and ny must be positive.");
                if (Utils::IsZero(bx - ax) || Utils::IsZero(by - ay)) return 0.0;

                double hx = (bx - ax) / nx;
                double hy = (by - ay) / ny;

                double sum = 0.0;
                for (int i = 0; i <= nx; i++)
                {
                    double x = ax + i * hx;
                    double wx = (i == 0 || i == nx) ? 0.5 : 1.0;

                    for (int j = 0; j <= ny; j++)
                    {
                        double y = ay + j * hy;
                        double wy = (j == 0 || j == ny) ? 0.5 : 1.0;

                        sum += wx * wy * f(x, y);
                    }
                }

                return sum * (hx * hy);
            }

        } // namespace Calculus
    } // namespace Math
} // namespace Agrine

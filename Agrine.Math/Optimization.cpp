#include "Optimization.h"
#include <cmath>

using namespace System;
using namespace Agrine::Math::Core;
using namespace Agrine::Math::Algebra;

namespace Agrine {
    namespace Math {
        namespace Advanced {

            // --------------------------
            // Private helpers
            // --------------------------
            double OptimizationTools::Dot(array<double>^ a, array<double>^ b)
            {
                CheckVectorSameLength(a, b, "a", "b");
                double s = 0.0;
                for (int i = 0; i < a->Length; ++i) s += a[i] * b[i];
                return s;
            }

            void OptimizationTools::Axpy(double alpha, array<double>^ x, array<double>^ y)
            {
                if (x == nullptr || y == nullptr) throw gcnew ArgumentNullException("x or y");
                if (x->Length != y->Length) throw gcnew InvalidArgumentException("Axpy: vector size mismatch.");
                for (int i = 0; i < x->Length; ++i) y[i] += alpha * x[i];
            }

            double OptimizationTools::Norm(array<double>^ a)
            {
                double s = 0.0;
                for (int i = 0; i < a->Length; ++i) s += a[i] * a[i];
                return System::Math::Sqrt(System::Math::Max(0.0, s));
            }

            array<double>^ OptimizationTools::CloneVec(array<double>^ v)
            {
                if (v == nullptr) return nullptr;
                array<double>^ out = gcnew array<double>(v->Length);
                for (int i = 0; i < v->Length; ++i) out[i] = v[i];
                return out;
            }

            void OptimizationTools::CheckVectorSameLength(array<double>^ a, array<double>^ b, String^ nameA, String^ nameB)
            {
                if (a == nullptr) throw gcnew ArgumentNullException(nameA);
                if (b == nullptr) throw gcnew ArgumentNullException(nameB);
                if (a->Length != b->Length) throw gcnew InvalidArgumentException(String::Format("{0} and {1} must have same length.", nameA, nameB));
            }

            // --------------------------
            // Golden section (1D)
            // --------------------------
            double OptimizationTools::GoldenSection(Func<double, double>^ f, double a, double b, double tol, int maxIter)
            {
                if (f == nullptr) throw gcnew ArgumentNullException("f");
                if (!(a < b)) throw gcnew InvalidArgumentException("GoldenSection: require a < b.");
                if (tol <= 0.0) tol = Core::Constants::Tolerance;
                const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
                double c = b - (b - a) / phi;
                double d = a + (b - a) / phi;
                double fc = f(c);
                double fd = f(d);
                int iter = 0;
                while ((b - a) > tol && iter < maxIter) {
                    if (fc < fd) {
                        b = d;
                        d = c;
                        fd = fc;
                        c = b - (b - a) / phi;
                        fc = f(c);
                    }
                    else {
                        a = c;
                        c = d;
                        fc = fd;
                        d = a + (b - a) / phi;
                        fd = f(d);
                    }
                    ++iter;
                }
                return 0.5 * (a + b);
            }

            // --------------------------
            // Newton for extremum (1D)
            // --------------------------
            double OptimizationTools::NewtonForExtremum(Func<double, double>^ f, Func<double, double>^ df, Func<double, double>^ df2, double x0, double tol, int maxIter)
            {
                if (f == nullptr) throw gcnew ArgumentNullException("f");
                if (df == nullptr) throw gcnew ArgumentNullException("df");
                if (tol <= 0.0) tol = Core::Constants::Tolerance;
                double x = x0;
                for (int iter = 0; iter < maxIter; ++iter) {
                    double fxp = df(x); // derivative
                    double d2;
                    if (df2 != nullptr) d2 = df2(x);
                    else {
                        // numerical second derivative (central)
                        double h = 1e-6;
                        double d1 = df(x + h);
                        double d0 = df(x - h);
                        d2 = (d1 - d0) / (2.0 * h);
                    }
                    if (Utils::IsZero(d2)) throw gcnew InvalidArgumentException("NewtonForExtremum: second derivative near zero.");
                    double x1 = x - fxp / d2;
                    if (System::Math::Abs(x1 - x) < tol) return x1;
                    x = x1;
                }
                return x;
            }

            // --------------------------
            // Gradient Descent (multivariate)
            // --------------------------
            array<double>^ OptimizationTools::GradientDescent(
                Func<array<double>^, double>^ f,
                Func<array<double>^, array<double>^>^ grad,
                array<double>^ x0,
                double alpha,
                bool useLineSearch,
                double tol,
                int maxIter)
            {
                if (f == nullptr) throw gcnew ArgumentNullException("f");
                if (grad == nullptr) throw gcnew ArgumentNullException("grad");
                if (x0 == nullptr) throw gcnew ArgumentNullException("x0");
                if (tol <= 0.0) tol = Core::Constants::Tolerance;

                int n = x0->Length;
                array<double>^ x = CloneVec(x0);
                array<double>^ g = grad(x);
                if (g == nullptr || g->Length != n) throw gcnew InvalidArgumentException("grad must return gradient of same length as x0.");

                for (int iter = 0; iter < maxIter; ++iter) {
                    double gnorm = Norm(g);
                    if (gnorm < tol) return x;

                    double step = alpha;
                    if (useLineSearch) {
                        // backtracking Armijo line search
                        double c = 1e-4;
                        double rho = 0.5;
                        double fx = f(x);
                        // direction is negative gradient
                        array<double>^ pk = CloneVec(g);
                        for (int i = 0; i < n; ++i) pk[i] = -pk[i];
                        // try step sizes
                        int ls_iter = 0;
                        while (ls_iter < 50) {
                            // x_new = x + step * pk
                            array<double>^ xnew = CloneVec(x);
                            for (int i = 0; i < n; ++i) xnew[i] += step * pk[i];
                            double fxnew = f(xnew);
                            // Armijo: f(x+alpha p) <= f(x) + c*alpha*grad^T p
                            double rhs = fx + c * step * Dot(g, pk);
                            if (fxnew <= rhs) { // sufficient decrease
                                // accept
                                x = xnew;
                                break;
                            }
                            step *= rho;
                            ++ls_iter;
                        }
                        if (ls_iter >= 50) {
                            // fallback to small step
                            step = 1e-6;
                            for (int i = 0; i < n; ++i) x[i] += step * (-g[i]);
                        }
                    }
                    else {
                        // simple fixed-step descent
                        for (int i = 0; i < n; ++i) x[i] -= step * g[i];
                    }

                    // update gradient
                    g = grad(x);
                    if (g == nullptr || g->Length != n) throw gcnew InvalidArgumentException("grad must return gradient of same length as x0.");
                }

                return x;
            }

            // --------------------------
            // Conjugate Gradient for linear systems A x = b (symmetric positive-definite)
            // --------------------------
            array<double>^ OptimizationTools::ConjugateGradientLinear(Matrix^ A, array<double>^ b, array<double>^ x0, double tol, int maxIter)
            {
                if (A == nullptr) throw gcnew ArgumentNullException("A");
                if (b == nullptr) throw gcnew ArgumentNullException("b");
                int n = b->Length;
                if (A->Rows != n || A->Cols != n) throw gcnew InvalidArgumentException("ConjugateGradientLinear: dimension mismatch.");

                array<double>^ x = (x0 == nullptr) ? gcnew array<double>(n) : CloneVec(x0);
                for (int i = 0; i < n; ++i) { if (x[i] != x[i]) x[i] = 0.0; } // sanitize NaN

                // r = b - A*x
                array<double>^ r = gcnew array<double>(n);
                array<double>^ Ap = gcnew array<double>(n);
                for (int i = 0; i < n; ++i) {
                    double s = 0.0;
                    for (int j = 0; j < n; ++j) s += A[i, j] * x[j];
                    r[i] = b[i] - s;
                }
                array<double>^ p = CloneVec(r);
                double rsold = Dot(r, r);
                if (System::Math::Sqrt(rsold) < tol) return x;

                for (int iter = 0; iter < maxIter; ++iter) {
                    // Ap = A * p
                    for (int i = 0; i < n; ++i) {
                        double s = 0.0;
                        for (int j = 0; j < n; ++j) s += A[i, j] * p[j];
                        Ap[i] = s;
                    }
                    double alpha = rsold / Dot(p, Ap);
                    // x = x + alpha * p
                    for (int i = 0; i < n; ++i) x[i] += alpha * p[i];
                    // r = r - alpha * Ap
                    for (int i = 0; i < n; ++i) r[i] -= alpha * Ap[i];
                    double rsnew = Dot(r, r);
                    if (System::Math::Sqrt(rsnew) < tol) return x;
                    double beta = rsnew / rsold;
                    for (int i = 0; i < n; ++i) p[i] = r[i] + beta * p[i];
                    rsold = rsnew;
                }

                return x;
            }

            // --------------------------
            // Conjugate gradient minimize quadratic f(x) = 0.5 x^T A x - b^T x
            // wrapper around linear CG solving A x = b
            // --------------------------
            array<double>^ OptimizationTools::ConjugateGradientMinimize(Matrix^ A, array<double>^ b, array<double>^ x0, double tol, int maxIter)
            {
                // for quadratic minimization, stationary condition: A x = b
                return ConjugateGradientLinear(A, b, x0, tol, maxIter);
            }

        } // namespace Advanced
    } // namespace Math
} // namespace Agrine

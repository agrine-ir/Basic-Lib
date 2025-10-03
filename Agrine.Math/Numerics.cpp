#include "Numerics.h"
#include <cmath>

using namespace System;
using namespace Agrine::Math::Core;
using namespace Agrine::Math::Algebra;
using namespace Agrine::Math::Advanced;

namespace Agrine {
    namespace Math {
        namespace Advanced {

            // ---------- Root finding ----------
            double Numerics::Bisection(Func<double, double>^ f, double a, double b, double tol, int maxIter)
            {
                double fa = f(a), fb = f(b);
                if (fa * fb > 0) throw gcnew InvalidArgumentException("Bisection: f(a) and f(b) must have opposite signs.");

                double mid = 0;
                for (int i = 0; i < maxIter; ++i) {
                    mid = 0.5 * (a + b);
                    double fm = f(mid);
                    if (System::Math::Abs(fm) < tol || (b - a) / 2 < tol) return mid;
                    if (fa * fm < 0) { b = mid; fb = fm; }
                    else { a = mid; fa = fm; }
                }
                return mid;
            }

            double Numerics::NewtonRaphson(Func<double, double>^ f, Func<double, double>^ df, double x0, double tol, int maxIter)
            {
                double x = x0;
                for (int i = 0; i < maxIter; ++i) {
                    double fx = f(x);
                    double dfx = df(x);
                    if (Utils::IsZero(dfx)) throw gcnew InvalidArgumentException("NewtonRaphson: derivative near zero.");
                    double x1 = x - fx / dfx;
                    if (System::Math::Abs(x1 - x) < tol) return x1;
                    x = x1;
                }
                return x;
            }

            double Numerics::Secant(Func<double, double>^ f, double x0, double x1, double tol, int maxIter)
            {
                double f0 = f(x0), f1 = f(x1);
                for (int i = 0; i < maxIter; ++i) {
                    if (Utils::AreEqual(f1, f0)) throw gcnew InvalidArgumentException("Secant: zero denominator.");
                    double x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
                    if (System::Math::Abs(x2 - x1) < tol) return x2;
                    x0 = x1; f0 = f1;
                    x1 = x2; f1 = f(x1);
                }
                return x1;
            }

            double Numerics::PolynomialRoot(Polynomial^ p, double x0, double tol, int maxIter)
            {
                if (p == nullptr) throw gcnew ArgumentNullException("p");
                Func<double, double>^ f = gcnew Func<double, double>([p](double x) { return p->Evaluate(x); });
                Func<double, double>^ df = gcnew Func<double, double>([p](double x) { return p->Derivative()->Evaluate(x); });
                return NewtonRaphson(f, df, x0, tol, maxIter);
            }

            // ---------- Linear systems ----------
            array<double>^ Numerics::SolveGaussian(Matrix^ A, array<double>^ b)
            {
                if (A == nullptr || b == nullptr) throw gcnew ArgumentNullException();
                int n = A->Rows;
                if (A->Cols != n || b->Length != n) throw gcnew InvalidArgumentException("SolveGaussian: dimension mismatch.");

                // Copy matrix and vector
                Matrix^ M = A->Clone();
                array<double>^ x = gcnew array<double>(n);
                array<double>^ rhs = (array<double>^)b->Clone();

                // Gaussian elimination with partial pivoting
                for (int k = 0; k < n; ++k) {
                    // pivot
                    int maxRow = k;
                    double maxVal = System::Math::Abs(M[k, k]);
                    for (int i = k + 1; i < n; ++i) {
                        double val = System::Math::Abs(M[i, k]);
                        if (val > maxVal) { maxVal = val; maxRow = i; }
                    }
                    if (Utils::IsZero(maxVal)) throw gcnew InvalidArgumentException("Matrix is singular.");

                    // swap rows
                    if (maxRow != k) {
                        for (int j = 0; j < n; ++j) {
                            double tmp = M[k, j]; M[k, j] = M[maxRow, j]; M[maxRow, j] = tmp;
                        }
                        double tmpR = rhs[k]; rhs[k] = rhs[maxRow]; rhs[maxRow] = tmpR;
                    }

                    // eliminate
                    for (int i = k + 1; i < n; ++i) {
                        double factor = M[i, k] / M[k, k];
                        for (int j = k; j < n; ++j) M[i, j] -= factor * M[k, j];
                        rhs[i] -= factor * rhs[k];
                    }
                }

                // back substitution
                for (int i = n - 1; i >= 0; --i) {
                    double sum = rhs[i];
                    for (int j = i + 1; j < n; ++j) sum -= M[i, j] * x[j];
                    x[i] = sum / M[i, i];
                }

                return x;
            }

            Tuple<Matrix^, Matrix^>^ Numerics::LUDecomposition(Matrix^ A)
            {
                if (A == nullptr) throw gcnew ArgumentNullException("A");
                int n = A->Rows;
                if (A->Cols != n) throw gcnew InvalidArgumentException("LUDecomposition: must be square.");

                Matrix^ L = gcnew Matrix(n, n);
                Matrix^ U = gcnew Matrix(n, n);

                for (int i = 0; i < n; ++i) {
                    // Upper
                    for (int k = i; k < n; ++k) {
                        double sum = 0;
                        for (int j = 0; j < i; ++j) sum += L[i, j] * U[j, k];
                        U[i, k] = A[i, k] - sum;
                    }

                    // Lower
                    for (int k = i; k < n; ++k) {
                        if (i == k) L[i, i] = 1.0;
                        else {
                            double sum = 0;
                            for (int j = 0; j < i; ++j) sum += L[k, j] * U[j, i];
                            if (Utils::IsZero(U[i, i])) throw gcnew InvalidArgumentException("LUDecomposition: singular matrix.");
                            L[k, i] = (A[k, i] - sum) / U[i, i];
                        }
                    }
                }
                return gcnew Tuple<Matrix^, Matrix^>(L, U);
            }

            array<double>^ Numerics::SolveLU(Matrix^ L, Matrix^ U, array<double>^ b)
            {
                int n = L->Rows;
                if (L->Cols != n || U->Rows != n || U->Cols != n || b->Length != n)
                    throw gcnew InvalidArgumentException("SolveLU: dimension mismatch.");

                array<double>^ y = gcnew array<double>(n);
                array<double>^ x = gcnew array<double>(n);

                // forward substitution Ly = b
                for (int i = 0; i < n; ++i) {
                    double sum = b[i];
                    for (int j = 0; j < i; ++j) sum -= L[i, j] * y[j];
                    y[i] = sum / L[i, i];
                }

                // back substitution Ux = y
                for (int i = n - 1; i >= 0; --i) {
                    double sum = y[i];
                    for (int j = i + 1; j < n; ++j) sum -= U[i, j] * x[j];
                    x[i] = sum / U[i, i];
                }

                return x;
            }

            // ---------- Differentiation ----------
            double Numerics::DerivativeForward(Func<double, double>^ f, double x, double h)
            {
                return (f(x + h) - f(x)) / h;
            }

            double Numerics::DerivativeCentral(Func<double, double>^ f, double x, double h)
            {
                return (f(x + h) - f(x - h)) / (2.0 * h);
            }

            // ---------- Integration ----------
            double Numerics::Trapezoidal(Func<double, double>^ f, double a, double b, int n)
            {
                if (a > b) { double tmp = a; a = b; b = tmp; }
                double h = (b - a) / n;
                double sum = 0.5 * (f(a) + f(b));
                for (int i = 1; i < n; ++i) sum += f(a + i * h);
                return sum * h;
            }

            double Numerics::Simpson(Func<double, double>^ f, double a, double b, int n)
            {
                if (a > b) { double tmp = a; a = b; b = tmp; }
                if (n % 2 != 0) ++n;
                double h = (b - a) / n;
                double sum = f(a) + f(b);
                for (int i = 1; i < n; ++i) {
                    double x = a + i * h;
                    sum += (i % 2 == 0) ? 2.0 * f(x) : 4.0 * f(x);
                }
                return sum * h / 3.0;
            }

        } // namespace Advanced
    } // namespace Math
} // namespace Agrine

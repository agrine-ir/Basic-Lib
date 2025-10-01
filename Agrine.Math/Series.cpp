#include "Series.h"
#include "Differentiation.h"

using namespace System;
using namespace Agrine::Math::Core;
using namespace Agrine::Math::Algebra;

namespace Agrine {
    namespace Math {
        namespace Calculus {

            // ===== Taylor / Maclaurin expansions =====
            Polynomial^ Series::Maclaurin(Func<double, double>^ f, int order, double h)
            {
                if (f == nullptr) throw gcnew ArgumentNullException("f");
                if (order < 0) throw gcnew InvalidArgumentException("Order must be non-negative.");

                array<double>^ coeffs = gcnew array<double>(order + 1);

                // compute derivatives at 0 using numerical central difference
                for (int k = 0; k <= order; k++)
                {
                    if (k == 0)
                    {
                        coeffs[0] = f(0.0);
                    }
                    else
                    {
                        // use Differentiation::Derivative numerically
                        coeffs[k] = Differentiation::Derivative(f, 0.0, k, h) / Utils::Factorial(k);
                    }
                }
                return gcnew Polynomial(coeffs);
            }

            Polynomial^ Series::Taylor(Func<double, double>^ f, double a, int order, double h)
            {
                if (f == nullptr) throw gcnew ArgumentNullException("f");
                if (order < 0) throw gcnew InvalidArgumentException("Order must be non-negative.");

                array<double>^ coeffs = gcnew array<double>(order + 1);

                for (int k = 0; k <= order; k++)
                {
                    if (k == 0)
                    {
                        coeffs[0] = f(a);
                    }
                    else
                    {
                        coeffs[k] = Differentiation::Derivative(f, a, k, h) / Utils::Factorial(k);
                    }
                }

                // Construct polynomial in (x - a)
                Polynomial^ P = gcnew Polynomial(coeffs);

                return P; // interpretation: P(x-a)
            }

            // ===== Series summation =====
            double Series::PartialSum(Func<int, double>^ term, int n)
            {
                if (term == nullptr) throw gcnew ArgumentNullException("term");
                if (n < 0) throw gcnew InvalidArgumentException("n must be >= 0");

                double sum = 0.0;
                for (int k = 0; k <= n; k++)
                    sum += term(k);
                return sum;
            }

            double Series::InfiniteSum(Func<int, double>^ term, int maxTerms, double tol)
            {
                if (term == nullptr) throw gcnew ArgumentNullException("term");
                if (maxTerms <= 0) throw gcnew InvalidArgumentException("maxTerms must be > 0");

                double sum = 0.0;
                double prev = 0.0;

                for (int k = 0; k < maxTerms; k++)
                {
                    sum += term(k);

                    if (System::Math::Abs(sum - prev) < tol)
                        return sum;

                    prev = sum;
                }
                return sum; // approximate
            }

            double Series::EulerTransform(Func<int, double>^ term, int maxTerms)
            {
                if (term == nullptr) throw gcnew ArgumentNullException("term");
                if (maxTerms <= 0) throw gcnew InvalidArgumentException("maxTerms must be > 0");

                // Euler transform works best for alternating series
                double sum = 0.0;
                for (int k = 0; k < maxTerms; k++)
                {
                    sum += term(k) / System::Math::Pow(2.0, k + 1);
                }
                return 2.0 * sum;
            }

            // ===== Convergence tests =====
            bool Series::RatioTest(Func<int, double>^ term, int n, double tol)
            {
                if (term == nullptr) throw gcnew ArgumentNullException("term");
                if (n <= 0) throw gcnew InvalidArgumentException("n must be > 0");

                double a_n = term(n);
                double a_next = term(n + 1);

                if (Utils::IsZero(a_n)) return false;

                double L = System::Math::Abs(a_next / a_n);
                return (L < 1.0 + tol);
            }

            bool Series::RootTest(Func<int, double>^ term, int n, double tol)
            {
                if (term == nullptr) throw gcnew ArgumentNullException("term");
                if (n < 1) throw gcnew InvalidArgumentException("n must be >= 1");

                double a_n = System::Math::Abs(term(n));
                double root = System::Math::Pow(a_n, 1.0 / n);
                return (root < 1.0 + tol);
            }

            bool Series::AlternatingSeriesTest(Func<int, double>^ term, int n, double tol)
            {
                if (term == nullptr) throw gcnew ArgumentNullException("term");
                if (n <= 0) throw gcnew InvalidArgumentException("n must be > 0");

                double a_n = System::Math::Abs(term(n));
                double a_next = System::Math::Abs(term(n + 1));

                bool decreasing = (a_next <= a_n + tol);
                bool limitZero = (System::Math::Abs(term(n + 50)) < tol); // crude limit check

                return decreasing && limitZero;
            }

        } // namespace Calculus
    } // namespace Math
} // namespace Agrine

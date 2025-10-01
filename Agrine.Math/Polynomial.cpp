#include "Polynomial.h"

using namespace System;
using namespace Agrine::Math::Core;

namespace Agrine {
    namespace Math {
        namespace Algebra {

            // ===== Private =====
            void Polynomial::Trim()
            {
                int deg = coeffs->Length - 1;
                while (deg > 0 && Utils::IsZero(coeffs[deg]))
                    deg--;

                if (deg < coeffs->Length - 1)
                {
                    array<double>^ newCoeffs = gcnew array<double>(deg + 1);
                    Array::Copy(coeffs, newCoeffs, deg + 1);
                    coeffs = newCoeffs;
                }
            }

            // ===== Constructors =====
            Polynomial::Polynomial(... array<double>^ coefficients)
            {
                if (coefficients == nullptr || coefficients->Length == 0)
                    throw gcnew ArgumentException("Coefficients cannot be null or empty.");
                coeffs = gcnew array<double>(coefficients->Length);
                Array::Copy(coefficients, coeffs, coefficients->Length);
                Trim();
            }

            Polynomial::Polynomial(array<double>^ coefficients)
            {
                if (coefficients == nullptr || coefficients->Length == 0)
                    throw gcnew ArgumentException("Coefficients cannot be null or empty.");
                coeffs = gcnew array<double>(coefficients->Length);
                Array::Copy(coefficients, coeffs, coefficients->Length);
                Trim();
            }

            // ===== Properties =====
            int Polynomial::Degree::get()
            {
                return coeffs->Length - 1;
            }

            double Polynomial::default::get(int i)
            {
                if (i < 0 || i >= coeffs->Length)
                    throw gcnew IndexOutOfRangeException("Coefficient index out of range.");
                return coeffs[i];
            }

            void Polynomial::default::set(int i, double value)
            {
                if (i < 0 || i >= coeffs->Length)
                    throw gcnew IndexOutOfRangeException("Coefficient index out of range.");
                coeffs[i] = value;
                Trim();
            }

            // ===== Utility =====
            Polynomial^ Polynomial::Clone()
            {
                return gcnew Polynomial(coeffs);
            }

            // ===== Basic operations =====
            double Polynomial::Evaluate(double x)
            {
                double result = 0.0;
                for (int i = Degree; i >= 0; i--)
                {
                    result = result * x + coeffs[i]; // Horner's method
                }
                return result;
            }

            Polynomial^ Polynomial::Add(Polynomial^ other)
            {
                int maxDeg = System::Math::Max(this->Degree, other->Degree);
                array<double>^ result = gcnew array<double>(maxDeg + 1);

                for (int i = 0; i <= maxDeg; i++)
                {
                    double a = (i <= this->Degree) ? this->coeffs[i] : 0.0;
                    double b = (i <= other->Degree) ? other->coeffs[i] : 0.0;
                    result[i] = a + b;
                }
                return gcnew Polynomial(result);
            }

            Polynomial^ Polynomial::Subtract(Polynomial^ other)
            {
                int maxDeg = System::Math::Max(this->Degree, other->Degree);
                array<double>^ result = gcnew array<double>(maxDeg + 1);

                for (int i = 0; i <= maxDeg; i++)
                {
                    double a = (i <= this->Degree) ? this->coeffs[i] : 0.0;
                    double b = (i <= other->Degree) ? other->coeffs[i] : 0.0;
                    result[i] = a - b;
                }
                return gcnew Polynomial(result);
            }

            Polynomial^ Polynomial::Multiply(Polynomial^ other)
            {
                array<double>^ result = gcnew array<double>(this->Degree + other->Degree + 1);

                for (int i = 0; i <= this->Degree; i++)
                    for (int j = 0; j <= other->Degree; j++)
                        result[i + j] += this->coeffs[i] * other->coeffs[j];

                return gcnew Polynomial(result);
            }

            void Polynomial::Divide(Polynomial^ divisor, Polynomial^% quotient, Polynomial^% remainder)
            {
                if (divisor->Degree < 0 || Utils::IsZero(divisor->coeffs[divisor->Degree]))
                    throw gcnew DivideByZeroException("Cannot divide by zero polynomial.");

                array<double>^ dividend = gcnew array<double>(this->coeffs->Length);
                Array::Copy(this->coeffs, dividend, this->coeffs->Length);

                array<double>^ q = gcnew array<double>(System::Math::Max(0, this->Degree - divisor->Degree + 1));

                for (int i = this->Degree; i >= divisor->Degree; i--)
                {
                    double factor = dividend[i] / divisor->coeffs[divisor->Degree];
                    q[i - divisor->Degree] = factor;

                    for (int j = 0; j <= divisor->Degree; j++)
                        dividend[i - j] -= factor * divisor->coeffs[divisor->Degree - j];
                }

                quotient = gcnew Polynomial(q);

                int remDeg = divisor->Degree - 1;
                if (remDeg < 0) remDeg = 0;

                array<double>^ r = gcnew array<double>(remDeg + 1);
                Array::Copy(dividend, r, r->Length);
                remainder = gcnew Polynomial(r);
            }

            // ===== Advanced operations =====
            Polynomial^ Polynomial::Derivative()
            {
                if (Degree == 0) return gcnew Polynomial(gcnew array<double>{0.0});

                array<double>^ result = gcnew array<double>(Degree);
                for (int i = 1; i <= Degree; i++)
                    result[i - 1] = coeffs[i] * i;

                return gcnew Polynomial(result);
            }

            Polynomial^ Polynomial::Integral(double constantTerm)
            {
                array<double>^ result = gcnew array<double>(Degree + 2);
                result[0] = constantTerm;

                for (int i = 0; i <= Degree; i++)
                    result[i + 1] = coeffs[i] / (i + 1);

                return gcnew Polynomial(result);
            }

            array<double>^ Polynomial::Roots()
            {
                if (Degree == 1)
                {
                    // ax + b = 0 -> root = -b/a
                    double a = coeffs[1], b = coeffs[0];
                    return gcnew array<double>{ -b / a };
                }
                else if (Degree == 2)
                {
                    double a = coeffs[2], b = coeffs[1], c = coeffs[0];
                    double disc = b * b - 4 * a * c;

                    if (disc < 0)
                        return gcnew array<double>(0); // no real roots

                    if (Utils::IsZero(disc))
                        return gcnew array<double>{ -b / (2 * a) };

                    double sqrtDisc = System::Math::Sqrt(disc);
                    return gcnew array<double>{
                        (-b + sqrtDisc) / (2 * a),
                            (-b - sqrtDisc) / (2 * a)
                    };
                }
                else
                {
                    // Placeholder: Newton's method for higher degrees could be added here
                    return gcnew array<double>(0); // Not implemented
                }
            }

            bool Polynomial::AreEqual(Polynomial^ a, Polynomial^ b, double tol)
            {
                if (a->Degree != b->Degree) return false;
                for (int i = 0; i <= a->Degree; i++)
                {
                    if (!Utils::AreEqual(a->coeffs[i], b->coeffs[i], tol))
                        return false;
                }
                return true;
            }

            // ===== ToString =====
            String^ Polynomial::ToString()
            {
                System::Text::StringBuilder^ sb = gcnew System::Text::StringBuilder();

                for (int i = Degree; i >= 0; i--)
                {
                    double coef = coeffs[i];
                    if (Utils::IsZero(coef)) continue;

                    if (sb->Length > 0)
                        sb->Append(coef >= 0 ? " + " : " - ");
                    else if (coef < 0)
                        sb->Append("-");

                    double absCoef = System::Math::Abs(coef);

                    if (!(Utils::AreEqual(absCoef, 1.0) && i > 0))
                        sb->Append(absCoef.ToString("G5"));

                    if (i > 0)
                        sb->Append("x");

                    if (i > 1)
                        sb->Append("^" + i);
                }

                if (sb->Length == 0)
                    sb->Append("0");

                return sb->ToString();
            }

        } // namespace Algebra
    } // namespace Math
} // namespace Agrine

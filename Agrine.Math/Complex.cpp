#include "Complex.h"

using namespace System;
using namespace Agrine::Math::Core;
using namespace Agrine::Math::Advanced;

namespace Agrine {
    namespace Math {
        namespace Advanced {

            // -------------------------
            // Constructors & properties
            // -------------------------
            Complex::Complex(double real, double imag)
            {
                this->Re = real;
                this->Im = imag;
            }

            Complex::Complex(double real)
            {
                this->Re = real;
                this->Im = 0.0;
            }

            double Complex::Real::get() { return Re; }
            void Complex::Real::set(double v) { Re = v; }

            double Complex::Imag::get() { return Im; }
            void Complex::Imag::set(double v) { Im = v; }

            // -------------------------
            // Constants
            // -------------------------
            Complex Complex::Zero() { return Complex(0.0, 0.0); }
            Complex Complex::One() { return Complex(1.0, 0.0); }
            Complex Complex::I() { return Complex(0.0, 1.0); }

            // -------------------------
            // Basic numeric helpers
            // -------------------------
            double Complex::Abs()
            {
                // use hypot for numeric robustness if available; fallback to sqrt
                return std::hypot(Re, Im);
            }

            double Complex::AbsSquared()
            {
                return Re * Re + Im * Im;
            }

            double Complex::Arg()
            {
                return std::atan2(Im, Re);
            }

            Complex Complex::Conjugate()
            {
                return Complex(Re, -Im);
            }

            bool Complex::IsZero()
            {
                return Utils::AreEqual(Re, 0.0) && Utils::AreEqual(Im, 0.0);
            }

            // -------------------------
            // Operators
            // -------------------------
            Complex Complex::operator+(Complex a, Complex b)
            {
                return Complex(a.Re + b.Re, a.Im + b.Im);
            }

            Complex Complex::operator-(Complex a, Complex b)
            {
                return Complex(a.Re - b.Re, a.Im - b.Im);
            }

            Complex Complex::operator-(Complex a)
            {
                return Complex(-a.Re, -a.Im);
            }

            Complex Complex::operator*(Complex a, Complex b)
            {
                // (a+ib)(c+id) = (ac - bd) + i(ad + bc)
                return Complex(a.Re * b.Re - a.Im * b.Im, a.Re * b.Im + a.Im * b.Re);
            }

            Complex Complex::operator/(Complex a, Complex b)
            {
                // stable complex division:
                double br = b.Re, bi = b.Im;
                if (std::fabs(br) >= std::fabs(bi)) {
                    if (Utils::AreEqual(br, 0.0) && Utils::AreEqual(bi, 0.0)) throw gcnew InvalidArgumentException("Division by zero complex.");
                    double r = bi / br;
                    double denom = br + r * bi;
                    double re = (a.Re + a.Im * r) / denom;
                    double im = (a.Im - a.Re * r) / denom;
                    return Complex(re, im);
                }
                else {
                    double r = br / bi;
                    double denom = bi + r * br;
                    double re = (a.Re * r + a.Im) / denom;
                    double im = (a.Im * r - a.Re) / denom;
                    return Complex(re, im);
                }
            }

            Complex Complex::operator*(Complex a, double k)
            {
                return Complex(a.Re * k, a.Im * k);
            }

            Complex Complex::operator*(double k, Complex a)
            {
                return Complex(a.Re * k, a.Im * k);
            }

            Complex Complex::operator/(Complex a, double k)
            {
                if (Utils::IsZero(k)) throw gcnew InvalidArgumentException("Division by zero scalar.");
                return Complex(a.Re / k, a.Im / k);
            }

            // -------------------------
            // Advanced functions
            // -------------------------
            Complex Complex::FromPolar(double r, double theta)
            {
                return Complex(r * std::cos(theta), r * std::sin(theta));
            }

            Complex Complex::Exp(Complex z)
            {
                // exp(x+iy) = e^x (cos y + i sin y)
                double ex = std::exp(z.Re);
                double c = std::cos(z.Im);
                double s = std::sin(z.Im);
                return Complex(ex * c, ex * s);
            }

            Complex Complex::Log(Complex z)
            {
                // principal log: ln|z| + i arg(z) , arg in (-pi, pi]
                double r = z.Abs();
                if (Utils::IsZero(r)) throw gcnew InvalidArgumentException("Log of zero complex number is undefined.");
                double theta = z.Arg();
                return Complex(std::log(r), theta);
            }

            Complex Complex::Pow(Complex z, Complex w)
            {
                // z^w = exp(w * log(z))
                if (z.IsZero() && (!Utils::IsZero(w.Re) || !Utils::IsZero(w.Im))) {
                    // 0^w -> 0 for Re(w)>0 ? but generally undefined for complex exponent; we define 0^0 = 1
                    if (Utils::IsZero(w.Re) && Utils::IsZero(w.Im)) return Complex::One();
                    // if z==0 and w!=0, return 0 (common convention when Re(w)>0); we keep simple: throw if ambiguous
                    if (w.Re > 0.0) return Complex(0.0, 0.0);
                    throw gcnew InvalidArgumentException("0^w is undefined for general complex exponent.");
                }
                Complex L = Complex::Log(z);
                Complex prod = Complex(L.Re * w.Re - L.Im * w.Im, L.Re * w.Im + L.Im * w.Re);
                return Complex::Exp(prod);
            }

            Complex Complex::PowReal(Complex z, double a)
            {
                // z^a = exp(a * log z)
                if (z.IsZero()) {
                    if (Utils::AreEqual(a, 0.0)) return Complex::One(); // 0^0 => 1 convention
                    if (a > 0.0) return Complex(0.0, 0.0);
                    throw gcnew InvalidArgumentException("0^a undefined for negative exponent.");
                }
                Complex L = Complex::Log(z);
                Complex prod = Complex(L.Re * a, L.Im * a);
                return Complex::Exp(prod);
            }

            Complex Complex::Sqrt(Complex z)
            {
                // principal sqrt: sqrt(r) * (cos(theta/2) + i sin(theta/2))
                double r = z.Abs();
                if (Utils::IsZero(r)) return Complex(0.0, 0.0);
                double re = std::sqrt((r + z.Re) / 2.0);
                double im = std::sqrt((r - z.Re) / 2.0);
                if (z.Im < 0) im = -im;
                return Complex(re, im);
            }

            Complex Complex::Sin(Complex z)
            {
                // sin(x+iy) = sin x cosh y + i cos x sinh y
                double x = z.Re, y = z.Im;
                double sx = std::sin(x);
                double cx = std::cos(x);
                double shy = std::sinh(y);
                double chy = std::cosh(y);
                return Complex(sx * chy, cx * shy);
            }

            Complex Complex::Cos(Complex z)
            {
                // cos(x+iy) = cos x cosh y - i sin x sinh y
                double x = z.Re, y = z.Im;
                double cx = std::cos(x);
                double sx = std::sin(x);
                double chy = std::cosh(y);
                double shy = std::sinh(y);
                return Complex(cx * chy, -sx * shy);
            }

            Complex Complex::Tan(Complex z)
            {
                // tan z = sin z / cos z
                Complex s = Complex::Sin(z);
                Complex c = Complex::Cos(z);
                return s / c;
            }

            // -------------------------
            // Roots of unity
            // -------------------------
            array<Complex>^ Complex::RootsOfUnity(int n)
            {
                if (n <= 0) throw gcnew InvalidArgumentException("n must be positive.");
                array<Complex>^ out = gcnew array<Complex>(n);
                double twoPi = 2.0 * Constants::Pi;
                for (int k = 0; k < n; ++k) {
                    double theta = twoPi * (double)k / (double)n;
                    out[k] = Complex::FromPolar(1.0, theta);
                }
                return out;
            }

            // -------------------------
            // Overrides
            // -------------------------
            String^ Complex::ToString()
            {
                // nice formatting with sign
                String^ re = Re.ToString("G6");
                String^ ims = System::Math::Abs(Im).ToString("G6");
                String^ sign = (Im >= 0.0) ? "+" : "-";
                return String::Format("{0} {1} {2}i", re, sign, ims);
            }

            bool Complex::Equals(Object^ obj)
            {
                if (obj == nullptr) return false;
                if (obj->GetType() != Complex::typeid) return false;
                Complex other = safe_cast<Complex>(obj);
                return Utils::AreEqual(this->Re, other.Re) && Utils::AreEqual(this->Im, other.Im);
            }

            int Complex::GetHashCode()
            {
                // combine bit patterns of doubles for a stable hash
                long long bitsRe = System::BitConverter::DoubleToInt64Bits(Re);
                long long bitsIm = System::BitConverter::DoubleToInt64Bits(Im);
                long long combined = bitsRe ^ (bitsIm << 1);
                // fold to int
                int h = (int)(combined ^ (combined >> 32));
                return h;
            }

        } // namespace Advanced
    } // namespace Math
} // namespace Agrine

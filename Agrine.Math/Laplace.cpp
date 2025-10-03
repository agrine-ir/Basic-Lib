#include "Laplace.h"
#include <cmath>

using namespace System;
using namespace Agrine::Math::Core;
using namespace Agrine::Math::Advanced;
using namespace Agrine::Math::Calculus;

namespace Agrine {
    namespace Math {
        namespace Advanced {

            Complex Laplace::LaplaceTransform(
                Func<double, double>^ f,
                Complex s,
                double T,
                int steps)
            {
                if (f == nullptr) throw gcnew ArgumentNullException("f");
                if (T <= 0.0) throw gcnew InvalidArgumentException("T must be > 0.");
                if (steps < 10) steps = 10;

                double h = T / steps;
                Complex sum = Complex::Zero();

                for (int k = 0; k <= steps; ++k) {
                    double t = k * h;
                    double ft = f(t);
                    // weight for trapezoidal rule
                    double w = (k == 0 || k == steps) ? 0.5 : 1.0;
                    Complex kernel = Complex::Exp(-s * t);
                    sum = sum + kernel * ft * w;
                }

                return sum * h;
            }

            double Laplace::InverseLaplaceTransform(
                Func<Complex, Complex>^ F,
                double t,
                double sigma,
                int N)
            {
                if (F == nullptr) throw gcnew ArgumentNullException("F");
                if (t < 0.0) throw gcnew InvalidArgumentException("t must be >= 0.");
                if (N < 10) N = 10;

                // Bromwich integral approximation:
                // f(t) = e^{sigma t}/(2pi) ∫_{-∞}^{∞} F(sigma + iω) e^{iωt} dω
                // approximate using trapezoidal on [-W, W]
                double W = 10.0; // truncation limit
                double h = 2 * W / N;
                Complex sum = Complex::Zero();

                for (int k = 0; k <= N; ++k) {
                    double w = -W + k * h;
                    Complex s = Complex(sigma, w);
                    Complex term = F(s) * Complex::Exp(s * t);
                    double weight = (k == 0 || k == N) ? 0.5 : 1.0;
                    sum = sum + term * weight;
                }

                Complex result = (sum * h) / (2.0 * Constants::Pi);
                return result.Real;
            }

            // Standard transforms
            Complex Laplace::TransformExpDecay(double a, Complex s)
            {
                // L{e^{-a t}} = 1 / (s + a)
                return Complex::One() / (s + Complex(a, 0.0));
            }

            Complex Laplace::TransformSin(double w, Complex s)
            {
                // L{sin(wt)} = w / (s^2 + w^2)
                return Complex(w, 0.0) / (s * s + Complex(w * w, 0.0));
            }

            Complex Laplace::TransformCos(double w, Complex s)
            {
                // L{cos(wt)} = s / (s^2 + w^2)
                return s / (s * s + Complex(w * w, 0.0));
            }

            Complex Laplace::TransformUnitStep(Complex s)
            {
                // L{1} = 1/s
                return Complex::One() / s;
            }

        } // namespace Advanced
    } // namespace Math
} // namespace Agrine

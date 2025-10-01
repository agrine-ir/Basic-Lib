#include "Distributions.h"
#include <math.h>

using namespace System;
using namespace Agrine::Math::Core;

namespace Agrine {
    namespace Math {
        namespace Probability {

            Random^ Distributions::rng = gcnew Random();

            // ---------------- Uniform(a,b) ----------------
            double Distributions::UniformPDF(double x, double a, double b)
            {
                if (a >= b) throw gcnew InvalidArgumentException("Uniform: a must be < b.");
                return (x >= a && x <= b) ? 1.0 / (b - a) : 0.0;
            }

            double Distributions::UniformCDF(double x, double a, double b)
            {
                if (a >= b) throw gcnew InvalidArgumentException("Uniform: a must be < b.");
                if (x < a) return 0.0;
                if (x > b) return 1.0;
                return (x - a) / (b - a);
            }

            double Distributions::UniformSample(double a, double b)
            {
                if (a >= b) throw gcnew InvalidArgumentException("Uniform: a must be < b.");
                double u = rng->NextDouble(); // [0,1)
                return a + u * (b - a);
            }

            // ---------------- Normal(mu,sigma) ----------------
            double Distributions::NormalPDF(double x, double mu, double sigma)
            {
                if (sigma <= 0) throw gcnew InvalidArgumentException("Normal: sigma must be > 0.");
                double z = (x - mu) / sigma;
                return (1.0 / (sigma * System::Math::Sqrt(2 * Constants::Pi))) * System::Math::Exp(-0.5 * z * z);
            }

            // Approximate Normal CDF using error function
            double Distributions::NormalCDF(double x, double mu, double sigma)
            {
                if (sigma <= 0) throw gcnew InvalidArgumentException("Normal: sigma must be > 0.");
                double z = (x - mu) / (sigma * System::Math::Sqrt(2.0));
                return 0.5 * (1.0 + System::Math::Erf(z));
            }

            double Distributions::NormalSample(double mu, double sigma)
            {
                if (sigma <= 0) throw gcnew InvalidArgumentException("Normal: sigma must be > 0.");
                // Box-Muller transform
                double u1 = rng->NextDouble();
                double u2 = rng->NextDouble();
                double r = System::Math::Sqrt(-2.0 * System::Math::Log(u1));
                double theta = 2.0 * Constants::Pi * u2;
                double z = r * System::Math::Cos(theta);
                return mu + sigma * z;
            }

            // ---------------- Binomial(n,p) ----------------
            double Distributions::BinomialPMF(int k, int n, double p)
            {
                if (n < 0) throw gcnew InvalidArgumentException("Binomial: n must be non-negative.");
                if (p < 0 || p > 1) throw gcnew InvalidArgumentException("Binomial: p must be in [0,1].");
                if (k < 0 || k > n) return 0.0;

                // Compute C(n,k)
                double logC = 0;
                for (int i = 1; i <= k; ++i) {
                    logC += System::Math::Log((n - k + i) / (double)i);
                }
                double coeff = System::Math::Exp(logC);
                return coeff * System::Math::Pow(p, k) * System::Math::Pow(1 - p, n - k);
            }

            double Distributions::BinomialCDF(int k, int n, double p)
            {
                double sum = 0.0;
                for (int i = 0; i <= k; ++i) sum += BinomialPMF(i, n, p);
                return sum;
            }

            int Distributions::BinomialSample(int n, double p)
            {
                if (n < 0) throw gcnew InvalidArgumentException("Binomial: n must be non-negative.");
                if (p < 0 || p > 1) throw gcnew InvalidArgumentException("Binomial: p must be in [0,1].");
                int count = 0;
                for (int i = 0; i < n; ++i) {
                    if (rng->NextDouble() < p) count++;
                }
                return count;
            }

            // ---------------- Poisson(lambda) ----------------
            double Distributions::PoissonPMF(int k, double lambda)
            {
                if (lambda <= 0) throw gcnew InvalidArgumentException("Poisson: lambda must be > 0.");
                if (k < 0) return 0.0;
                return System::Math::Exp(-lambda) * System::Math::Pow(lambda, k) / Utils::Factorial(k);
            }

            double Distributions::PoissonCDF(int k, double lambda)
            {
                if (k < 0) return 0.0;
                double sum = 0.0;
                for (int i = 0; i <= k; ++i) sum += PoissonPMF(i, lambda);
                return sum;
            }

            int Distributions::PoissonSample(double lambda)
            {
                if (lambda <= 0) throw gcnew InvalidArgumentException("Poisson: lambda must be > 0.");
                // Knuth's algorithm
                double L = System::Math::Exp(-lambda);
                int k = 0;
                double p = 1.0;
                do {
                    k++;
                    p *= rng->NextDouble();
                } while (p > L);
                return k - 1;
            }

        } // namespace Probability
    } // namespace Math
} // namespace Agrine

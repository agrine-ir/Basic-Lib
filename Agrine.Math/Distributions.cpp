#include "Distributions.h"

#include <cmath>    // std::erf, std::lgamma, std::exp, std::log, std::sqrt, std::pow

using namespace System;
using namespace Agrine::Math::Core;

namespace Agrine {
    namespace Math {
        namespace Probability {

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
                if (sigma <= 0.0) throw gcnew InvalidArgumentException("Normal: sigma must be > 0.");
                double z = (x - mu) / sigma;
                double coeff = 1.0 / (sigma * std::sqrt(2.0 * Constants::Pi));
                return coeff * std::exp(-0.5 * z * z);
            }

            double Distributions::NormalCDF(double x, double mu, double sigma)
            {
                if (sigma <= 0.0) throw gcnew InvalidArgumentException("Normal: sigma must be > 0.");
                // use std::erf (available via <cmath>)
                double z = (x - mu) / (sigma * std::sqrt(2.0));
                return 0.5 * (1.0 + std::erf(z));
            }

            double Distributions::NormalSample(double mu, double sigma)
            {
                if (sigma <= 0.0) throw gcnew InvalidArgumentException("Normal: sigma must be > 0.");
                // Box-Muller
                double u1 = rng->NextDouble();
                double u2 = rng->NextDouble();
                if (u1 < 1e-300) u1 = 1e-300;
                double r = std::sqrt(-2.0 * std::log(u1));
                double theta = 2.0 * Constants::Pi * u2;
                double z = r * std::cos(theta);
                return mu + sigma * z;
            }

            // ---------------- Binomial(n,p) ----------------
            double Distributions::BinomialPMF(int k, int n, double p)
            {
                if (n < 0) throw gcnew InvalidArgumentException("Binomial: n must be non-negative.");
                if (p < 0.0 || p > 1.0) throw gcnew InvalidArgumentException("Binomial: p must be in [0,1].");
                if (k < 0 || k > n) return 0.0;

                // Use log-gamma for robust combination calculation: C(n,k) = exp(lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1))
                double logC = std::lgamma((double)n + 1.0) - std::lgamma((double)k + 1.0) - std::lgamma((double)(n - k) + 1.0);
                double coeff = std::exp(logC);
                return coeff * std::pow(p, (double)k) * std::pow(1.0 - p, (double)(n - k));
            }

            double Distributions::BinomialCDF(int k, int n, double p)
            {
                if (n < 0) throw gcnew InvalidArgumentException("Binomial: n must be non-negative.");
                if (p < 0.0 || p > 1.0) throw gcnew InvalidArgumentException("Binomial: p must be in [0,1].");
                if (k < 0) return 0.0;
                if (k >= n) return 1.0;
                double sum = 0.0;
                for (int i = 0; i <= k; ++i) sum += BinomialPMF(i, n, p);
                return sum;
            }

            int Distributions::BinomialSample(int n, double p)
            {
                if (n < 0) throw gcnew InvalidArgumentException("Binomial: n must be non-negative.");
                if (p < 0.0 || p > 1.0) throw gcnew InvalidArgumentException("Binomial: p must be in [0,1].");
                int count = 0;
                for (int i = 0; i < n; ++i) {
                    if (rng->NextDouble() < p) ++count;
                }
                return count;
            }

            // ---------------- Poisson(lambda) ----------------
            double Distributions::PoissonPMF(int k, double lambda)
            {
                if (lambda <= 0.0) throw gcnew InvalidArgumentException("Poisson: lambda must be > 0.");
                if (k < 0) return 0.0;
                // use log form: exp(-lambda + k*log(lambda) - lgamma(k+1))
                return std::exp(-lambda + k * std::log(lambda) - std::lgamma((double)k + 1.0));
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
                if (lambda <= 0.0) throw gcnew InvalidArgumentException("Poisson: lambda must be > 0.");
                // Knuth's algorithm
                double L = std::exp(-lambda);
                int k = 0;
                double p = 1.0;
                do {
                    ++k;
                    p *= rng->NextDouble();
                } while (p > L);
                return k - 1;
            }

        } // namespace Probability
    } // namespace Math
} // namespace Agrine

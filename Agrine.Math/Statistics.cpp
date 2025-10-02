#include "Statistics.h"
#include <algorithm> // for std::sort
#include <numeric>   // for accumulate
#include <cmath>

using namespace System;
using namespace System::Collections::Generic;
using namespace Agrine::Math::Probability;
using namespace Agrine::Math::Core;

namespace Agrine {
    namespace Math {
        namespace Statistics {

            // -----------------------
            // Helpers
            // -----------------------
            void StatisticsTools::ValidateNonEmpty(array<double>^ a, String^ name)
            {
                if (a == nullptr) throw gcnew ArgumentNullException(name);
                if (a->Length == 0) throw gcnew InvalidArgumentException(name + " must contain at least one element.");
            }

            // convert matrix-of-observations to column vector (rowIndex) - not used widely but kept for extensibility
            array<double>^ StatisticsTools::ToDoubleArrayFromMatrixRow(array<array<double>^>^ data, int rowIndex)
            {
                if (data == nullptr) throw gcnew ArgumentNullException("data");
                int n = data->Length;
                array<double>^ out = gcnew array<double>(n);
                for (int i = 0; i < n; ++i) {
                    if (data[i] == nullptr) throw gcnew InvalidArgumentException("Each observation must be an array.");
                    if (rowIndex < 0 || rowIndex >= data[i]->Length) throw gcnew InvalidArgumentException("rowIndex out of range for some observation.");
                    out[i] = data[i][rowIndex];
                }
                return out;
            }

            // Numerically stable mean (Kahan compensated)
            double StatisticsTools::Mean(array<double>^ data)
            {
                ValidateNonEmpty(data, "data");
                double sum = 0.0;
                double c = 0.0; // compensation
                for each (double v in data) {
                    double y = v - c;
                    double t = sum + y;
                    c = (t - sum) - y;
                    sum = t;
                }
                return sum / data->Length;
            }

            double StatisticsTools::Variance(array<double>^ data, bool sample)
            {
                ValidateNonEmpty(data, "data");
                int n = data->Length;
                if (n == 1) return 0.0;
                double mean = Mean(data);
                // use two-pass algorithm for better numeric stability
                double ssd = 0.0; // sum squared deviations
                for each (double v in data) {
                    double d = v - mean;
                    ssd += d * d;
                }
                return ssd / (sample ? (n - 1) : n);
            }

            double StatisticsTools::StdDev(array<double>^ data, bool sample)
            {
                double var = Variance(data, sample);
                return System::Math::Sqrt(System::Math::Max(0.0, var));
            }

            double StatisticsTools::WeightedMean(array<double>^ data, array<double>^ weights)
            {
                if (data == nullptr) throw gcnew ArgumentNullException("data");
                if (weights == nullptr) throw gcnew ArgumentNullException("weights");
                int n = data->Length;
                if (weights->Length != n) throw gcnew InvalidArgumentException("weights length must match data length.");
                double wsum = 0.0;
                double sum = 0.0;
                for (int i = 0; i < n; ++i) {
                    double w = weights[i];
                    if (System::Double::IsNaN(w) || System::Double::IsInfinity(w)) throw gcnew InvalidArgumentException("weights contain invalid value.");
                    sum += data[i] * w;
                    wsum += w;
                }
                if (Utils::IsZero(wsum)) throw gcnew InvalidArgumentException("sum of weights must be non-zero.");
                return sum / wsum;
            }

            double StatisticsTools::WeightedVariance(array<double>^ data, array<double>^ weights, bool sample)
            {
                if (data == nullptr) throw gcnew ArgumentNullException("data");
                if (weights == nullptr) throw gcnew ArgumentNullException("weights");
                int n = data->Length;
                if (weights->Length != n) throw gcnew InvalidArgumentException("weights length must match data length.");
                double wsum = 0.0;
                for (int i = 0; i < n; ++i) {
                    wsum += weights[i];
                }
                if (Utils::IsZero(wsum)) throw gcnew InvalidArgumentException("sum of weights must be non-zero.");
                double mu = WeightedMean(data, weights);

                double num = 0.0;
                for (int i = 0; i < n; ++i) {
                    double d = data[i] - mu;
                    num += weights[i] * d * d;
                }

                if (sample) {
                    // effective degrees of freedom approximation (Satterthwaite-like)
                    double w2 = 0.0;
                    for (int i = 0; i < n; ++i) w2 += weights[i] * weights[i];
                    double denom = wsum * wsum - w2;
                    if (Utils::IsZero(denom)) throw gcnew InvalidArgumentException("weights result in zero degrees of freedom for sample variance.");
                    return (wsum / denom) * num;
                }
                else {
                    return num / wsum;
                }
            }

            // -----------------------
            // Covariance and correlation
            // -----------------------
            double StatisticsTools::Covariance(array<double>^ x, array<double>^ y, bool sample)
            {
                ValidateNonEmpty(x, "x");
                if (y == nullptr) throw gcnew ArgumentNullException("y");
                if (y->Length != x->Length) throw gcnew InvalidArgumentException("x and y must have same length.");
                int n = x->Length;
                if (n == 1) return 0.0;
                double mx = Mean(x);
                double my = Mean(y);
                double s = 0.0;
                for (int i = 0; i < n; ++i) {
                    s += (x[i] - mx) * (y[i] - my);
                }
                return s / (sample ? (n - 1) : n);
            }

            array<array<double>^>^ StatisticsTools::CovarianceMatrix(array<array<double>^>^ data, bool sample)
            {
                if (data == nullptr) throw gcnew ArgumentNullException("data");
                int n = data->Length;
                if (n == 0) return gcnew array<array<double>^>(0);
                int d = data[0]->Length;
                for (int i = 0; i < n; ++i) {
                    if (data[i] == nullptr) throw gcnew InvalidArgumentException("each observation must be an array");
                    if (data[i]->Length != d) throw gcnew InvalidArgumentException("all observations must have the same dimension");
                }

                // compute means per dimension
                array<double>^ mean = gcnew array<double>(d);
                for (int j = 0; j < d; ++j) mean[j] = 0.0;
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < d; ++j) mean[j] += data[i][j];
                for (int j = 0; j < d; ++j) mean[j] /= (double)n;

                array<array<double>^>^ cov = gcnew array<array<double>^>(d);
                for (int i = 0; i < d; ++i) {
                    cov[i] = gcnew array<double>(d);
                    for (int j = 0; j < d; ++j) cov[i][j] = 0.0;
                }

                for (int t = 0; t < n; ++t) {
                    for (int i = 0; i < d; ++i) {
                        for (int j = 0; j < d; ++j) {
                            double xi = data[t][i] - mean[i];
                            double xj = data[t][j] - mean[j];
                            cov[i][j] += xi * xj;
                        }
                    }
                }

                double denom = sample ? (double)(n - 1) : (double)n;
                for (int i = 0; i < d; ++i)
                    for (int j = 0; j < d; ++j)
                        cov[i][j] /= denom;

                return cov;
            }

            double StatisticsTools::PearsonCorrelation(array<double>^ x, array<double>^ y)
            {
                ValidateNonEmpty(x, "x");
                if (y == nullptr) throw gcnew ArgumentNullException("y");
                if (y->Length != x->Length) throw gcnew InvalidArgumentException("x and y must have same length.");
                double cov = Covariance(x, y, true);
                double sx = StdDev(x, true);
                double sy = StdDev(y, true);
                if (Utils::IsZero(sx) || Utils::IsZero(sy)) throw gcnew InvalidArgumentException("standard deviation is zero for one of the series.");
                return cov / (sx * sy);
            }

            // -----------------------
            // Hypothesis tests and intervals
            // -----------------------
            double StatisticsTools::ZTestOneSample(array<double>^ sample, double mu0, double sigma)
            {
                ValidateNonEmpty(sample, "sample");
                if (sigma <= 0.0) throw gcnew InvalidArgumentException("sigma must be > 0.");
                int n = sample->Length;
                double m = Mean(sample);
                double z = (m - mu0) / (sigma / std::sqrt((double)n));
                return z;
            }

            double StatisticsTools::TTestOneSample(array<double>^ sample, double mu0)
            {
                ValidateNonEmpty(sample, "sample");
                int n = sample->Length;
                if (n < 2) throw gcnew InvalidArgumentException("at least two observations required for t-test.");
                double m = Mean(sample);
                double s = StdDev(sample, true);
                double se = s / std::sqrt((double)n);
                return (m - mu0) / se;
            }

            double StatisticsTools::TTestTwoSampleWelch(array<double>^ a, array<double>^ b)
            {
                ValidateNonEmpty(a, "a");
                ValidateNonEmpty(b, "b");
                double ma = Mean(a), mb = Mean(b);
                double va = Variance(a, true), vb = Variance(b, true);
                double na = a->Length, nb = b->Length;
                double se = std::sqrt(va / na + vb / nb);
                if (Utils::IsZero(se)) throw gcnew InvalidArgumentException("standard error is zero.");
                return (ma - mb) / se;
            }

            // Mean confidence interval (two-sided)
            KeyValuePair<double, double> StatisticsTools::MeanConfidenceInterval(array<double>^ sample, double confLevel, bool useT)
            {
                ValidateNonEmpty(sample, "sample");
                if (confLevel <= 0.0 || confLevel >= 1.0) throw gcnew InvalidArgumentException("confLevel must be in (0,1).");
                int n = sample->Length;
                double m = Mean(sample);
                double s = StdDev(sample, useT);
                double alpha = 1.0 - confLevel;
                double tail = 1.0 - alpha / 2.0;
                double quantile;
                if (useT) {
                    // For now approximate t quantile using normal if df large; a proper inverse-t implementation can be added
                    if (n - 1 >= 30) quantile = QuantileNormal(tail);
                    else quantile = QuantileNormal(tail); // approximation
                }
                else {
                    quantile = QuantileNormal(tail);
                }
                double se = s / std::sqrt((double)n);
                double lo = m - quantile * se;
                double hi = m + quantile * se;
                return KeyValuePair<double, double>(lo, hi);
            }

            // OLS simple linear regression
            Tuple<double, double, double, double, double>^ StatisticsTools::OLS_Simple(array<double>^ x, array<double>^ y)
            {
                ValidateNonEmpty(x, "x");
                if (y == nullptr) throw gcnew ArgumentNullException("y");
                if (y->Length != x->Length) throw gcnew InvalidArgumentException("x and y must have same length.");
                int n = x->Length;
                if (n < 2) throw gcnew InvalidArgumentException("at least two observations required for OLS.");

                double xm = Mean(x);
                double ym = Mean(y);
                double num = 0.0, den = 0.0;
                for (int i = 0; i < n; ++i) {
                    num += (x[i] - xm) * (y[i] - ym);
                    den += (x[i] - xm) * (x[i] - xm);
                }
                if (Utils::IsZero(den)) throw gcnew InvalidArgumentException("x has zero variance.");

                double beta = num / den;
                double alpha = ym - beta * xm;

                // estimate residual standard error
                double ssr = 0.0;
                for (int i = 0; i < n; ++i) {
                    double pred = alpha + beta * x[i];
                    double r = y[i] - pred;
                    ssr += r * r;
                }
                double sigma2 = ssr / (n - 2); // unbiased
                double stderr_beta = std::sqrt(sigma2 / den);
                double stderr_alpha = std::sqrt(sigma2 * (1.0 / n + xm * xm / den));

                // R^2
                double sst = 0.0;
                for (int i = 0; i < n; ++i) sst += (y[i] - ym) * (y[i] - ym);
                double r2 = 1.0 - ssr / sst;

                return gcnew Tuple<double, double, double, double, double>(alpha, beta, stderr_beta, stderr_alpha, r2);
            }

            // Bootstrap mean CI (percentile)
            KeyValuePair<double, double> StatisticsTools::BootstrapMeanCI(array<double>^ sample, double confLevel, int nBootstrap, unsigned int seed)
            {
                ValidateNonEmpty(sample, "sample");
                if (confLevel <= 0.0 || confLevel >= 1.0) throw gcnew InvalidArgumentException("confLevel must be in (0,1).");
                if (nBootstrap <= 0) throw gcnew InvalidArgumentException("nBootstrap must be positive.");

                int n = sample->Length;
                MersenneTwister mt(seed == 0 ? (unsigned int)Environment::TickCount : seed);
                List<double>^ boots = gcnew List<double>(nBootstrap);
                for (int b = 0; b < nBootstrap; ++b) {
                    // resample with replacement
                    double sum = 0.0;
                    for (int i = 0; i < n; ++i) {
                        int idx = mt.NextInt(n); // 0..n-1
                        sum += sample[idx];
                    }
                    boots->Add(sum / n);
                }
                // sort bootstrap means
                boots->Sort();
                double alpha = 1.0 - confLevel;
                int loIndex = (int)System::Math::Floor((alpha / 2.0) * nBootstrap);
                int hiIndex = (int)System::Math::Ceiling((1.0 - alpha / 2.0) * nBootstrap) - 1;
                loIndex = System::Math::Max(0, System::Math::Min(nBootstrap - 1, loIndex));
                hiIndex = System::Math::Max(0, System::Math::Min(nBootstrap - 1, hiIndex));
                return KeyValuePair<double, double>(boots[loIndex], boots[hiIndex]);
            }

            double StatisticsTools::SampleSizeForMean(double sigma, double marginError, double confLevel, bool useZ)
            {
                if (sigma <= 0.0) throw gcnew InvalidArgumentException("sigma must be > 0.");
                if (marginError <= 0.0) throw gcnew InvalidArgumentException("marginError must be > 0.");
                if (confLevel <= 0.0 || confLevel >= 1.0) throw gcnew InvalidArgumentException("confLevel must be in (0,1).");

                double alpha = 1.0 - confLevel;
                double tail = 1.0 - alpha / 2.0;
                double z = QuantileNormal(tail);
                // if using t, cannot compute exact without df; user likely expects z
                double q = z;
                double n = (q * sigma / marginError);
                n = n * n;
                return std::ceil(n);
            }

            // -----------------------
            // Quantile helpers (approximations)
            // -----------------------
            // Inverse normal CDF using probit approximation (Beasley-Springer-Moro / simpler)
            double StatisticsTools::QuantileNormal(double p)
            {
                if (p <= 0.0 || p >= 1.0) throw gcnew InvalidArgumentException("p must be in (0,1).");
                // use std::erf inverse via approximation: inverse erf available via std::erf? not standardized
                // We'll use approximation: quantile = sqrt(2) * erfinv(2p-1)
                // Implement erfinv using a rational approximation (Winitzki) for speed/robustness

                double a = 0.147; // constant for Winitzki approximation
                double x = 2.0 * p - 1.0;
                double sign = (x < 0) ? -1.0 : 1.0;
                double lnTerm = std::log(1.0 - x * x);
                double first = 2.0 / (M_E * a) + lnTerm / 2.0;
                double second = lnTerm / a;
                double inside = std::sqrt(first * first - second);
                double erfinv = sign * std::sqrt(inside - first);
                return std::sqrt(2.0) * erfinv;
            }

            // Placeholder: QuantileT - for simplicity we approximate t quantile by normal for df>30,
            // otherwise we fallback to normal approximation (user should be warned if strict accuracy is required).
            double StatisticsTools::QuantileT(double p, int df)
            {
                // for production-quality code replace with precise inverse-t implementation
                if (df <= 0) throw gcnew InvalidArgumentException("df must be positive.");
                return QuantileNormal(p);
            }

        } // namespace Statistics
    } // namespace Math
} // namespace Agrine

#pragma once

// Statistics.h - Statistical utilities and tests (univariate & multivariate)
// Part of Agrine.Math (C++/CLI, .NET Framework 4.6.2)

#include "Constants.h"
#include "Utils.h"
#include "Exceptions.h"
#include "Randgen.h"   // use Randgen::MersenneTwister for bootstrap
#include <cmath>

using namespace System;
using namespace System::Collections::Generic;
using namespace Agrine::Math::Probability;

namespace Agrine {
    namespace Math {
        namespace Statistics {

            /// <summary>
            /// Statistics: collection of statistical utilities (mean, variance, covariance, tests, bootstrap, OLS).
            /// All methods are static.
            /// </summary>
            public ref class StatisticsTools abstract sealed
            {
            public:
                // -------- Basic descriptive (univariate) --------
                static double Mean(array<double>^ data);
                // sampleVariance: uses (n-1) denominator if sample==true, otherwise population (n)
                static double Variance(array<double>^ data, bool sample);
                static double StdDev(array<double>^ data, bool sample);

                // weighted mean/variance (weights must sum > 0)
                static double WeightedMean(array<double>^ data, array<double>^ weights);
                static double WeightedVariance(array<double>^ data, array<double>^ weights, bool sample);

                // -------- Multivariate summary --------
                // data: array of observations, each observation is array<double> of length d
                // returns covariance matrix as array of rows (d arrays)
                static array<array<double>^>^ CovarianceMatrix(array<array<double>^>^ data, bool sample);

                // covariance between two arrays (univariate)
                static double Covariance(array<double>^ x, array<double>^ y, bool sample);

                // Pearson correlation coefficient
                static double PearsonCorrelation(array<double>^ x, array<double>^ y);

                // -------- Hypothesis tests & intervals --------
                // z-test (one-sample) - assume known sigma
                // returns z-statistic
                static double ZTestOneSample(array<double>^ sample, double mu0, double sigma);

                // t-test (one-sample) - returns t-statistic (df = n-1)
                static double TTestOneSample(array<double>^ sample, double mu0);

                // two-sample unpaired t-test (Welch) - returns t-statistic
                static double TTestTwoSampleWelch(array<double>^ a, array<double>^ b);

                // compute (approx) two-sided confidence interval for mean
                // if useT==true uses t quantile with df=n-1, otherwise z (normal)
                static KeyValuePair<double, double> MeanConfidenceInterval(array<double>^ sample, double confLevel, bool useT);

                // -------- Regression (simple OLS) --------
                // simple linear regression y = alpha + beta*x
                // returns tuple: (alpha, beta, stderr_beta, stderr_alpha, R2)
                static Tuple<double, double, double, double, double>^ OLS_Simple(array<double>^ x, array<double>^ y);

                // -------- Bootstrap --------
                // bootstrapMeanCI: bootstrap sample means, returns two-sided CI
                // nBootstrap default 1000
                static KeyValuePair<double, double> BootstrapMeanCI(array<double>^ sample, double confLevel, int nBootstrap, unsigned int seed);

                // -------- Utility --------
                static double SampleSizeForMean(double sigma, double marginError, double confLevel, bool useZ);

            private:
                // helpers
                static void ValidateNonEmpty(array<double>^ a, String^ name);
                static array<double>^ ToDoubleArrayFromMatrixRow(array<array<double>^>^ data, int rowIndex);
                static double QuantileNormal(double p); // inverse CDF normal (approx) using std::erf inverse via approximation
                static double QuantileT(double p, int df); // placeholder: use approximation via inverse student's t (simple Newton on CDF could be added)
            };
        }
    }
}

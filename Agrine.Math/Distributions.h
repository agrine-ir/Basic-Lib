#pragma once

// Distributions.h - Probability distributions (PDF, CDF, sampling)
// Part of Agrine.Math (C++/CLI, .NET Framework 4.6.2)

#include "Constants.h"
#include "Utils.h"
#include "Exceptions.h"

using namespace System;

namespace Agrine {
    namespace Math {
        namespace Probability {

            /// <summary>
            /// Distributions - provides PDF/CDF/Sampling for common probability distributions.
            /// </summary>
            public ref class Distributions abstract sealed
            {
            public:
                // -------- Uniform(a,b) --------
                static double UniformPDF(double x, double a, double b);
                static double UniformCDF(double x, double a, double b);
                static double UniformSample(double a, double b);

                // -------- Normal(mu, sigma) --------
                static double NormalPDF(double x, double mu, double sigma);
                static double NormalCDF(double x, double mu, double sigma);
                static double NormalSample(double mu, double sigma);

                // -------- Binomial(n,p) --------
                static double BinomialPMF(int k, int n, double p);
                static double BinomialCDF(int k, int n, double p);
                static int    BinomialSample(int n, double p);

                // -------- Poisson(lambda) --------
                static double PoissonPMF(int k, double lambda);
                static double PoissonCDF(int k, double lambda);
                static int    PoissonSample(double lambda);

            private:
                // random generator (lazy initialized)
                static Random^ rng;
                static Random^ GetRng();
            };

        } // namespace Probability
    } // namespace Math
} // namespace Agrine

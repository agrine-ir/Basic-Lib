#pragma once
#include "MathCore.h"
#include <vector>

namespace Agrine
{
    namespace Math
    {
        namespace Statistics
        {
            // ================================
            // Statistics Module
            // ================================
            // This module provides basic and advanced statistical operations
            // for data analysis, including descriptive statistics,
            // probability distributions, and statistical tests.

            // ----------------------------
            // Descriptive Statistics
            // ----------------------------
            AGRINE_API double Mean(const std::vector<double>& data);
            AGRINE_API double Median(std::vector<double> data);
            AGRINE_API double Variance(const std::vector<double>& data, bool sample = true);
            AGRINE_API double StandardDeviation(const std::vector<double>& data, bool sample = true);
            AGRINE_API double Skewness(const std::vector<double>& data);
            AGRINE_API double Kurtosis(const std::vector<double>& data);

            // ----------------------------
            // Probability Distributions
            // ----------------------------
            AGRINE_API double NormalPDF(double x, double mean, double stddev);
            AGRINE_API double NormalCDF(double x, double mean, double stddev);
            AGRINE_API double UniformPDF(double x, double a, double b);
            AGRINE_API double UniformCDF(double x, double a, double b);

            // ----------------------------
            // Correlation and Regression
            // ----------------------------
            AGRINE_API double PearsonCorrelation(const std::vector<double>& x, const std::vector<double>& y);
            AGRINE_API std::pair<double, double> LinearRegression(const std::vector<double>& x, const std::vector<double>& y);

            // ----------------------------
            // Hypothesis Testing
            // ----------------------------
            AGRINE_API double TTest(const std::vector<double>& sample1, const std::vector<double>& sample2, bool equalVariance = true);

            // ----------------------------
            // Utility Functions
            // ----------------------------
            AGRINE_API double Sum(const std::vector<double>& data);
        }
    }
}

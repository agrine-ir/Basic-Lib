#include "Statistics.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <stdexcept>

using namespace Agrine::Math;
using namespace Agrine::Math::Constants;

// ----------------------------
// Descriptive Statistics
// ----------------------------
double Statistics::Mean(const std::vector<double>& data)
{
    if (data.empty()) throw std::invalid_argument("Data cannot be empty.");
    return std::accumulate(data.begin(), data.end(), 0.0) / data.size();
}

double Statistics::Median(std::vector<double> data)
{
    if (data.empty()) throw std::invalid_argument("Data cannot be empty.");
    std::sort(data.begin(), data.end());
    size_t n = data.size();
    if (n % 2 == 0) return (data[n / 2 - 1] + data[n / 2]) / 2.0;
    return data[n / 2];
}

double Statistics::Variance(const std::vector<double>& data, bool sample)
{
    if (data.size() < 2) throw std::invalid_argument("Data must contain at least two elements.");
    double mean = Mean(data);
    double sumSq = 0;
    for (double x : data) sumSq += (x - mean) * (x - mean);
    return sumSq / (data.size() - (sample ? 1 : 0));
}

double Statistics::StandardDeviation(const std::vector<double>& data, bool sample)
{
    return std::sqrt(Variance(data, sample));
}

double Statistics::Skewness(const std::vector<double>& data)
{
    double mean = Mean(data);
    double stddev = StandardDeviation(data);
    double sumCubed = 0;
    for (double x : data) sumCubed += std::pow(x - mean, 3);
    return (sumCubed / data.size()) / std::pow(stddev, 3);
}

double Statistics::Kurtosis(const std::vector<double>& data)
{
    double mean = Mean(data);
    double stddev = StandardDeviation(data);
    double sumQuad = 0;
    for (double x : data) sumQuad += std::pow(x - mean, 4);
    return (sumQuad / data.size()) / std::pow(stddev, 4) - 3.0; // Excess kurtosis
}

// ----------------------------
// Probability Distributions
// ----------------------------
double Statistics::NormalPDF(double x, double mean, double stddev)
{
    double var = stddev * stddev;
    return (1.0 / (stddev * std::sqrt(2 * PI))) * std::exp(-0.5 * std::pow((x - mean) / stddev, 2));
}

double Statistics::NormalCDF(double x, double mean, double stddev)
{
    return 0.5 * (1 + std::erf((x - mean) / (stddev * std::sqrt(2))));
}

double Statistics::UniformPDF(double x, double a, double b)
{
    if (a >= b) throw std::invalid_argument("Invalid range for Uniform distribution.");
    return (x >= a && x <= b) ? 1.0 / (b - a) : 0.0;
}

double Statistics::UniformCDF(double x, double a, double b)
{
    if (a >= b) throw std::invalid_argument("Invalid range for Uniform distribution.");
    if (x < a) return 0.0;
    if (x > b) return 1.0;
    return (x - a) / (b - a);
}

// ----------------------------
// Correlation and Regression
// ----------------------------
double Statistics::PearsonCorrelation(const std::vector<double>& x, const std::vector<double>& y)
{
    if (x.size() != y.size() || x.empty()) throw std::invalid_argument("Vectors must have the same size and not be empty.");
    double meanX = Mean(x), meanY = Mean(y);
    double num = 0, denX = 0, denY = 0;
    for (size_t i = 0; i < x.size(); i++)
    {
        num += (x[i] - meanX) * (y[i] - meanY);
        denX += (x[i] - meanX) * (x[i] - meanX);
        denY += (y[i] - meanY) * (y[i] - meanY);
    }
    return num / std::sqrt(denX * denY);
}

std::pair<double, double> Statistics::LinearRegression(const std::vector<double>& x, const std::vector<double>& y)
{
    if (x.size() != y.size() || x.empty()) throw std::invalid_argument("Vectors must have the same size and not be empty.");
    double meanX = Mean(x), meanY = Mean(y);
    double num = 0, den = 0;
    for (size_t i = 0; i < x.size(); i++)
    {
        num += (x[i] - meanX) * (y[i] - meanY);
        den += (x[i] - meanX) * (x[i] - meanX);
    }
    double slope = num / den;
    double intercept = meanY - slope * meanX;
    return { slope, intercept };
}

// ----------------------------
// Hypothesis Testing
// ----------------------------
double Statistics::TTest(const std::vector<double>& sample1, const std::vector<double>& sample2, bool equalVariance)
{
    double mean1 = Mean(sample1), mean2 = Mean(sample2);
    double var1 = Variance(sample1), var2 = Variance(sample2);
    size_t n1 = sample1.size(), n2 = sample2.size();

    if (equalVariance)
    {
        double pooledVar = ((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2);
        return (mean1 - mean2) / std::sqrt(pooledVar * (1.0 / n1 + 1.0 / n2));
    }
    else
    {
        return (mean1 - mean2) / std::sqrt(var1 / n1 + var2 / n2);
    }
}

// ----------------------------
// Utility Functions
// ----------------------------
double Statistics::Sum(const std::vector<double>& data)
{
    return std::accumulate(data.begin(), data.end(), 0.0);
}

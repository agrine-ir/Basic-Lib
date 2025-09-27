#include "Utils.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <random>
#include <stdexcept>

using namespace Agrine::Math;
using namespace Agrine::Math::Constants;

//
// ----------------------------
// Basic Math Helpers
// ----------------------------
int Utils::GCD(int a, int b)
{
    while (b != 0)
    {
        int temp = b;
        b = a % b;
        a = temp;
    }
    return std::abs(a);
}

int Utils::LCM(int a, int b)
{
    if (a == 0 || b == 0) return 0;
    return std::abs(a * b) / GCD(a, b);
}

long long Utils::Factorial(int n)
{
    if (n < 0) throw std::invalid_argument("Factorial is not defined for negative numbers.");
    long long result = 1;
    for (int i = 2; i <= n; i++) result *= i;
    return result;
}

double Utils::Power(double base, int exponent)
{
    double result = 1.0;
    if (exponent >= 0)
    {
        for (int i = 0; i < exponent; i++) result *= base;
    }
    else
    {
        for (int i = 0; i < -exponent; i++) result /= base;
    }
    return result;
}

//
// ----------------------------
// Number Utilities
// ----------------------------
bool Utils::IsPrime(int n)
{
    if (n <= 1) return false;
    if (n <= 3) return true;
    if (n % 2 == 0 || n % 3 == 0) return false;
    for (int i = 5; i * i <= n; i += 6)
    {
        if (n % i == 0 || n % (i + 2) == 0) return false;
    }
    return true;
}

std::vector<int> Utils::PrimeFactors(int n)
{
    std::vector<int> factors;
    if (n <= 1) return factors;
    while (n % 2 == 0)
    {
        factors.push_back(2);
        n /= 2;
    }
    for (int i = 3; i * i <= n; i += 2)
    {
        while (n % i == 0)
        {
            factors.push_back(i);
            n /= i;
        }
    }
    if (n > 2) factors.push_back(n);
    return factors;
}

double Utils::Clamp(double value, double minVal, double maxVal)
{
    return std::max(minVal, std::min(value, maxVal));
}

bool Utils::AlmostEqual(double a, double b, double epsilon)
{
    return std::fabs(a - b) <= epsilon;
}

//
// ----------------------------
// Random Utilities
// ----------------------------
double Utils::RandomDouble(double minVal, double maxVal)
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(minVal, maxVal);
    return dist(gen);
}

int Utils::RandomInt(int minVal, int maxVal)
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(minVal, maxVal);
    return dist(gen);
}

double Utils::RandomChoice(const std::vector<double>& data)
{
    if (data.empty()) throw std::invalid_argument("Cannot choose from empty vector.");
    int idx = RandomInt(0, static_cast<int>(data.size()) - 1);
    return data[idx];
}

//
// ----------------------------
// Vector Utilities
// ----------------------------
std::vector<double> Utils::Normalize(const std::vector<double>& data)
{
    if (data.empty()) throw std::invalid_argument("Data cannot be empty.");
    std::pair<double, double> mm = MinMax(data);
    double minVal = mm.first;
    double maxVal = mm.second;
    std::vector<double> result;
    result.reserve(data.size());
    for (double x : data)
    {
        result.push_back((x - minVal) / (maxVal - minVal));
    }
    return result;
}

double Utils::DotProduct(const std::vector<double>& a, const std::vector<double>& b)
{
    if (a.size() != b.size()) throw std::invalid_argument("Vectors must have the same size.");
    return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
}

std::pair<double, double> Utils::MinMax(const std::vector<double>& data)
{
    if (data.empty()) throw std::invalid_argument("Data cannot be empty.");
    auto mm = std::minmax_element(data.begin(), data.end());
    auto minIt = mm.first;
    auto maxIt = mm.second;
    return { *minIt, *maxIt };
}

//
// ----------------------------
// Matrix Utilities
// ----------------------------
std::vector<std::vector<double>> Utils::IdentityMatrix(size_t n)
{
    std::vector<std::vector<double>> mat(n, std::vector<double>(n, 0.0));
    for (size_t i = 0; i < n; i++) mat[i][i] = 1.0;
    return mat;
}

std::vector<std::vector<double>> Utils::ZeroMatrix(size_t rows, size_t cols)
{
    return std::vector<std::vector<double>>(rows, std::vector<double>(cols, 0.0));
}

std::vector<std::vector<double>> Utils::RandomMatrix(size_t rows, size_t cols,
    double minVal, double maxVal)
{
    std::vector<std::vector<double>> mat(rows, std::vector<double>(cols));
    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < cols; j++)
        {
            mat[i][j] = RandomDouble(minVal, maxVal);
        }
    }
    return mat;
}

std::vector<std::vector<double>> Utils::Transpose(const std::vector<std::vector<double>>& mat)
{
    if (mat.empty()) return {};
    size_t rows = mat.size();
    size_t cols = mat[0].size();
    std::vector<std::vector<double>> result(cols, std::vector<double>(rows));
    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < cols; j++)
        {
            result[j][i] = mat[i][j];
        }
    }
    return result;
}

bool Utils::SameDimensions(const std::vector<std::vector<double>>& a,
    const std::vector<std::vector<double>>& b)
{
    if (a.size() != b.size()) return false;
    if (!a.empty() && !b.empty() && a[0].size() != b[0].size()) return false;
    return true;
}

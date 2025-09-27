#include "Linalg.h"
#include <cmath>
#include <stdexcept>

using namespace Agrine::Math;

// ----------------------------
// Vector Operations
// ----------------------------
double Linalg::VectorNorm(const std::vector<double>& v)
{
    double sum = 0;
    for (auto val : v) sum += val * val;
    return std::sqrt(sum);
}

std::vector<double> Linalg::VectorNormalize(const std::vector<double>& v)
{
    double norm = VectorNorm(v);
    if (norm == 0) return std::vector<double>(v.size(), 0.0);
    std::vector<double> result(v.size());
    for (size_t i = 0; i < v.size(); i++) result[i] = v[i] / norm;
    return result;
}

double Linalg::VectorDot(const std::vector<double>& v1, const std::vector<double>& v2)
{
    if (v1.size() != v2.size()) throw std::invalid_argument("Vector sizes do not match.");
    double sum = 0;
    for (size_t i = 0; i < v1.size(); i++) sum += v1[i] * v2[i];
    return sum;
}

std::vector<double> Linalg::VectorCross(const std::vector<double>& v1, const std::vector<double>& v2)
{
    if (v1.size() != 3 || v2.size() != 3) throw std::invalid_argument("Cross product is defined only for 3D vectors.");
    return { v1[1] * v2[2] - v1[2] * v2[1],
             v1[2] * v2[0] - v1[0] * v2[2],
             v1[0] * v2[1] - v1[1] * v2[0] };
}

// ----------------------------
// Matrix Operations
// ----------------------------
std::vector<std::vector<double>> Linalg::MatrixTranspose(const std::vector<std::vector<double>>& mat)
{
    if (mat.empty()) return {};
    size_t rows = mat.size(), cols = mat[0].size();
    std::vector<std::vector<double>> result(cols, std::vector<double>(rows));
    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < cols; j++)
            result[j][i] = mat[i][j];
    return result;
}

std::vector<std::vector<double>> Linalg::MatrixMultiply(const std::vector<std::vector<double>>& A,
    const std::vector<std::vector<double>>& B)
{
    if (A.empty() || B.empty()) return {};
    if (A[0].size() != B.size()) throw std::invalid_argument("Matrix dimensions do not match for multiplication.");
    size_t rows = A.size(), cols = B[0].size(), inner = B.size();
    std::vector<std::vector<double>> result(rows, std::vector<double>(cols, 0.0));
    for (size_t i = 0; i < rows; i++)
        for (size_t j = 0; j < cols; j++)
            for (size_t k = 0; k < inner; k++)
                result[i][j] += A[i][k] * B[k][j];
    return result;
}

double Linalg::MatrixDeterminant(const std::vector<std::vector<double>>& mat)
{
    // This will only implement for 2x2 and 3x3 matrices for now
    if (mat.size() != mat[0].size()) throw std::invalid_argument("Matrix must be square.");
    if (mat.size() == 2)
        return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
    if (mat.size() == 3)
        return mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1])
        - mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0])
        + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
    throw std::runtime_error("Determinant calculation not implemented for this size.");
}

std::vector<std::vector<double>> Linalg::MatrixInverse(const std::vector<std::vector<double>>& mat, bool& success)
{
    success = false;
    // Only simple case: 2x2 matrix inversion
    if (mat.size() == 2 && mat[0].size() == 2)
    {
        double det = MatrixDeterminant(mat);
        if (det == 0) return {};
        success = true;
        double invDet = 1.0 / det;
        return {
            { mat[1][1] * invDet, -mat[0][1] * invDet },
            { -mat[1][0] * invDet, mat[0][0] * invDet }
        };
    }
    return {};
}

// ----------------------------
// Decompositions
// ----------------------------
bool Linalg::LUDecomposition(const std::vector<std::vector<double>>& mat,
    std::vector<std::vector<double>>& L,
    std::vector<std::vector<double>>& U)
{
    size_t n = mat.size();
    if (n != mat[0].size()) return false;
    L = std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));
    U = std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));

    for (size_t i = 0; i < n; i++)
    {
        for (size_t k = i; k < n; k++)
        {
            double sum = 0.0;
            for (size_t j = 0; j < i; j++)
                sum += L[i][j] * U[j][k];
            U[i][k] = mat[i][k] - sum;
        }
        for (size_t k = i; k < n; k++)
        {
            if (i == k) L[i][i] = 1.0;
            else
            {
                double sum = 0.0;
                for (size_t j = 0; j < i; j++)
                    sum += L[k][j] * U[j][i];
                L[k][i] = (mat[k][i] - sum) / U[i][i];
            }
        }
    }
    return true;
}

// Placeholder for QR decomposition
bool Linalg::QRDecomposition(const std::vector<std::vector<double>>& mat,
    std::vector<std::vector<double>>& Q,
    std::vector<std::vector<double>>& R)
{
    return false; // Implement if needed
}

// Placeholder for eigen calculations
std::vector<double> Linalg::Eigenvalues(const std::vector<std::vector<double>>& mat)
{
    return {}; // Implement if needed
}

std::vector<std::vector<double>> Linalg::Eigenvectors(const std::vector<std::vector<double>>& mat)
{
    return {}; // Implement if needed
}

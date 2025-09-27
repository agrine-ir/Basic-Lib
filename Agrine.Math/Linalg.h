#pragma once
#include "MathCore.h"
#include <vector>
#include <array>

namespace Agrine
{
    namespace Math
    {
        namespace Linalg
        {
            // ================================
            // Linear Algebra Module
            // ================================
            // This module provides advanced linear algebra operations,
            // including matrix operations, vector operations, and decompositions.

            // ----------------------------
            // Vector Operations
            // ----------------------------
            AGRINE_API double VectorNorm(const std::vector<double>& v);
            AGRINE_API std::vector<double> VectorNormalize(const std::vector<double>& v);
            AGRINE_API double VectorDot(const std::vector<double>& v1, const std::vector<double>& v2);
            AGRINE_API std::vector<double> VectorCross(const std::vector<double>& v1, const std::vector<double>& v2);

            // ----------------------------
            // Matrix Operations
            // ----------------------------
            AGRINE_API std::vector<std::vector<double>> MatrixTranspose(const std::vector<std::vector<double>>& mat);
            AGRINE_API std::vector<std::vector<double>> MatrixMultiply(const std::vector<std::vector<double>>& A,
                const std::vector<std::vector<double>>& B);
            AGRINE_API std::vector<std::vector<double>> MatrixInverse(const std::vector<std::vector<double>>& mat, bool& success);

            AGRINE_API double MatrixDeterminant(const std::vector<std::vector<double>>& mat);

            // ----------------------------
            // Decompositions
            // ----------------------------
            AGRINE_API bool LUDecomposition(const std::vector<std::vector<double>>& mat,
                std::vector<std::vector<double>>& L,
                std::vector<std::vector<double>>& U);

            AGRINE_API bool QRDecomposition(const std::vector<std::vector<double>>& mat,
                std::vector<std::vector<double>>& Q,
                std::vector<std::vector<double>>& R);

            // ----------------------------
            // Eigenvalues and Eigenvectors
            // ----------------------------
            AGRINE_API std::vector<double> Eigenvalues(const std::vector<std::vector<double>>& mat);
            AGRINE_API std::vector<std::vector<double>> Eigenvectors(const std::vector<std::vector<double>>& mat);
        }
    }
}

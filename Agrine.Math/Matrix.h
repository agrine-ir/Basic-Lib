#pragma once

// Matrix.h - Algebra module
// Part of Agrine.Math (C++/CLI, .NET Framework 4.6.2)

#include "Vector.h"

using namespace System;
using namespace Agrine::Math::Core;

namespace Agrine {
    namespace Math {
        namespace Algebra {

            /// <summary>
            /// Represents a mathematical matrix (m x n).
            /// Provides basic and advanced linear algebra operations.
            /// </summary>
            public ref class Matrix
            {
            private:
                int rows;
                int cols;
                array<double, 2>^ data; // 2D managed array

            public:
                // ===== Properties =====
                property int Rows { int get(); }
                property int Cols { int get(); }

                property double default[int, int] {
                    double get(int r, int c);
                    void set(int r, int c, double value);
                    }

                    // ===== Constructors =====
                Matrix(int rows, int cols);
                Matrix(array<double, 2>^ values);

                // ===== Utility =====
                Matrix^ Clone();
                static Matrix^ Identity(int n);

                // ===== Basic operations =====
                Matrix^ Add(Matrix^ other);
                Matrix^ Subtract(Matrix^ other);
                Matrix^ Multiply(double scalar);
                Vector^ Multiply(Vector^ v);
                Matrix^ Multiply(Matrix^ other);

                Matrix^ Transpose();

                // ===== Advanced operations =====
                double Determinant();     // Only for square matrices
                Matrix^ Inverse();        // Only for square matrices
                int Rank();

                static bool AreEqual(Matrix^ a, Matrix^ b, double tol = Constants::Tolerance);

                // ===== ToString override =====
                virtual String^ ToString() override;
            };
        }
    }
}

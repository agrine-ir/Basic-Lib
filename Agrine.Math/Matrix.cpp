#include "Matrix.h"

using namespace System;
using namespace Agrine::Math::Core;

namespace Agrine {
    namespace Math {
        namespace Algebra {

            // ===== Constructors =====
            Matrix::Matrix(int r, int c)
            {
                if (r <= 0 || c <= 0)
                    throw gcnew ArgumentException("Matrix dimensions must be positive.");

                this->rows = r;
                this->cols = c;
                data = gcnew array<double, 2>(r, c);
            }

            Matrix::Matrix(array<double, 2>^ values)
            {
                if (values == nullptr || values->Length == 0)
                    throw gcnew ArgumentException("Values array cannot be null or empty.");

                this->rows = values->GetLength(0);
                this->cols = values->GetLength(1);
                data = gcnew array<double, 2>(rows, cols);

                for (int i = 0; i < rows; i++)
                    for (int j = 0; j < cols; j++)
                        data[i, j] = values[i, j];
            }

            // ===== Properties =====
            int Matrix::Rows::get() { return rows; }
            int Matrix::Cols::get() { return cols; }

            double Matrix::default::get(int r, int c)
            {
                if (r < 0 || r >= rows || c < 0 || c >= cols)
                    throw gcnew IndexOutOfRangeException("Matrix index out of range.");
                return data[r, c];
            }

            void Matrix::default::set(int r, int c, double value)
            {
                if (r < 0 || r >= rows || c < 0 || c >= cols)
                    throw gcnew IndexOutOfRangeException("Matrix index out of range.");
                data[r, c] = value;
            }

            // ===== Utility =====
            Matrix^ Matrix::Clone()
            {
                return gcnew Matrix(data);
            }

            Matrix^ Matrix::Identity(int n)
            {
                if (n <= 0)
                    throw gcnew ArgumentException("Size must be positive.");

                Matrix^ I = gcnew Matrix(n, n);
                for (int i = 0; i < n; i++)
                    I[i, i] = 1.0;
                return I;
            }

            // ===== Basic operations =====
            Matrix^ Matrix::Add(Matrix^ other)
            {
                if (this->rows != other->rows || this->cols != other->cols)
                    throw gcnew DimensionMismatchException("Matrix dimensions must match for addition.");

                Matrix^ result = gcnew Matrix(rows, cols);
                for (int i = 0; i < rows; i++)
                    for (int j = 0; j < cols; j++)
                        result[i, j] = this->data[i, j] + other[i, j];
                return result;
            }

            Matrix^ Matrix::Subtract(Matrix^ other)
            {
                if (this->rows != other->rows || this->cols != other->cols)
                    throw gcnew DimensionMismatchException("Matrix dimensions must match for subtraction.");

                Matrix^ result = gcnew Matrix(rows, cols);
                for (int i = 0; i < rows; i++)
                    for (int j = 0; j < cols; j++)
                        result[i, j] = this->data[i, j] - other[i, j];
                return result;
            }

            Matrix^ Matrix::Multiply(double scalar)
            {
                Matrix^ result = gcnew Matrix(rows, cols);
                for (int i = 0; i < rows; i++)
                    for (int j = 0; j < cols; j++)
                        result[i, j] = this->data[i, j] * scalar;
                return result;
            }

            Vector^ Matrix::Multiply(Vector^ v)
            {
                if (this->cols != v->Dimension)
                    throw gcnew DimensionMismatchException("Matrix columns must match vector dimension.");

                array<double>^ result = gcnew array<double>(rows);
                for (int i = 0; i < rows; i++)
                {
                    double sum = 0.0;
                    for (int j = 0; j < cols; j++)
                        sum += this->data[i, j] * v[j];
                    result[i] = sum;
                }
                return gcnew Vector(result);
            }

            Matrix^ Matrix::Multiply(Matrix^ other)
            {
                if (this->cols != other->rows)
                    throw gcnew DimensionMismatchException("Matrix dimensions must match for multiplication.");

                Matrix^ result = gcnew Matrix(this->rows, other->cols);
                for (int i = 0; i < this->rows; i++)
                {
                    for (int j = 0; j < other->cols; j++)
                    {
                        double sum = 0.0;
                        for (int k = 0; k < this->cols; k++)
                            sum += this->data[i, k] * other[k, j];
                        result[i, j] = sum;
                    }
                }
                return result;
            }

            Matrix^ Matrix::Transpose()
            {
                Matrix^ result = gcnew Matrix(cols, rows);
                for (int i = 0; i < rows; i++)
                    for (int j = 0; j < cols; j++)
                        result[j, i] = this->data[i, j];
                return result;
            }

            // ===== Advanced operations =====

            // Determinant using Gaussian elimination (O(n^3))
            double Matrix::Determinant()
            {
                if (rows != cols)
                    throw gcnew InvalidArgumentException("Determinant is only defined for square matrices.");

                int n = rows;
                array<double, 2>^ temp = gcnew array<double, 2>(n, n);
                Array::Copy(data, temp, data->Length);

                double det = 1.0;
                for (int i = 0; i < n; i++)
                {
                    // Pivoting
                    int pivot = i;
                    for (int j = i + 1; j < n; j++)
                        if (System::Math::Abs(temp[j, i]) > System::Math::Abs(temp[pivot, i]))
                            pivot = j;

                    if (Utils::IsZero(temp[pivot, i]))
                        return 0.0;

                    // Swap rows if needed
                    if (pivot != i)
                    {
                        for (int k = 0; k < n; k++)
                        {
                            double t = temp[i, k];
                            temp[i, k] = temp[pivot, k];
                            temp[pivot, k] = t;
                        }
                        det = -det;
                    }

                    det *= temp[i, i];
                    for (int j = i + 1; j < n; j++)
                        temp[i, j] /= temp[i, i];

                    for (int j = i + 1; j < n; j++)
                    {
                        for (int k = i + 1; k < n; k++)
                            temp[j, k] -= temp[j, i] * temp[i, k];
                    }
                }
                return det;
            }

            // Inverse using Gauss-Jordan elimination
            Matrix^ Matrix::Inverse()
            {
                if (rows != cols)
                    throw gcnew InvalidArgumentException("Inverse is only defined for square matrices.");

                int n = rows;
                Matrix^ A = this->Clone();
                Matrix^ I = Identity(n);

                for (int i = 0; i < n; i++)
                {
                    // Find pivot
                    double pivot = A[i, i];
                    if (Utils::IsZero(pivot))
                        throw gcnew DivideByZeroException("Matrix is singular and cannot be inverted.");

                    // Scale row
                    for (int j = 0; j < n; j++)
                    {
                        A[i, j] /= pivot;
                        I[i, j] /= pivot;
                    }

                    // Eliminate other rows
                    for (int k = 0; k < n; k++)
                    {
                        if (k == i) continue;
                        double factor = A[k, i];
                        for (int j = 0; j < n; j++)
                        {
                            A[k, j] -= factor * A[i, j];
                            I[k, j] -= factor * I[i, j];
                        }
                    }
                }
                return I;
            }

            // Rank using row echelon form
            int Matrix::Rank()
            {
                array<double, 2>^ temp = gcnew array<double, 2>(rows, cols);
                Array::Copy(data, temp, data->Length);

                int rank = 0;
                int row = 0;

                for (int col = 0; col < cols && row < rows; col++)
                {
                    // Find pivot
                    int pivot = row;
                    for (int i = row + 1; i < rows; i++)
                        if (System::Math::Abs(temp[i, col]) > System::Math::Abs(temp[pivot, col]))
                            pivot = i;

                    if (Utils::IsZero(temp[pivot, col]))
                        continue;

                    // Swap
                    for (int j = 0; j < cols; j++)
                    {
                        double t = temp[row, j];
                        temp[row, j] = temp[pivot, j];
                        temp[pivot, j] = t;
                    }

                    // Normalize pivot row
                    double pv = temp[row, col];
                    for (int j = 0; j < cols; j++)
                        temp[row, j] /= pv;

                    // Eliminate
                    for (int i = 0; i < rows; i++)
                    {
                        if (i == row) continue;
                        double factor = temp[i, col];
                        for (int j = 0; j < cols; j++)
                            temp[i, j] -= factor * temp[row, j];
                    }

                    row++;
                    rank++;
                }

                return rank;
            }

            bool Matrix::AreEqual(Matrix^ a, Matrix^ b, double tol)
            {
                if (a->Rows != b->Rows || a->Cols != b->Cols)
                    return false;

                for (int i = 0; i < a->Rows; i++)
                    for (int j = 0; j < a->Cols; j++)
                        if (!Utils::AreEqual(a[i, j], b[i, j], tol))
                            return false;

                return true;
            }

            // ===== ToString =====
            String^ Matrix::ToString()
            {
                array<String^>^ rowsStr = gcnew array<String^>(rows);
                for (int i = 0; i < rows; i++)
                {
                    array<String^>^ elems = gcnew array<String^>(cols);
                    for (int j = 0; j < cols; j++)
                        elems[j] = data[i, j].ToString("G5");
                    rowsStr[i] = "[" + String::Join(", ", elems) + "]";
                }
                return "Matrix(" + String::Join("; ", rowsStr) + ")";
            }

        } // namespace Algebra
    } // namespace Math
} // namespace Agrine

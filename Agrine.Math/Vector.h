#pragma once

// Vector.h - Algebra module
// Part of Agrine.Math (C++/CLI, .NET Framework 4.6.2)


#include "Constants.h"
#include "Utils.h"
#include "Exceptions.h"

using namespace System;
using namespace System::Collections::Generic;

namespace Agrine {
    namespace Math {
        namespace Algebra {

            /// <summary>
            /// Represents a mathematical vector in n-dimensional space.
            /// Provides basic and advanced vector operations.
            /// </summary>
            public ref class Vector
            {
            private:
                array<double>^ data;

            public:
                /// <summary>
                /// Gets the dimension of the vector.
                /// </summary>
                property int Dimension {
                    int get() { return data->Length; }
                }

                /// <summary>
                /// Indexer for vector elements.
                /// </summary>
                property double default[int] {
                    double get(int i) {
                        if (i < 0 || i >= data->Length)
                            throw gcnew IndexOutOfRangeException("Vector index out of range.");
                        return data[i];
                    }
                    void set(int i, double value) {
                        if (i < 0 || i >= data->Length)
                            throw gcnew IndexOutOfRangeException("Vector index out of range.");
                        data[i] = value;
                    }
                    }

                    /// <summary>
                    /// Initializes a zero vector of given dimension.
                    /// </summary>
                    Vector(int dim) {
                    if (dim <= 0)
                        throw gcnew ArgumentException("Dimension must be positive.");
                    data = gcnew array<double>(dim);
                }

                /// <summary>
                /// Initializes vector from a managed array.
                /// </summary>
                Vector(array<double>^ values) {
                    if (values == nullptr || values->Length == 0)
                        throw gcnew ArgumentException("Values array cannot be null or empty.");
                    data = gcnew array<double>(values->Length);
                    Array::Copy(values, data, values->Length);
                }

                /// <summary>
                /// Creates a deep copy of the vector.
                /// </summary>
                Vector^ Clone() {
                    return gcnew Vector(data);
                }

                // ==============================
                // Basic operations
                // ==============================

                /// <summary>
                /// Adds two vectors of the same dimension.
                /// </summary>
                static Vector^ Add(Vector^ a, Vector^ b) {
                    if (a->Dimension != b->Dimension)
                        throw gcnew Agrine::Math::Core::DimensionMismatchException("Vector dimensions must match for addition.");

                    array<double>^ result = gcnew array<double>(a->Dimension);
                    for (int i = 0; i < a->Dimension; i++)
                        result[i] = a[i] + b[i];
                    return gcnew Vector(result);
                }

                /// <summary>
                /// Subtracts two vectors of the same dimension.
                /// </summary>
                static Vector^ Subtract(Vector^ a, Vector^ b) {
                    if (a->Dimension != b->Dimension)
                        throw gcnew Agrine::Math::Core::DimensionMismatchException("Vector dimensions must match for subtraction.");

                    array<double>^ result = gcnew array<double>(a->Dimension);
                    for (int i = 0; i < a->Dimension; i++)
                        result[i] = a[i] - b[i];
                    return gcnew Vector(result);
                }

                /// <summary>
                /// Scalar multiplication.
                /// </summary>
                static Vector^ Multiply(Vector^ v, double scalar) {
                    array<double>^ result = gcnew array<double>(v->Dimension);
                    for (int i = 0; i < v->Dimension; i++)
                        result[i] = v[i] * scalar;
                    return gcnew Vector(result);
                }

                /// <summary>
                /// Dot product of two vectors.
                /// </summary>
                static double Dot(Vector^ a, Vector^ b) {
                    if (a->Dimension != b->Dimension)
                        throw gcnew Agrine::Math::Core::DimensionMismatchException("Vector dimensions must match for dot product.");

                    double sum = 0.0;
                    for (int i = 0; i < a->Dimension; i++)
                        sum += a[i] * b[i];
                    return sum;
                }

                // ==============================
                // Advanced operations
                // ==============================

                /// <summary>
                /// Returns the Euclidean norm (magnitude) of the vector.
                /// </summary>
                double Norm() {
                    double sumSq = 0.0;
                    for (int i = 0; i < Dimension; i++)
                        sumSq += data[i] * data[i];
                    return System::Math::Sqrt(sumSq);
                }

                /// <summary>
                /// Normalizes the vector to unit length.
                /// </summary>
                Vector^ Normalize() {
                    double n = Norm();
                    if (Core::Utils::IsZero(n))
                        throw gcnew DivideByZeroException("Cannot normalize a zero vector.");
                    return Multiply(this, 1.0 / n);
                }

                /// <summary>
                /// Computes the angle (in radians) between two vectors.
                /// </summary>
                static double Angle(Vector^ a, Vector^ b) {
                    double dot = Dot(a, b);
                    double denom = a->Norm() * b->Norm();
                    if (Core::Utils::IsZero(denom))
                        throw gcnew DivideByZeroException("Cannot compute angle with zero vector.");
                    double cosTheta = dot / denom;

                    // Clamp to [-1, 1] to avoid NaN due to rounding
                    cosTheta = Core::Utils::Clamp(cosTheta, -1.0, 1.0);
                    return System::Math::Acos(cosTheta);
                }

                /// <summary>
                /// Checks if two vectors are approximately equal within tolerance.
                /// </summary>
                static bool AreEqual(Vector^ a, Vector^ b, double tol = Agrine::Math::Core::Constants::Tolerance) {
                    if (a->Dimension != b->Dimension) return false;
                    for (int i = 0; i < a->Dimension; i++) {
                        if (!Core::Utils::AreEqual(a[i], b[i], tol))
                            return false;
                    }
                    return true;
                }

                /// <summary>
                /// Returns string representation of vector.
                /// </summary>
                virtual String^ ToString() override {
                    return String::Format("Vector({0})", String::Join(", ", data));
                }
            };
        }
    }
}

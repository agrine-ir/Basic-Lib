#include "Vector.h"

using namespace System;
using namespace Agrine::Math::Core;

namespace Agrine {
    namespace Math {
        namespace Algebra {

            // Constructor: initialize vector with given dimension
            Vector::Vector(int dimension)
            {
                if (dimension <= 0)
                    throw gcnew ArgumentException("Dimension must be positive.");

                this->dimension = dimension;
                data = gcnew array<double>(dimension);

                for (int i = 0; i < dimension; i++)
                    data[i] = 0.0;
            }

            // Constructor: initialize vector with given array
            Vector::Vector(array<double>^ values)
            {
                if (values == nullptr || values->Length == 0)
                    throw gcnew ArgumentException("Values array cannot be null or empty.");

                dimension = values->Length;
                data = gcnew array<double>(dimension);
                Array::Copy(values, data, dimension);
            }

            // Get dimension of vector
            int Vector::Dimension::get()
            {
                return dimension;
            }

            // Indexer
            double Vector::default::get(int index)
            {
                if (index < 0 || index >= dimension)
                    throw gcnew IndexOutOfRangeException("Index out of range.");
                return data[index];
            }

            void Vector::default::set(int index, double value)
            {
                if (index < 0 || index >= dimension)
                    throw gcnew IndexOutOfRangeException("Index out of range.");
                data[index] = value;
            }

            // Vector addition
            Vector^ Vector::Add(Vector^ other)
            {
                if (this->dimension != other->dimension)
                    throw gcnew DimensionMismatchException("Vector dimensions must match for addition.");

                array<double>^ result = gcnew array<double>(dimension);
                for (int i = 0; i < dimension; i++)
                    result[i] = this->data[i] + other[i];

                return gcnew Vector(result);
            }

            // Vector subtraction
            Vector^ Vector::Subtract(Vector^ other)
            {
                if (this->dimension != other->dimension)
                    throw gcnew DimensionMismatchException("Vector dimensions must match for subtraction.");

                array<double>^ result = gcnew array<double>(dimension);
                for (int i = 0; i < dimension; i++)
                    result[i] = this->data[i] - other[i];

                return gcnew Vector(result);
            }

            // Dot product
            double Vector::Dot(Vector^ other)
            {
                if (this->dimension != other->dimension)
                    throw gcnew DimensionMismatchException("Vector dimensions must match for dot product.");

                double sum = 0.0;
                for (int i = 0; i < dimension; i++)
                    sum += this->data[i] * other[i];

                return sum;
            }

            // Norm (magnitude)
            double Vector::Norm()
            {
                double sum = 0.0;
                for (int i = 0; i < dimension; i++)
                    sum += data[i] * data[i];

                return Math::Sqrt(sum);
            }

            // Equality check with tolerance
            bool Vector::AreEqual(Vector^ a, Vector^ b, double tol)
            {
                if (a->Dimension != b->Dimension)
                    return false;

                for (int i = 0; i < a->Dimension; i++)
                {
                    if (Math::Abs(a[i] - b[i]) > tol)
                        return false;
                }

                return true;
            }

            // ToString override for nice output
            String^ Vector::ToString()
            {
                array<String^>^ elems = gcnew array<String^>(dimension);
                for (int i = 0; i < dimension; i++)
                    elems[i] = data[i].ToString("G5");

                return "[" + String::Join(", ", elems) + "]";
            }

        } // namespace Algebra
    } // namespace Math
} // namespace Agrine

#include "Vector.h"

using namespace System;
using namespace Agrine::Math::Core;

namespace Agrine {
    namespace Math {
        namespace Algebra {

            // ===== Constructors =====
            Vector::Vector(int dim)
            {
                if (dim <= 0)
                    throw gcnew ArgumentException("Dimension must be positive.");

                this->dimension = dim;
                this->data = gcnew array<double>(dim);
            }

            Vector::Vector(array<double>^ values)
            {
                if (values == nullptr || values->Length == 0)
                    throw gcnew ArgumentException("Values array cannot be null or empty.");

                this->dimension = values->Length;
                this->data = gcnew array<double>(dimension);
                Array::Copy(values, this->data, dimension);
            }

            // ===== Properties =====
            int Vector::Dimension::get()
            {
                return this->dimension;
            }

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

            // ===== Utility =====
            Vector^ Vector::Clone()
            {
                return gcnew Vector(this->data);
            }

            // ===== Basic operations =====
            Vector^ Vector::Add(Vector^ other)
            {
                if (this->dimension != other->dimension)
                    throw gcnew DimensionMismatchException("Vector dimensions must match for addition.");

                array<double>^ result = gcnew array<double>(dimension);
                for (int i = 0; i < dimension; i++)
                    result[i] = this->data[i] + other[i];

                return gcnew Vector(result);
            }

            Vector^ Vector::Subtract(Vector^ other)
            {
                if (this->dimension != other->dimension)
                    throw gcnew DimensionMismatchException("Vector dimensions must match for subtraction.");

                array<double>^ result = gcnew array<double>(dimension);
                for (int i = 0; i < dimension; i++)
                    result[i] = this->data[i] - other[i];

                return gcnew Vector(result);
            }

            Vector^ Vector::Multiply(double scalar)
            {
                array<double>^ result = gcnew array<double>(dimension);
                for (int i = 0; i < dimension; i++)
                    result[i] = this->data[i] * scalar;

                return gcnew Vector(result);
            }

            double Vector::Dot(Vector^ other)
            {
                if (this->dimension != other->dimension)
                    throw gcnew DimensionMismatchException("Vector dimensions must match for dot product.");

                double sum = 0.0;
                for (int i = 0; i < dimension; i++)
                    sum += this->data[i] * other[i];

                return sum;
            }

            // ===== Advanced operations =====
            double Vector::Norm()
            {
                double sum = 0.0;
                for (int i = 0; i < dimension; i++)
                    sum += data[i] * data[i];

                return System::Math::Sqrt(sum);
            }

            Vector^ Vector::Normalize()
            {
                double n = Norm();
                if (Utils::IsZero(n))
                    throw gcnew DivideByZeroException("Cannot normalize a zero vector.");

                return Multiply(1.0 / n);
            }

            double Vector::Angle(Vector^ a, Vector^ b)
            {
                double dot = a->Dot(b);
                double denom = a->Norm() * b->Norm();
                if (Utils::IsZero(denom))
                    throw gcnew DivideByZeroException("Cannot compute angle with zero vector.");

                double cosTheta = dot / denom;
                cosTheta = Utils::Clamp(cosTheta, -1.0, 1.0); // safe clamp
                return System::Math::Acos(cosTheta);
            }

            bool Vector::AreEqual(Vector^ a, Vector^ b, double tol)
            {
                if (a->Dimension != b->Dimension)
                    return false;

                for (int i = 0; i < a->Dimension; i++)
                {
                    if (!Utils::AreEqual(a[i], b[i], tol))
                        return false;
                }
                return true;
            }

            // ===== ToString =====
            String^ Vector::ToString()
            {
                array<String^>^ elems = gcnew array<String^>(dimension);
                for (int i = 0; i < dimension; i++)
                    elems[i] = data[i].ToString("G5");

                return String::Format("Vector({0})", String::Join(", ", elems));
            }

        } // namespace Algebra
    } // namespace Math
} // namespace Agrine

#pragma once

// Vector.h - Algebra module
// Part of Agrine.Math (C++/CLI, .NET Framework 4.6.2)

#include "Types.h"
#include "Constants.h"
#include "Utils.h"
#include "Exceptions.h"

using namespace System;
using namespace Agrine::Math::Core;

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
                int dimension;
                array<double>^ data;

            public:
                // ===== Properties =====
                property int Dimension{
                    int get();
                }

                    property double default[int] {
                    double get(int index);
                    void set(int index, double value);
                    }

                    // ===== Constructors =====
                Vector(int dim);
                Vector(array<double>^ values);

                // ===== Utility =====
                Vector^ Clone();

                // ===== Basic operations =====
                Vector^ Add(Vector^ other);
                Vector^ Subtract(Vector^ other);
                Vector^ Multiply(double scalar);
                double Dot(Vector^ other);

                // ===== Advanced operations =====
                double Norm();
                Vector^ Normalize();
                static double Angle(Vector^ a, Vector^ b);
                static bool AreEqual(Vector^ a, Vector^ b, double tol = Constants::Tolerance);

                // ===== ToString override =====
                virtual String^ ToString() override;
            };
        }
    }
}

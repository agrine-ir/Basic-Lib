#pragma once

// Core types and constants for Agrine.Math
// Project: Agrine.Math (C++/CLI, .NET Framework 4.6.2)

using namespace System;

namespace Agrine {
    namespace Math {
        namespace Core {

            /// <summary>
            /// Common mathematical constants.
            /// </summary>
            public ref class Constants abstract sealed
            {
            public:
                literal double Pi = 3.14159265358979323846;
                literal double E = 2.71828182845904523536;
                literal double GoldenR = 1.61803398874989484820; // φ
                literal double Sqrt2 = 1.41421356237309504880;
                literal double Epsilon = 1e-12; // Precision threshold
            };

            /// <summary>
            /// Basic 2D point structure.
            /// Useful for geometry and vector operations.
            /// </summary>
            public value struct Point2D
            {
            public:
                double X;
                double Y;

                Point2D(double x, double y)
                    : X(x), Y(y) {
                }

                virtual String^ ToString() override {
                    return String::Format("({0}, {1})", X, Y);
                }
            };

            /// <summary>
            /// Basic 3D point structure.
            /// </summary>
            public value struct Point3D
            {
            public:
                double X;
                double Y;
                double Z;

                Point3D(double x, double y, double z)
                    : X(x), Y(y), Z(z) {
                }

                virtual String^ ToString() override {
                    return String::Format("({0}, {1}, {2})", X, Y, Z);
                }
            };

            /// <summary>
            /// A lightweight Complex number struct.
            /// Avoids dependency on System::Numerics::Complex.
            /// </summary>
            public value struct Complex
            {
            public:
                double Real;
                double Imag;

                Complex(double r, double i)
                    : Real(r), Imag(i) {
                }

                // Magnitude (absolute value)
                double Magnitude() {
                    return System::Math::Sqrt(Real * Real + Imag * Imag);
                }

                // Conjugate of the complex number
                Complex Conjugate() {
                    return Complex(Real, -Imag);
                }

                virtual String^ ToString() override {
                    return String::Format("{0} + {1}i", Real, Imag);
                }
            };
        }
    }
}

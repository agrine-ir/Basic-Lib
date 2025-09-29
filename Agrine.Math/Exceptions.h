#pragma once

// Custom math exceptions for Agrine.Math
// Project: Agrine.Math (C++/CLI, .NET Framework 4.6.2)

using namespace System;

namespace Agrine {
    namespace Math {
        namespace Core {

            /// <summary>
            /// Base class for all custom math exceptions in Agrine.Math.
            /// </summary>
            [Serializable]
            public ref class MathException : public Exception
            {
            public:
                MathException() : Exception("A math error occurred.") {}
                MathException(String^ message) : Exception(message) {}
                MathException(String^ message, Exception^ inner) : Exception(message, inner) {}
            };

            /// <summary>
            /// Exception thrown when dividing by zero.
            /// </summary>
            [Serializable]
                public ref class DivisionByZeroException : public MathException
            {
            public:
                DivisionByZeroException()
                    : MathException("Division by zero is not allowed.") {
                }
                DivisionByZeroException(String^ message)
                    : MathException(message) {
                }
            };

            /// <summary>
            /// Exception thrown when a matrix or vector dimension is invalid.
            /// </summary>
            [Serializable]
                public ref class InvalidDimensionException : public MathException
            {
            public:
                InvalidDimensionException()
                    : MathException("Invalid dimensions for this operation.") {
                }
                InvalidDimensionException(String^ message)
                    : MathException(message) {
                }
            };

            /// <summary>
            /// Exception thrown when a numerical method fails to converge.
            /// </summary>
            [Serializable]
                public ref class ConvergenceException : public MathException
            {
            public:
                ConvergenceException()
                    : MathException("The numerical method failed to converge.") {
                }
                ConvergenceException(String^ message)
                    : MathException(message) {
                }
            };

            /// <summary>
            /// Exception thrown when a calculation overflows or exceeds limits.
            /// </summary>
            [Serializable]
                public ref class OverflowException : public MathException
            {
            public:
                OverflowException()
                    : MathException("Mathematical overflow occurred.") {
                }
                OverflowException(String^ message)
                    : MathException(message) {
                }
            };
        }
    }
}

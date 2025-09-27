#pragma once

// ============================================================
// MathCore.h - Core Header for Agrine.Math Library
// ============================================================
// This header defines the base macros, common includes,
// and gathers all the public module headers (Algebra, Geometry, ...).
// It ensures that functions can be exported/imported correctly
// when building or consuming the DLL in C++ or C# (.NET).
// ============================================================

//
// DLL Export/Import Macro
// ------------------------
// AGRINE_API ensures functions are correctly exported when building the DLL
// and imported when using the DLL from external projects (e.g., C#, C++).
//
// - When building the DLL, the symbol AGRINE_MATH_EXPORTS must be defined
//   (Visual Studio sets this automatically in project properties).
// - When consuming the DLL, this symbol should NOT be defined.
//
#ifdef AGRINE_MATH_EXPORTS
#define AGRINE_API extern "C" __declspec(dllexport)
#else
#define AGRINE_API extern "C" __declspec(dllimport)
#endif

//
// Common includes
// These headers are available to all modules in the library.
//
#include <cstdint>   // For fixed-width integer types
#include <cstddef>   // For size_t
#include <cmath>     // For mathematical functions

//
// Module includes
// Each of these headers declares functions for a specific area of math.
//
#include "Constants.h"
#include "Algebra.h"
#include "Geometry.h"
#include "Linalg.h"       // or LinearAlgebra.h if you keep the old name
#include "Statistics.h"
#include "Complex.h"
#include "Utils.h"

//
// Common types and error codes
// These are shared across all modules.
//
enum class MathError
{
    None = 0,          // No error
    DivisionByZero,    // Division by zero occurred
    Overflow,          // Arithmetic overflow
    InvalidArgument    // Invalid function argument
};

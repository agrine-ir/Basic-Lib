#pragma once

// ============================================================
// Geometry.h - Geometry Module for Agrine.Math Library
// ============================================================
//
// This header defines advanced geometry functions, including
// distance calculations, area computations, vector operations,
// coordinate transformations, and more.
//
// Uses mathematical constants from Agrine.Math.Constants.
// ============================================================

#include "MathCore.h"
#include "Constants.h"
#include <cmath>
#include <array>

namespace Agrine
{
    namespace Math
    {
        namespace Geometry
        {
            // ----------------------------
            // Basic Geometry Functions
            // ----------------------------

            /// <summary>
            /// Calculate the Euclidean distance between two points in 2D space.
            /// </summary>
            AGRINE_API double Distance2D(double x1, double y1, double x2, double y2);

            /// <summary>
            /// Calculate the Euclidean distance between two points in 3D space.
            /// </summary>
            AGRINE_API double Distance3D(double x1, double y1, double z1,
                double x2, double y2, double z2);

            /// <summary>
            /// Calculate the area of a triangle given its three vertices in 2D.
            /// </summary>
            AGRINE_API double TriangleArea(double x1, double y1,
                double x2, double y2,
                double x3, double y3);

            // ----------------------------
            // Advanced Geometry Functions
            // ----------------------------

            /// <summary>
            /// Rotate a 2D point (x, y) around the origin by a given angle in degrees.
            /// </summary>
            AGRINE_API std::array<double, 2> RotatePoint2D(double x, double y, double angleDegrees);

            /// <summary>
            /// Convert Cartesian coordinates (x, y, z) to spherical coordinates (r, theta, phi).
            /// Theta: inclination, Phi: azimuth.
            /// </summary>
            AGRINE_API std::array<double, 3> CartesianToSpherical(double x, double y, double z);

            /// <summary>
            /// Convert spherical coordinates (r, theta, phi) to Cartesian coordinates (x, y, z).
            /// </summary>
            AGRINE_API std::array<double, 3> SphericalToCartesian(double r, double theta, double phi);

            /// <summary>
            /// Calculate the angle between two vectors in degrees.
            /// </summary>
            AGRINE_API double AngleBetweenVectors(const std::array<double, 3>& v1,
                const std::array<double, 3>& v2);
        }
    }
}

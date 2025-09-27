#pragma once
#include "MathCore.h"
#include "Constants.h"
#include <array>

namespace Agrine
{
    namespace Math
    {
        namespace Geometry
        {
            // ================================
            // Geometry - Advanced 2D & 3D Operations
            // ================================

            // ----------------------------
            // Distance Calculations
            // ----------------------------
            AGRINE_API double Distance2D(double x1, double y1, double x2, double y2);
            AGRINE_API double Distance3D(double x1, double y1, double z1,
                double x2, double y2, double z2);

            // ----------------------------
            // Triangle Calculations
            // ----------------------------
            AGRINE_API double TriangleArea(double a, double b, double c); // Using Heron's formula

            // ----------------------------
            // Volume and Surface Area Calculations
            // ----------------------------
            AGRINE_API double CubeVolume(double side);
            AGRINE_API double CubeSurfaceArea(double side);

            AGRINE_API double SphereVolume(double radius);
            AGRINE_API double SphereSurfaceArea(double radius);

            AGRINE_API double CylinderVolume(double radius, double height);
            AGRINE_API double CylinderSurfaceArea(double radius, double height);

            AGRINE_API double ConeVolume(double radius, double height);
            AGRINE_API double ConeSurfaceArea(double radius, double height);

            AGRINE_API double RectangularPrismVolume(double length, double width, double height);
            AGRINE_API double RectangularPrismSurfaceArea(double length, double width, double height);

            // ----------------------------
            // Vector and 3D Operations
            // ----------------------------
            AGRINE_API std::array<double, 3> CrossProduct(const std::array<double, 3>& v1,
                const std::array<double, 3>& v2);

            AGRINE_API std::array<double, 3> NormalizeVector(const std::array<double, 3>& v);

            AGRINE_API double DotProduct(const std::array<double, 3>& v1,
                const std::array<double, 3>& v2);

            AGRINE_API std::array<double, 3> ProjectVector(const std::array<double, 3>& v,
                const std::array<double, 3>& onto);

            // ----------------------------
            // 3D Rotations
            // ----------------------------
            AGRINE_API std::array<double, 3> RotateAroundX(const std::array<double, 3>& point, double angleDegrees);
            AGRINE_API std::array<double, 3> RotateAroundY(const std::array<double, 3>& point, double angleDegrees);
            AGRINE_API std::array<double, 3> RotateAroundZ(const std::array<double, 3>& point, double angleDegrees);

            // ----------------------------
            // Plane Operations
            // ----------------------------
            AGRINE_API double DistancePointToPlane(const std::array<double, 3>& point,
                const std::array<double, 3>& planePoint,
                const std::array<double, 3>& planeNormal);

        } // namespace Geometry
    }
}

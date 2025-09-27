#include "Geometry.h"
#include <cmath>

using namespace Agrine::Math;
using namespace Agrine::Math::Constants;

// ----------------------------
// Distance Calculations
// ----------------------------
double Geometry::Distance2D(double x1, double y1, double x2, double y2)
{
    return std::sqrt(std::pow(x2 - x1, 2) + std::pow(y2 - y1, 2));
}

double Geometry::Distance3D(double x1, double y1, double z1,
    double x2, double y2, double z2)
{
    return std::sqrt(std::pow(x2 - x1, 2) +
        std::pow(y2 - y1, 2) +
        std::pow(z2 - z1, 2));
}

// ----------------------------
// Triangle Calculations
// ----------------------------
double Geometry::TriangleArea(double a, double b, double c)
{
    double s = (a + b + c) / 2.0;
    double areaSquared = s * (s - a) * (s - b) * (s - c);
    if (areaSquared < 0) return 0.0; // Handle invalid triangles
    return std::sqrt(areaSquared);
}

// ----------------------------
// Volume and Surface Area Calculations
// ----------------------------
double Geometry::CubeVolume(double side) { return side * side * side; }
double Geometry::CubeSurfaceArea(double side) { return 6 * side * side; }

double Geometry::SphereVolume(double radius) { return 4.0 / 3.0 * Constants::PI * std::pow(radius, 3); }
double Geometry::SphereSurfaceArea(double radius) { return 4 * Constants::PI * radius * radius; }

double Geometry::CylinderVolume(double radius, double height) { return Constants::PI * radius * radius * height; }
double Geometry::CylinderSurfaceArea(double radius, double height) { return 2 * Constants::PI * radius * (radius + height); }

double Geometry::ConeVolume(double radius, double height) { return Constants::PI * radius * radius * height / 3.0; }
double Geometry::ConeSurfaceArea(double radius, double height) { return Constants::PI * radius * (radius + std::sqrt(height * height + radius * radius)); }

double Geometry::RectangularPrismVolume(double length, double width, double height) { return length * width * height; }
double Geometry::RectangularPrismSurfaceArea(double length, double width, double height) { return 2 * (length * width + width * height + length * height); }

// ----------------------------
// Vector and 3D Operations
// ----------------------------
std::array<double, 3> Geometry::CrossProduct(const std::array<double, 3>& v1,
    const std::array<double, 3>& v2)
{
    return {
        v1[1] * v2[2] - v1[2] * v2[1],
        v1[2] * v2[0] - v1[0] * v2[2],
        v1[0] * v2[1] - v1[1] * v2[0]
    };
}

std::array<double, 3> Geometry::NormalizeVector(const std::array<double, 3>& v)
{
    double mag = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (mag == 0) return { 0.0, 0.0, 0.0 };
    return { v[0] / mag, v[1] / mag, v[2] / mag };
}

double Geometry::DotProduct(const std::array<double, 3>& v1,
    const std::array<double, 3>& v2)
{
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

std::array<double, 3> Geometry::ProjectVector(const std::array<double, 3>& v,
    const std::array<double, 3>& onto)
{
    double dot = DotProduct(v, onto);
    double magSq = DotProduct(onto, onto);
    if (magSq == 0) return { 0.0, 0.0, 0.0 };
    return { (dot / magSq) * onto[0], (dot / magSq) * onto[1], (dot / magSq) * onto[2] };
}

// ----------------------------
// 3D Rotations
// ----------------------------
std::array<double, 3> Geometry::RotateAroundX(const std::array<double, 3>& point, double angleDegrees)
{
    double rad = angleDegrees * Constants::DEG_TO_RAD;
    double cosA = std::cos(rad), sinA = std::sin(rad);
    return { point[0], point[1] * cosA - point[2] * sinA, point[1] * sinA + point[2] * cosA };
}

std::array<double, 3> Geometry::RotateAroundY(const std::array<double, 3>& point, double angleDegrees)
{
    double rad = angleDegrees * Constants::DEG_TO_RAD;
    double cosA = std::cos(rad), sinA = std::sin(rad);
    return { point[0] * cosA + point[2] * sinA, point[1], -point[0] * sinA + point[2] * cosA };
}

std::array<double, 3> Geometry::RotateAroundZ(const std::array<double, 3>& point, double angleDegrees)
{
    double rad = angleDegrees * Constants::DEG_TO_RAD;
    double cosA = std::cos(rad), sinA = std::sin(rad);
    return { point[0] * cosA - point[1] * sinA, point[0] * sinA + point[1] * cosA, point[2] };
}

// ----------------------------
// Plane Operations
// ----------------------------
double Geometry::DistancePointToPlane(const std::array<double, 3>& point,
    const std::array<double, 3>& planePoint,
    const std::array<double, 3>& planeNormal)
{
    auto normalizedNormal = NormalizeVector(planeNormal);
    return std::abs(DotProduct(normalizedNormal, std::array<double, 3>{
        point[0] - planePoint[0],
            point[1] - planePoint[1],
            point[2] - planePoint[2]
    }));
}

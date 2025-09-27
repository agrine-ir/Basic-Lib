#include "Geometry.h"

using namespace Agrine::Math;
using namespace Agrine::Math::Constants;

// ----------------------------
// Basic Geometry Implementations
// ----------------------------

double Geometry::Distance2D(double x1, double y1, double x2, double y2)
{
    return std::sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

double Geometry::Distance3D(double x1, double y1, double z1,
    double x2, double y2, double z2)
{
    return std::sqrt((x2 - x1) * (x2 - x1) +
        (y2 - y1) * (y2 - y1) +
        (z2 - z1) * (z2 - z1));
}

double Geometry::TriangleArea(double x1, double y1,
    double x2, double y2,
    double x3, double y3)
{
    return std::abs((x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2.0);
}

// ----------------------------
// Advanced Geometry Implementations
// ----------------------------

std::array<double, 2> Geometry::RotatePoint2D(double x, double y, double angleDegrees)
{
    double radians = angleDegrees * Constants::DEG_TO_RAD;
    double cosA = std::cos(radians);
    double sinA = std::sin(radians);
    return { x * cosA - y * sinA, x * sinA + y * cosA };
}

std::array<double, 3> Geometry::CartesianToSpherical(double x, double y, double z)
{
    double r = std::sqrt(x * x + y * y + z * z);
    double theta = std::acos(z / r); // inclination
    double phi = std::atan2(y, x);   // azimuth
    return { r, theta, phi };
}

std::array<double, 3> Geometry::SphericalToCartesian(double r, double theta, double phi)
{
    double sinTheta = std::sin(theta);
    return {
        r * sinTheta * std::cos(phi),
        r * sinTheta * std::sin(phi),
        r * std::cos(theta)
    };
}

double Geometry::AngleBetweenVectors(const std::array<double, 3>& v1,
    const std::array<double, 3>& v2)
{
    double dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    double mag1 = std::sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
    double mag2 = std::sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
    return std::acos(dot / (mag1 * mag2)) * Constants::RAD_TO_DEG;
}

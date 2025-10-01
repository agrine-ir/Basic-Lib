#include "Euclidean.h"

using namespace System;
using namespace Agrine::Math::Core;

namespace Agrine {
    namespace Math {
        namespace Geometry {

            // -----------------------------
            // Basic distances
            // -----------------------------
            double Distance2D(Point2D a, Point2D b)
            {
                double dx = a.X - b.X;
                double dy = a.Y - b.Y;
                return System::Math::Sqrt(dx * dx + dy * dy);
            }

            double Distance3D(Point3D a, Point3D b)
            {
                double dx = a.X - b.X;
                double dy = a.Y - b.Y;
                double dz = a.Z - b.Z;
                return System::Math::Sqrt(dx * dx + dy * dy + dz * dz);
            }

            // -----------------------------
            // Closest point on segment (2D) and distance
            // -----------------------------
            Point2D ClosestPointOnSegment2D(Point2D p, Point2D a, Point2D b)
            {
                // If a == b return a
                if (Utils::AreEqual(a.X, b.X) && Utils::AreEqual(a.Y, b.Y))
                    return a;

                double vx = b.X - a.X;
                double vy = b.Y - a.Y;
                double wx = p.X - a.X;
                double wy = p.Y - a.Y;

                double c1 = vx * wx + vy * wy;
                double c2 = vx * vx + vy * vy;

                double t = c1 / c2;
                if (t <= 0.0) return a;
                if (t >= 1.0) return b;

                return Point2D(a.X + t * vx, a.Y + t * vy);
            }

            double DistancePointToSegment2D(Point2D p, Point2D a, Point2D b)
            {
                Point2D cp = ClosestPointOnSegment2D(p, a, b);
                return Distance2D(p, cp);
            }

            // -----------------------------
            // 2D triangle & polygon area (shoelace)
            // -----------------------------
            double TriangleArea2D(Point2D a, Point2D b, Point2D c)
            {
                // signed area = 0.5 * ((x1(y2 - y3) + x2(y3 - y1) + x3(y1 - y2)))
                double val = a.X * (b.Y - c.Y) + b.X * (c.Y - a.Y) + c.X * (a.Y - b.Y);
                return System::Math::Abs(0.5 * val);
            }

            double PolygonArea2D(array<Point2D>^ polygon)
            {
                if (polygon == nullptr) throw gcnew ArgumentNullException("polygon");
                int n = polygon->Length;
                if (n < 3) return 0.0;

                double sum = 0.0;
                for (int i = 0; i < n; i++)
                {
                    Point2D p1 = polygon[i];
                    Point2D p2 = polygon[(i + 1) % n];
                    sum += p1.X * p2.Y - p2.X * p1.Y;
                }
                return System::Math::Abs(0.5 * sum);
            }

            Point2D PolygonCentroid2D(array<Point2D>^ polygon)
            {
                if (polygon == nullptr) throw gcnew ArgumentNullException("polygon");
                int n = polygon->Length;
                if (n == 0) throw gcnew InvalidArgumentException("polygon must have at least one vertex.");
                if (n == 1) return polygon[0];
                if (n == 2) return Point2D((polygon[0].X + polygon[1].X) / 2.0, (polygon[0].Y + polygon[1].Y) / 2.0);

                double A = 0.0; // signed area * 2
                double Cx = 0.0, Cy = 0.0;

                for (int i = 0; i < n; i++)
                {
                    Point2D p0 = polygon[i];
                    Point2D p1 = polygon[(i + 1) % n];
                    double cross = p0.X * p1.Y - p1.X * p0.Y;
                    A += cross;
                    Cx += (p0.X + p1.X) * cross;
                    Cy += (p0.Y + p1.Y) * cross;
                }

                double area = 0.5 * A;
                if (Utils::IsZero(area))
                {
                    // degenerate polygon: fallback to average of points
                    double sx = 0.0, sy = 0.0;
                    for each (Point2D p in polygon) { sx += p.X; sy += p.Y; }
                    return Point2D(sx / n, sy / n);
                }

                Cx /= (6.0 * area);
                Cy /= (6.0 * area);
                return Point2D(Cx, Cy);
            }

            // -----------------------------
            // Point in polygon: ray casting
            // returns true if point inside or on edge
            // -----------------------------
            bool PointInPolygon(array<Point2D>^ polygon, Point2D p)
            {
                if (polygon == nullptr) throw gcnew ArgumentNullException("polygon");
                int n = polygon->Length;
                if (n < 3) return false;

                bool inside = false;
                for (int i = 0, j = n - 1; i < n; j = i++)
                {
                    Point2D pi = polygon[i];
                    Point2D pj = polygon[j];

                    // check if point is on the edge (between pi and pj)
                    double minx = System::Math::Min(pi.X, pj.X);
                    double maxx = System::Math::Max(pi.X, pj.X);
                    double miny = System::Math::Min(pi.Y, pj.Y);
                    double maxy = System::Math::Max(pi.Y, pj.Y);
                    // collinearity test
                    double cross = (pj.X - pi.X) * (p.Y - pi.Y) - (pj.Y - pi.Y) * (p.X - pi.X);
                    if (Utils::IsZero(cross) && p.X >= minx - Constants::Epsilon && p.X <= maxx + Constants::Epsilon
                        && p.Y >= miny - Constants::Epsilon && p.Y <= maxy + Constants::Epsilon)
                        return true; // on edge

                    bool intersect = ((pi.Y > p.Y) != (pj.Y > p.Y)) &&
                        (p.X < (pj.X - pi.X) * (p.Y - pi.Y) / (pj.Y - pi.Y + 1e-20) + pi.X);
                    if (intersect) inside = !inside;
                }
                return inside;
            }

            // -----------------------------
            // 3D triangle area & tetrahedron volume
            // -----------------------------
            // cross product (a x b)
            static void Cross3D(double ax, double ay, double az, double bx, double by, double bz, double& rx, double& ry, double& rz)
            {
                rx = ay * bz - az * by;
                ry = az * bx - ax * bz;
                rz = ax * by - ay * bx;
            }

            double TriangleArea3D(Point3D a, Point3D b, Point3D c)
            {
                double ux = b.X - a.X, uy = b.Y - a.Y, uz = b.Z - a.Z;
                double vx = c.X - a.X, vy = c.Y - a.Y, vz = c.Z - a.Z;
                double rx, ry, rz;
                Cross3D(ux, uy, uz, vx, vy, vz, rx, ry, rz);
                double norm = System::Math::Sqrt(rx * rx + ry * ry + rz * rz);
                return 0.5 * norm;
            }

            double TetrahedronVolume(Point3D a, Point3D b, Point3D c, Point3D d)
            {
                // volume = | dot( (b-a) x (c-a), (d-a) ) | / 6
                double ux = b.X - a.X, uy = b.Y - a.Y, uz = b.Z - a.Z;
                double vx = c.X - a.X, vy = c.Y - a.Y, vz = c.Z - a.Z;
                double wx = d.X - a.X, wy = d.Y - a.Y, wz = d.Z - a.Z;
                double rx, ry, rz;
                Cross3D(ux, uy, uz, vx, vy, vz, rx, ry, rz);
                double dot = rx * wx + ry * wy + rz * wz;
                return System::Math::Abs(dot) / 6.0;
            }

            // -----------------------------
            // Plane implementation
            // -----------------------------
            void Plane::Normalize()
            {
                double norm = System::Math::Sqrt(A * A + B * B + C * C);
                if (!Utils::IsZero(norm))
                {
                    A /= norm; B /= norm; C /= norm; D /= norm;
                }
            }

            Plane::Plane(double a, double b, double c, double d)
            {
                if (Utils::IsZero(a) && Utils::IsZero(b) && Utils::IsZero(c))
                    throw gcnew InvalidArgumentException("Plane normal cannot be zero vector.");

                this->A = a; this->B = b; this->C = c; this->D = d;
                Normalize();
            }

            Plane::Plane(Point3D p1, Point3D p2, Point3D p3)
            {
                // normal = (p2 - p1) x (p3 - p1)
                double ux = p2.X - p1.X, uy = p2.Y - p1.Y, uz = p2.Z - p1.Z;
                double vx = p3.X - p1.X, vy = p3.Y - p1.Y, vz = p3.Z - p1.Z;
                double nx, ny, nz;
                Cross3D(ux, uy, uz, vx, vy, vz, nx, ny, nz);

                if (Utils::IsZero(nx) && Utils::IsZero(ny) && Utils::IsZero(nz))
                    throw gcnew InvalidArgumentException("Three points are collinear; cannot define a plane.");

                this->A = nx; this->B = ny; this->C = nz;
                // plane: A(x - x0) + B(y - y0) + C(z - z0) = 0 => A x + B y + C z + D = 0
                this->D = -(A * p1.X + B * p1.Y + C * p1.Z);
                Normalize();
            }

            double Plane::DistanceToPoint(Point3D p)
            {
                double num = System::Math::Abs(A * p.X + B * p.Y + C * p.Z + D);
                // since normalized denom is 1
                double denom = System::Math::Sqrt(A * A + B * B + C * C);
                if (Utils::IsZero(denom)) throw gcnew InvalidArgumentException("Invalid plane (zero normal).");
                return num / denom;
            }

            Point3D Plane::ProjectPointToPlane(Point3D p)
            {
                double denom = System::Math::Sqrt(A * A + B * B + C * C);
                if (Utils::IsZero(denom)) throw gcnew InvalidArgumentException("Invalid plane (zero normal).");

                double signedDist = (A * p.X + B * p.Y + C * p.Z + D) / denom;
                double nx = A / denom, ny = B / denom, nz = C / denom;
                return Point3D(p.X - signedDist * nx, p.Y - signedDist * ny, p.Z - signedDist * nz);
            }

            bool Plane::Contains(Point3D p)
            {
                double val = A * p.X + B * p.Y + C * p.Z + D;
                return Utils::IsZero(val);
            }

            double Plane::CoefA::get() { return A; }
            double Plane::CoefB::get() { return B; }
            double Plane::CoefC::get() { return C; }
            double Plane::CoefD::get() { return D; }

            String^ Plane::ToString()
            {
                return String::Format("{0}x + {1}y + {2}z + {3} = 0", A.ToString("G6"), B.ToString("G6"), C.ToString("G6"), D.ToString("G6"));
            }

            double DistancePointToPlane(Point3D p, Plane^ plane)
            {
                if (plane == nullptr) throw gcnew ArgumentNullException("plane");
                return plane->DistanceToPoint(p);
            }

            Point3D ClosestPointOnPlane(Point3D p, Plane^ plane)
            {
                if (plane == nullptr) throw gcnew ArgumentNullException("plane");
                return plane->ProjectPointToPlane(p);
            }

        } // namespace Geometry
    } // namespace Math
} // namespace Agrine

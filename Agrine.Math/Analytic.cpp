#include "Analytic.h"

using namespace System;
using namespace Agrine::Math::Core;

namespace Agrine {
    namespace Math {
        namespace Geometry {

            // -------------------------
            // Free helpers
            // -------------------------
            double Distance(Point2D a, Point2D b)
            {
                double dx = a.X - b.X;
                double dy = a.Y - b.Y;
                return System::Math::Sqrt(dx * dx + dy * dy);
            }

            Point2D Midpoint(Point2D a, Point2D b)
            {
                return Point2D((a.X + b.X) / 2.0, (a.Y + b.Y) / 2.0);
            }

            // -------------------------
            // Line implementation
            // -------------------------
            void Line::Normalize()
            {
                double norm = System::Math::Sqrt(A * A + B * B);
                if (!Utils::IsZero(norm)) {
                    A /= norm;
                    B /= norm;
                    C /= norm;
                }
            }

            Line::Line(double a, double b, double c)
            {
                if (Utils::IsZero(a) && Utils::IsZero(b))
                    throw gcnew InvalidArgumentException("Invalid line coefficients: A and B cannot both be zero.");

                this->A = a; this->B = b; this->C = c;
                Normalize();
            }

            Line::Line(Point2D p1, Point2D p2)
            {
                if (Utils::AreEqual(p1.X, p2.X) && Utils::AreEqual(p1.Y, p2.Y))
                    throw gcnew InvalidArgumentException("Two distinct points required to define a line.");

                // line through p1 and p2: (y2 - y1)x - (x2 - x1)y + (x2*y1 - y2*x1) = 0
                double a = p2.Y - p1.Y;
                double b = -(p2.X - p1.X);
                double c = (p2.X * p1.Y - p2.Y * p1.X);
                this->A = a; this->B = b; this->C = c;
                Normalize();
            }

            Line^ Line::FromPointSlope(Point2D p, double slope)
            {
                // y - y0 = slope (x - x0) => slope*x - y + (y0 - slope*x0) = 0
                double a = slope;
                double b = -1.0;
                double c = p.Y - slope * p.X;
                return gcnew Line(a, b, c);
            }

            double Line::CoefA::get() { return A; }
            double Line::CoefB::get() { return B; }
            double Line::CoefC::get() { return C; }

            double Line::DistanceToPoint(Point2D p)
            {
                // distance = |A x0 + B y0 + C| / sqrt(A^2 + B^2)
                // since normalized, denom = 1
                double num = System::Math::Abs(A * p.X + B * p.Y + C);
                // if not normalized (should be normalized), divide by norm
                double denom = System::Math::Sqrt(A * A + B * B);
                if (Utils::IsZero(denom)) throw gcnew InvalidArgumentException("Invalid line (zero normal).");
                return num / denom;
            }

            Point2D Line::ProjectPoint(Point2D p)
            {
                // projection formula: point on line closest to p
                // If line normalized (A,B) unit: projection = p - (A*x0 + B*y0 + C) * (A,B)
                double denom = System::Math::Sqrt(A * A + B * B);
                if (Utils::IsZero(denom)) throw gcnew InvalidArgumentException("Invalid line (zero normal).");

                double distSigned = (A * p.X + B * p.Y + C) / denom; // signed distance
                // unit normal:
                double nx = A / denom;
                double ny = B / denom;
                // projected point = p - distSigned * normal
                return Point2D(p.X - distSigned * nx, p.Y - distSigned * ny);
            }

            bool Line::Contains(Point2D p)
            {
                double val = A * p.X + B * p.Y + C;
                return Utils::IsZero(val);
            }

            bool Line::IntersectionWith(Line^ other, Point2D% intersection)
            {
                if (other == nullptr) throw gcnew ArgumentNullException("other");

                // Solve:
                // A1 x + B1 y = -C1
                // A2 x + B2 y = -C2
                double a1 = this->A, b1 = this->B, c1 = this->C;
                double a2 = other->A, b2 = other->B, c2 = other->C;
                double det = a1 * b2 - a2 * b1;

                if (Utils::IsZero(det)) {
                    // parallel or coincident
                    return false;
                }

                double x = (b2 * (-c1) - b1 * (-c2)) / det;
                double y = (a1 * (-c2) - a2 * (-c1)) / det;

                intersection = Point2D(x, y);
                return true;
            }

            double Line::AngleWith(Line^ other)
            {
                if (other == nullptr) throw gcnew ArgumentNullException("other");
                // angle between normals -> same as angle between lines
                // using dot product of direction vectors: direction = ( -B, A ) for line Ax+By+C=0
                double ux = -this->B; double uy = this->A;
                double vx = -other->B; double vy = other->A;

                double dot = ux * vx + uy * vy;
                double nu = System::Math::Sqrt(ux * ux + uy * uy);
                double nv = System::Math::Sqrt(vx * vx + vy * vy);
                if (Utils::IsZero(nu) || Utils::IsZero(nv))
                    throw gcnew InvalidArgumentException("Invalid line direction.");

                double cosTheta = dot / (nu * nv);
                cosTheta = Utils::Clamp(cosTheta, -1.0, 1.0);
                double angle = System::Math::Acos(cosTheta);
                // return the acute angle between lines
                if (angle > (System::Math::PI / 2.0)) angle = System::Math::PI - angle;
                return angle;
            }

            String^ Line::ToString()
            {
                return String::Format("{0}x + {1}y + {2} = 0", A.ToString("G6"), B.ToString("G6"), C.ToString("G6"));
            }

            // -------------------------
            // Circle implementation
            // -------------------------
            Circle::Circle(Point2D center, double radius)
            {
                if (radius < 0.0) throw gcnew InvalidArgumentException("Radius must be non-negative.");
                this->center = center;
                this->radius = radius;
            }

            Point2D Circle::Center::get() { return center; }
            double Circle::Radius::get() { return radius; }

            bool Circle::Contains(Point2D p)
            {
                double d = Distance(center, p);
                return d <= radius + Constants::Epsilon;
            }

            double Circle::Area() { return Constants::Pi * radius * radius; }
            double Circle::Circumference() { return 2.0 * Constants::Pi * radius; }

            array<Point2D>^ Circle::IntersectWithLine(Line^ line)
            {
                if (line == nullptr) throw gcnew ArgumentNullException("line");

                // Project center onto line, get distance
                Point2D proj = line->ProjectPoint(center);
                double d = Distance(center, proj);

                if (d > radius + Constants::Epsilon) {
                    // no intersection
                    return gcnew array<Point2D>(0);
                }
                else if (Utils::AreEqual(d, radius, Constants::Tolerance)) {
                    // tangent: one point
                    array<Point2D>^ out = gcnew array<Point2D>(1);
                    out[0] = proj;
                    return out;
                }
                else {
                    // two intersection points: move from projection along line direction
                    // direction vector of the line (unit)
                    // line normal (A,B). direction = (-B, A)
                    double denom = System::Math::Sqrt(line->CoefA * line->CoefA + line->CoefB * line->CoefB);
                    double dx = -line->CoefB / denom;
                    double dy = line->CoefA / denom;

                    double h = System::Math::Sqrt(radius * radius - d * d); // distance from proj to intersection along dir

                    array<Point2D>^ out = gcnew array<Point2D>(2);
                    out[0] = Point2D(proj.X + dx * h, proj.Y + dy * h);
                    out[1] = Point2D(proj.X - dx * h, proj.Y - dy * h);
                    return out;
                }
            }

            array<Point2D>^ Circle::IntersectWithCircle(Circle^ other)
            {
                if (other == nullptr) throw gcnew ArgumentNullException("other");

                double d = Distance(this->center, other->center);

                // no intersection
                if (d > this->radius + other->radius + Constants::Epsilon) return gcnew array<Point2D>(0);
                // one inside another or coincident? check
                if (d < System::Math::Abs(this->radius - other->radius) - Constants::Epsilon) return gcnew array<Point2D>(0);
                // coincident circles (infinite intersections) - treat as none (user must handle)
                if (Utils::IsZero(d) && Utils::AreEqual(this->radius, other->radius)) return gcnew array<Point2D>(0);

                // tangent externally or internally -> one point
                if (Utils::AreEqual(d, this->radius + other->radius, Constants::Tolerance) ||
                    Utils::AreEqual(d, System::Math::Abs(this->radius - other->radius), Constants::Tolerance))
                {
                    // single point: compute along the line centers
                    double a = (this->radius * this->radius - other->radius * other->radius + d * d) / (2.0 * d);
                    double px = this->center.X + a * (other->center.X - this->center.X) / d;
                    double py = this->center.Y + a * (other->center.Y - this->center.Y) / d;
                    array<Point2D>^ out = gcnew array<Point2D>(1);
                    out[0] = Point2D(px, py);
                    return out;
                }

                // two intersection points
                double a_val = (this->radius * this->radius - other->radius * other->radius + d * d) / (2.0 * d);
                double h = System::Math::Sqrt(System::Math::Max(0.0, this->radius * this->radius - a_val * a_val));

                double xm = this->center.X + a_val * (other->center.X - this->center.X) / d;
                double ym = this->center.Y + a_val * (other->center.Y - this->center.Y) / d;

                double rx = -(other->center.Y - this->center.Y) * (h / d);
                double ry = (other->center.X - this->center.X) * (h / d);

                array<Point2D>^ out = gcnew array<Point2D>(2);
                out[0] = Point2D(xm + rx, ym + ry);
                out[1] = Point2D(xm - rx, ym - ry);
                return out;
            }

            String^ Circle::ToString()
            {
                return String::Format("Circle(center={0}, r={1})", center.ToString(), radius.ToString("G6"));
            }

        } // namespace Geometry
    } // namespace Math
} // namespace Agrine

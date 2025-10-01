#pragma once

// Analytic.h - Geometry/Analytic module
// Part of Agrine.Math (C++/CLI, .NET Framework 4.6.2)
// Provides analytic geometry primitives and operations (2D)

#include "Types.h"       // Point2D
#include "Constants.h"
#include "Utils.h"
#include "Exceptions.h"

using namespace System;
using namespace System::Collections::Generic;
using namespace Agrine::Math::Core;

namespace Agrine {
    namespace Math {
        namespace Geometry {

            // Simple helpers (free functions)
            double Distance(Point2D a, Point2D b);
            Point2D Midpoint(Point2D a, Point2D b);

            /// <summary>
            /// Represents a 2D line in form: A*x + B*y + C = 0
            /// Internally (A,B) normalized so that A^2 + B^2 = 1 (if possible).
            /// </summary>
            public ref class Line
            {
            private:
                double A;
                double B;
                double C;

                // normalize A,B,C such that sqrt(A^2+B^2) == 1 (unless A=B=0)
                void Normalize();

            public:
                // Constructors
                Line(double a, double b, double c);                 // explicit coefficients
                Line(Point2D p1, Point2D p2);                       // from two points (p1 != p2)
                static Line^ FromPointSlope(Point2D p, double slope); // y = slope * (x - x0) + y0

                // Accessors
                property double CoefA { double get(); }
                property double CoefB { double get(); }
                property double CoefC { double get(); }

                // Methods
                double DistanceToPoint(Point2D p);
                Point2D ProjectPoint(Point2D p);                     // orthogonal projection of p onto line
                bool Contains(Point2D p);                            // point exactly on line (within tolerance)

                // Intersection: returns true if intersection exists (not parallel) and out param is the point
                bool IntersectionWith(Line^ other, Point2D% intersection);

                // Angle between two lines (in radians, in [0, pi/2] by absolute value)
                double AngleWith(Line^ other);

                virtual String^ ToString() override;
            };

            /// <summary>
            /// Represents a 2D circle with center and radius.
            /// </summary>
            public ref class Circle
            {
            private:
                Point2D center;
                double radius;

            public:
                Circle(Point2D center, double radius);

                property Point2D Center { Point2D get(); }
                property double Radius { double get(); }

                bool Contains(Point2D p);   // contains (within tolerance)
                double Area();
                double Circumference();

                // Intersections:
                // Intersect line -> return array of 0,1,2 Point2D
                array<Point2D>^ IntersectWithLine(Line^ line);

                // Intersect circle -> return array of 0,1,2 Point2D
                array<Point2D>^ IntersectWithCircle(Circle^ other);

                virtual String^ ToString() override;
            };

        } // namespace Geometry
    } // namespace Math
} // namespace Agrine

#pragma once

// Computational.h - Geometry / Computational geometry algorithms
// Part of Agrine.Math (C++/CLI, .NET Framework 4.6.2)

#include "Types.h"       // Point2D
#include "Constants.h"
#include "Utils.h"
#include "Exceptions.h"
#include "Euclidean.h"

using namespace System;
using namespace System::Collections::Generic;
using namespace Agrine::Math::Core;
using namespace Agrine::Math::Geometry;

namespace Agrine {
    namespace Math {
        namespace Geometry {
            /// <summary>
            /// Computational geometry algorithms (2D).
            /// </summary>
            public ref class Computational abstract sealed
            {
            public:
                // Convex hull (Andrew monotone chain)
                // Input: array of points (may contain duplicates). Output: convex hull in CCW order (no duplicate endpoint).
                static array<Point2D>^ ConvexHull(array<Point2D>^ points);

                // Segment intersection: returns array of 0,1,2 points
                // - empty => no intersection
                // - length==1 => single intersection point (proper or tangent)
                // - length==2 => overlapping segment interval (endpoints of overlap)
                static array<Point2D>^ IntersectSegments(Point2D p1, Point2D p2, Point2D q1, Point2D q2);

                // Closest pair of points (returns minimal distance)
                // If fewer than 2 points -> throws InvalidArgumentException
                static double ClosestPairDistance(array<Point2D>^ points);

                // Diameter (farthest pair) using rotating calipers; returns pair as array of 2 points (p1,p2)
                // If fewer than 2 points -> empty array
                static array<Point2D>^ Diameter(array<Point2D>^ points);

            private:
                // helpers
                static int CompareByXThenY(Point2D a, Point2D b);
                static double Cross(const Point2D& O, const Point2D& A, const Point2D& B);
                static bool OnSegment(Point2D a, Point2D b, Point2D p);
                static int Sign(double x);
                static array<Point2D>^ UniqueSorted(array<Point2D>^ pts);
                // closest-pair internals
                static double ClosestPairRecursive(List<Point2D>^ ptsSortedX, int l, int r, List<Point2D>^ aux);
            };
        }
    }
}

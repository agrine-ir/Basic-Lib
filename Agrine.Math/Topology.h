#pragma once

// Topology.h - simple, robust topology utilities for 2D point-clouds
// Part of Agrine.Math (C++/CLI, .NET Framework 4.6.2)

#include "Types.h"       // Point2D
#include "Constants.h"
#include "Utils.h"
#include "Exceptions.h"
#include "Euclidean.h" // Distance2D, etc.

using namespace System;
using namespace System::Collections::Generic;
using namespace Agrine::Math::Core;
using namespace Agrine::Math::Geometry;

namespace Agrine {
    namespace Math {
        namespace Topology {

            /// <summary>
            /// TopologyTools - working on finite point sets in R^2 (point-cloud).
            /// All neighborhood ops are relative to an ambient set provided by the caller.
            /// </summary>
            public ref class TopologyTools abstract sealed
            {
            public:
                // Neighborhoods & set operators (2D)
                static array<Point2D>^ EpsilonNeighborhood2D(array<Point2D>^ S, array<Point2D>^ ambient, double eps);
                static array<Point2D>^ Closure2D(array<Point2D>^ S, array<Point2D>^ ambient, double eps);
                static array<Point2D>^ Interior2D(array<Point2D>^ S, array<Point2D>^ ambient, double eps);
                static array<Point2D>^ Boundary2D(array<Point2D>^ S, array<Point2D>^ ambient, double eps);

                static bool IsOpen2D(array<Point2D>^ S, array<Point2D>^ ambient, double eps);
                static bool IsClosed2D(array<Point2D>^ S, array<Point2D>^ ambient, double eps);

                // Point classification
                static bool IsLimitPoint2D(Point2D p, array<Point2D>^ S, double eps);
                static bool IsIsolatedPoint2D(Point2D p, array<Point2D>^ S, double eps);

                // Connectivity
                static array<array<Point2D>^>^ ConnectedComponents2D(array<Point2D>^ S, double eps);
                static bool IsConnected2D(array<Point2D>^ S, double eps);

                // Distances
                static double HausdorffDistance2D(array<Point2D>^ A, array<Point2D>^ B);

                // Rips 0D persistence (simple)
                // returns intervals (birth, death) as KeyValuePair<double,double>
                static array<KeyValuePair<double, double>>^ RipsPersistence0D(array<Point2D>^ points);

            private:
                // internal helpers (Point2D-specific)
                static array<Point2D>^ UniqueArray(array<Point2D>^ arr);
                static bool ContainsPoint(array<Point2D>^ arr, Point2D p);
                static int CompareEdgeTuples(Tuple<int, int, double>^ a, Tuple<int, int, double>^ b);
            };

        } // namespace Topology
    } // namespace Math
} // namespace Agrine

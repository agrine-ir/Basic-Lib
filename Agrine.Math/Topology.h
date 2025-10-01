#pragma once

// Topology.h - Geometry / Topology utilities (point-cloud based, 2D/3D)
// Part of Agrine.Math (C++/CLI, .NET Framework 4.6.2)

#include "Types.h"       // Point2D, Point3D
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
        namespace Topology {

            /// <summary>
            /// TopologyTools - utilities for working with point-clouds as metric-space approximations.
            /// All methods that require neighborhoods use a radius 'eps'.
            /// </summary>
            public ref class TopologyTools abstract sealed
            {
            public:
                // ===== Neighborhoods & set operators (2D) =====
                // Epsilon-neighborhood of a set S inside ambient (points considered as universe)
                static array<Point2D>^ EpsilonNeighborhood2D(array<Point2D>^ S, array<Point2D>^ ambient, double eps);

                // Closure, Interior, Boundary relative to ambient point cloud
                static array<Point2D>^ Closure2D(array<Point2D>^ S, array<Point2D>^ ambient, double eps);
                static array<Point2D>^ Interior2D(array<Point2D>^ S, array<Point2D>^ ambient, double eps);
                static array<Point2D>^ Boundary2D(array<Point2D>^ S, array<Point2D>^ ambient, double eps);

                // Open/Closed (relative to ambient and eps)
                static bool IsOpen2D(array<Point2D>^ S, array<Point2D>^ ambient, double eps);
                static bool IsClosed2D(array<Point2D>^ S, array<Point2D>^ ambient, double eps);

                // Point classification
                static bool IsLimitPoint2D(Point2D p, array<Point2D>^ S, double eps);
                static bool IsIsolatedPoint2D(Point2D p, array<Point2D>^ S, double eps);

                // ===== Connectivity (2D) =====
                // Connected components of set S using eps-neighborhood graph (edge if dist <= eps)
                static array<array<Point2D>^>^ ConnectedComponents2D(array<Point2D>^ S, double eps);
                static bool IsConnected2D(array<Point2D>^ S, double eps);

                // ===== Distances =====
                // Hausdorff distance between two finite point sets (2D)
                static double HausdorffDistance2D(array<Point2D>^ A, array<Point2D>^ B);

                // ===== Persistence (0D Vietoris-Rips) =====
                // Returns list of (birth, death) intervals as KeyValuePair<double,double>.
                // For 0D persistence births are 0; death = edge distance when components merge; surviving component death = +inf.
                static array<KeyValuePair<double, double>>^ RipsPersistence0D(array<Point2D>^ points);

                // ===== 3D versions (same semantics) =====
                static array<Point3D>^ EpsilonNeighborhood3D(array<Point3D>^ S, array<Point3D>^ ambient, double eps);
                static array<Point3D>^ Closure3D(array<Point3D>^ S, array<Point3D>^ ambient, double eps);
                static array<Point3D>^ Interior3D(array<Point3D>^ S, array<Point3D>^ ambient, double eps);
                static array<Point3D>^ Boundary3D(array<Point3D>^ S, array<Point3D>^ ambient, double eps);
                static array<array<Point3D>^>^ ConnectedComponents3D(array<Point3D>^ S, double eps);
                static bool IsConnected3D(array<Point3D>^ S, double eps);
                static double HausdorffDistance3D(array<Point3D>^ A, array<Point3D>^ B);

            private:
                // small helpers
                static array<T>^ UniqueArray<T>(array<T>^ arr);
                generic<typename T>
                where T : value class
                    static List<T>^ ToList(array<T>^ arr);
            };
        }
    }
}

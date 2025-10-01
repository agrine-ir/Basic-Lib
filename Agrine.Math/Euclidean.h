#pragma once

// Euclidean.h - Geometry / Euclidean utilities (2D & 3D)
// Part of Agrine.Math (C++/CLI, .NET Framework 4.6.2)

#include "Types.h"        // Point2D, Point3D
#include "Constants.h"
#include "Utils.h"
#include "Exceptions.h"

using namespace System;
using namespace System::Collections::Generic;
using namespace Agrine::Math::Core;

namespace Agrine {
	namespace Math {
		namespace Geometry {

			// -----------------------
			// 2D helpers
			// -----------------------
			double Distance2D(Point2D a, Point2D b);
			double Distance3D(Point3D a, Point3D b);

			// point to segment (2D) - shortest distance
			double DistancePointToSegment2D(Point2D p, Point2D a, Point2D b);
			Point2D ClosestPointOnSegment2D(Point2D p, Point2D a, Point2D b);

			// triangle area (2D) using shoelace (signed area)
			double TriangleArea2D(Point2D a, Point2D b, Point2D c);
			double PolygonArea2D(array<Point2D>^ polygon); // polygon vertices in order (CW or CCW)
			Point2D PolygonCentroid2D(array<Point2D>^ polygon); // NOTE: placeholder to be replaced

			// point in polygon: ray-casting winding rule (true if inside or on edge)
			bool PointInPolygon(array<Point2D>^ polygon, Point2D p);

			// -----------------------
			// 3D helpers
			// -----------------------
			double TriangleArea3D(Point3D a, Point3D b, Point3D c);
			double TetrahedronVolume(Point3D a, Point3D b, Point3D c, Point3D d);

			// Plane class: Ax + By + Cz + D = 0 (A,B,C normalized)
			public ref class Plane
			{
			private:
				double A, B, C, D;
				void Normalize();

			public:
				// constructors
				Plane(double a, double b, double c, double d);
				Plane(Point3D p1, Point3D p2, Point3D p3); // through three non-collinear points

				// accessors
				property double CoefA { double get(); }
				property double CoefB { double get(); }
				property double CoefC { double get(); }
				property double CoefD { double get(); }

				// operations
				double DistanceToPoint(Point3D p);
				Point3D ProjectPointToPlane(Point3D p);
				bool Contains(Point3D p); // within tolerance

				virtual String^ ToString() override;
			};

			// point to plane distance / projection convenience functions
			double DistancePointToPlane(Point3D p, Plane^ plane);
			Point3D ClosestPointOnPlane(Point3D p, Plane^ plane);

		} // namespace Geometry
	} // namespace Math
}

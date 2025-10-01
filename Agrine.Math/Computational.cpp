#include "Computational.h"

using namespace System;
using namespace System::Collections::Generic;
using namespace Agrine::Math::Core;
using namespace Agrine::Math::Geometry;

namespace Agrine {
    namespace Math {
        namespace Geometry {

            // ---------------------------------------
            // Helpers
            // ---------------------------------------
            // Compare by X then Y
            int Computational::CompareByXThenY(Point2D a, Point2D b)
            {
                if (!Utils::AreEqual(a.X, b.X))
                    return (a.X < b.X) ? -1 : 1;
                if (!Utils::AreEqual(a.Y, b.Y))
                    return (a.Y < b.Y) ? -1 : 1;
                return 0;
            }

            // Cross product (OA x OB) where O,A,B are points
            double Computational::Cross(const Point2D& O, const Point2D& A, const Point2D& B)
            {
                double ax = A.X - O.X;
                double ay = A.Y - O.Y;
                double bx = B.X - O.X;
                double by = B.Y - O.Y;
                return ax * by - ay * bx;
            }

            // On-segment: check if p lies on segment ab (assumes collinear)
            bool Computational::OnSegment(Point2D a, Point2D b, Point2D p)
            {
                double minx = System::Math::Min(a.X, b.X) - Constants::Epsilon;
                double maxx = System::Math::Max(a.X, b.X) + Constants::Epsilon;
                double miny = System::Math::Min(a.Y, b.Y) - Constants::Epsilon;
                double maxy = System::Math::Max(a.Y, b.Y) + Constants::Epsilon;
                return p.X >= minx && p.X <= maxx && p.Y >= miny && p.Y <= maxy;
            }

            int Computational::Sign(double x)
            {
                if (Utils::IsZero(x)) return 0;
                return (x > 0.0) ? 1 : -1;
            }

            // Remove duplicates and sort points by X then Y
            array<Point2D>^ Computational::UniqueSorted(array<Point2D>^ pts)
            {
                if (pts == nullptr) throw gcnew ArgumentNullException("pts");
                int n = pts->Length;
                if (n == 0) return gcnew array<Point2D>(0);

                // Convert to List for sort
                List<Point2D>^ list = gcnew List<Point2D>(pts);
                // sort using comparison delegate
                list->Sort(gcnew Comparison<Point2D>(&Computational::CompareByXThenY));

                // unique
                List<Point2D>^ uniq = gcnew List<Point2D>();
                Point2D last = list[0];
                uniq->Add(last);
                for (int i = 1; i < list->Count; i++) {
                    Point2D cur = list[i];
                    if (CompareByXThenY(cur, last) != 0) {
                        uniq->Add(cur);
                        last = cur;
                    }
                }

                return uniq->ToArray();
            }

            // ---------------------------------------
            // Convex Hull - Andrew Monotone Chain
            // ---------------------------------------
            array<Point2D>^ Computational::ConvexHull(array<Point2D>^ points)
            {
                if (points == nullptr) throw gcnew ArgumentNullException("points");
                int n = points->Length;
                if (n == 0) return gcnew array<Point2D>(0);
                if (n == 1) return gcnew array<Point2D>(1) { points[0] };

                array<Point2D>^ pts = UniqueSorted(points);
                int m = pts->Length;
                if (m == 1) return gcnew array<Point2D>(1) { pts[0] };

                List<Point2D>^ lower = gcnew List<Point2D>();
                for (int i = 0; i < m; i++) {
                    while (lower->Count >= 2 && Cross(lower[lower->Count - 2], lower[lower->Count - 1], pts[i]) <= 0.0)
                        lower->RemoveAt(lower->Count - 1);
                    lower->Add(pts[i]);
                }

                List<Point2D>^ upper = gcnew List<Point2D>();
                for (int i = m - 1; i >= 0; i--) {
                    while (upper->Count >= 2 && Cross(upper[upper->Count - 2], upper[upper->Count - 1], pts[i]) <= 0.0)
                        upper->RemoveAt(upper->Count - 1);
                    upper->Add(pts[i]);
                }

                // Concatenate lower and upper excluding last point of each (it's the start of the other)
                List<Point2D>^ hull = gcnew List<Point2D>(lower);
                // skip first and last of upper? standard: skip first and last of upper because they are endpoints
                for (int i = 1; i + 1 < upper->Count; i++)
                    hull->Add(upper[i]);

                return hull->ToArray();
            }

            // ---------------------------------------
            // Segment intersection: robust handling
            // ---------------------------------------
            array<Point2D>^ Computational::IntersectSegments(Point2D p1, Point2D p2, Point2D q1, Point2D q2)
            {
                // compute orientations
                double d1 = Cross(p1, p2, q1);
                double d2 = Cross(p1, p2, q2);
                double d3 = Cross(q1, q2, p1);
                double d4 = Cross(q1, q2, p2);

                bool col12 = Utils::IsZero(d1) && Utils::IsZero(d2) && Utils::IsZero(d3) && Utils::IsZero(d4);
                List<Point2D>^ out = gcnew List<Point2D>();

                // General proper intersection
                if (Sign(d1) * Sign(d2) < 0 && Sign(d3) * Sign(d4) < 0) {
                    // compute intersection point of two lines
                    // line p: p1 + t*(p2-p1), line q: q1 + u*(q2-q1)
                    double A1 = p2.Y - p1.Y;
                    double B1 = p1.X - p2.X;
                    double C1 = A1 * p1.X + B1 * p1.Y;

                    double A2 = q2.Y - q1.Y;
                    double B2 = q1.X - q2.X;
                    double C2 = A2 * q1.X + B2 * q1.Y;

                    double det = A1 * B2 - A2 * B1;
                    if (Utils::IsZero(det)) {
                        // numerical degenerate: return empty
                        return gcnew array<Point2D>(0);
                    }
                    double x = (B2 * C1 - B1 * C2) / det;
                    double y = (A1 * C2 - A2 * C1) / det;
                    out->Add(Point2D(x, y));
                    return out->ToArray();
                }

                // Special cases: endpoints on segments
                if (Utils::IsZero(d1) && OnSegment(p1, p2, q1)) out->Add(q1);
                if (Utils::IsZero(d2) && OnSegment(p1, p2, q2)) out->Add(q2);
                if (Utils::IsZero(d3) && OnSegment(q1, q2, p1)) out->Add(p1);
                if (Utils::IsZero(d4) && OnSegment(q1, q2, p2)) out->Add(p2);

                // If multiple points added due to collinearity, we must return the overlap segment endpoints (unique)
                if (out->Count == 0) return gcnew array<Point2D>(0);

                // make unique
                List<Point2D>^ uniq = gcnew List<Point2D>();
                for each (Point2D pt in out) {
                    bool found = false;
                    for each (Point2D up in uniq) {
                        if (Utils::AreEqual(pt.X, up.X) && Utils::AreEqual(pt.Y, up.Y)) { found = true; break; }
                    }
                    if (!found) uniq->Add(pt);
                }

                // If collinear overlap (we may have 2 endpoints), sort them along the segment and return endpoints of overlap
                if (uniq->Count > 1) {
                    // sort by projection along p1->p2
                    uniq->Sort(gcnew Comparison<Point2D>([&](Point2D a, Point2D b) {
                        double ta = ((a.X - p1.X) * (p2.X - p1.X) + (a.Y - p1.Y) * (p2.Y - p1.Y));
                        double tb = ((b.X - p1.X) * (p2.X - p1.X) + (b.Y - p1.Y) * (p2.Y - p1.Y));
                        if (!Utils::AreEqual(ta, tb)) return (ta < tb) ? -1 : 1;
                        return 0;
                        }));
                    return uniq->ToArray();
                }

                // single point
                return uniq->ToArray();
            }

            // ---------------------------------------
            // Closest pair - divide & conquer (O(n log n))
            // ---------------------------------------
            double Computational::ClosestPairDistance(array<Point2D>^ points)
            {
                if (points == nullptr) throw gcnew ArgumentNullException("points");
                int n = points->Length;
                if (n < 2) throw gcnew InvalidArgumentException("At least two points are required.");

                // prepare sorted by X list
                List<Point2D>^ pts = gcnew List<Point2D>(points);
                pts->Sort(gcnew Comparison<Point2D>(&Computational::CompareByXThenY));

                // aux buffer for merge by Y
                List<Point2D>^ aux = gcnew List<Point2D>(n);
                for (int i = 0; i < n; i++) aux->Add(pts[i]);

                double d = ClosestPairRecursive(pts, 0, n - 1, aux);
                return d;
            }

            // recursive: ptsSortedX is sorted by X (in List)
            double Computational::ClosestPairRecursive(List<Point2D>^ ptsSortedX, int l, int r, List<Point2D>^ aux)
            {
                int len = r - l + 1;
                if (len <= 3) {
                    // brute force
                    double best = System::Double::PositiveInfinity;
                    for (int i = l; i <= r; i++)
                        for (int j = i + 1; j <= r; j++) {
                            double dx = ptsSortedX[i].X - ptsSortedX[j].X;
                            double dy = ptsSortedX[i].Y - ptsSortedX[j].Y;
                            double dist = System::Math::Sqrt(dx * dx + dy * dy);
                            if (dist < best) best = dist;
                        }
                    // sort by Y for parent convenience (simple insertion sort)
                    for (int i = l + 1; i <= r; i++) {
                        Point2D key = ptsSortedX[i];
                        int j = i - 1;
                        while (j >= l && ptsSortedX[j].Y > key.Y) {
                            ptsSortedX[j + 1] = ptsSortedX[j];
                            j--;
                        }
                        ptsSortedX[j + 1] = key;
                    }
                    return best;
                }

                int mid = (l + r) / 2;
                double midx = ptsSortedX[mid].X;

                double dl = ClosestPairRecursive(ptsSortedX, l, mid, aux);
                double dr = ClosestPairRecursive(ptsSortedX, mid + 1, r, aux);
                double d = System::Math::Min(dl, dr);

                // merge sort by Y into aux
                int i = l, j = mid + 1, k = l;
                while (i <= mid && j <= r) {
                    if (ptsSortedX[i].Y < ptsSortedX[j].Y) aux[k++] = ptsSortedX[i++];
                    else aux[k++] = ptsSortedX[j++];
                }
                while (i <= mid) aux[k++] = ptsSortedX[i++];
                while (j <= r) aux[k++] = ptsSortedX[j++];
                for (int s = l; s <= r; s++) ptsSortedX[s] = aux[s];

                // build strip of points within d of midx
                List<Point2D>^ strip = gcnew List<Point2D>();
                for (int s = l; s <= r; s++) {
                    if (System::Math::Abs(ptsSortedX[s].X - midx) < d + Constants::Epsilon)
                        strip->Add(ptsSortedX[s]);
                }

                // check up to next 7 points in strip
                for (int s = 0; s < strip->Count; s++) {
                    for (int t = s + 1; t < strip->Count && (strip[t].Y - strip[s].Y) < d + Constants::Epsilon; t++) {
                        double dx = strip[s].X - strip[t].X;
                        double dy = strip[s].Y - strip[t].Y;
                        double dist = System::Math::Sqrt(dx * dx + dy * dy);
                        if (dist < d) d = dist;
                    }
                }

                return d;
            }

            // ---------------------------------------
            // Diameter (farthest pair) via rotating calipers
            // ---------------------------------------
            array<Point2D>^ Computational::Diameter(array<Point2D>^ points)
            {
                if (points == nullptr) throw gcnew ArgumentNullException("points");
                int n = points->Length;
                if (n < 2) return gcnew array<Point2D>(0);

                array<Point2D>^ hull = ConvexHull(points);
                int m = hull->Length;
                if (m == 1) return gcnew array<Point2D>(1) { hull[0] };
                if (m == 2) return gcnew array<Point2D>(2) { hull[0], hull[1] };

                // find farthest pair using rotating calipers
                int k = 1;
                // find k so that area between 0-1-k is maximized
                while (System::Math::Abs(Cross(hull[0], hull[1], hull[(k + 1) % m])) > System::Math::Abs(Cross(hull[0], hull[1], hull[k]))) k++;

                double bestDist = 0.0;
                Point2D bestA = hull[0], bestB = hull[0];

                int i = 0, j = k;
                while (i <= k && j < m) {
                    // update best
                    double dx = hull[i].X - hull[j].X;
                    double dy = hull[i].Y - hull[j].Y;
                    double d = dx * dx + dy * dy;
                    if (d > bestDist) {
                        bestDist = d;
                        bestA = hull[i];
                        bestB = hull[j];
                    }

                    // advance j while area increases
                    int ni = (i + 1) % m;
                    int nj = (j + 1) % m;
                    double crossVal = Cross(hull[i], hull[ni], hull[nj]) - Cross(hull[i], hull[ni], hull[j]);
                    if (crossVal > 0) j = nj;
                    else i = ni;
                    if (i == 0) break; // one full loop
                }

                return gcnew array<Point2D>(2) { bestA, bestB };
            }

        } // namespace Geometry
    } // namespace Math
}

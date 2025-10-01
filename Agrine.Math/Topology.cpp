#include "Topology.h"

using namespace System;
using namespace System::Collections::Generic;
using namespace Agrine::Math::Core;
using namespace Agrine::Math::Geometry;

namespace Agrine {
    namespace Math {
        namespace Topology {

            // -------------------------
            // Helpers (Point2D-specific)
            // -------------------------
            array<Point2D>^ TopologyTools::UniqueArray(array<Point2D>^ arr)
            {
                if (arr == nullptr) return gcnew array<Point2D>(0);
                List<Point2D>^ uniq = gcnew List<Point2D>();
                for each (Point2D p in arr) {
                    bool found = false;
                    for each (Point2D q in uniq) {
                        if (Utils::AreEqual(p.X, q.X) && Utils::AreEqual(p.Y, q.Y)) { found = true; break; }
                    }
                    if (!found) uniq->Add(p);
                }
                return uniq->ToArray();
            }

            bool TopologyTools::ContainsPoint(array<Point2D>^ arr, Point2D p)
            {
                if (arr == nullptr) return false;
                for each (Point2D q in arr) {
                    if (Utils::AreEqual(p.X, q.X) && Utils::AreEqual(p.Y, q.Y)) return true;
                }
                return false;
            }

            int TopologyTools::CompareEdgeTuples(Tuple<int, int, double>^ a, Tuple<int, int, double>^ b)
            {
                if (a->Item3 < b->Item3) return -1;
                if (a->Item3 > b->Item3) return 1;
                return 0;
            }

            // -------------------------
            // Neighborhoods & set ops
            // -------------------------
            array<Point2D>^ TopologyTools::EpsilonNeighborhood2D(array<Point2D>^ S, array<Point2D>^ ambient, double eps)
            {
                if (S == nullptr) throw gcnew ArgumentNullException("S");
                if (ambient == nullptr) throw gcnew ArgumentNullException("ambient");
                if (eps < 0.0) throw gcnew InvalidArgumentException("eps must be non-negative.");

                List<Point2D>^ out = gcnew List<Point2D>();
                for each (Point2D a in ambient) {
                    for each (Point2D s in S) {
                        if (Distance2D(a, s) <= eps + Constants::Epsilon) { out->Add(a); break; }
                    }
                }
                return UniqueArray(out->ToArray());
            }

            array<Point2D>^ TopologyTools::Closure2D(array<Point2D>^ S, array<Point2D>^ ambient, double eps)
            {
                if (S == nullptr) throw gcnew ArgumentNullException("S");
                if (ambient == nullptr) throw gcnew ArgumentNullException("ambient");
                if (eps < 0.0) throw gcnew InvalidArgumentException("eps must be non-negative.");

                List<Point2D>^ out = gcnew List<Point2D>();
                // include original S
                for each (Point2D p in S) out->Add(p);
                // add ambient points within eps of any s
                for each (Point2D a in ambient) {
                    for each (Point2D s in S) {
                        if (Distance2D(a, s) <= eps + Constants::Epsilon) { out->Add(a); break; }
                    }
                }
                return UniqueArray(out->ToArray());
            }

            array<Point2D>^ TopologyTools::Interior2D(array<Point2D>^ S, array<Point2D>^ ambient, double eps)
            {
                if (S == nullptr) throw gcnew ArgumentNullException("S");
                if (ambient == nullptr) throw gcnew ArgumentNullException("ambient");
                if (eps < 0.0) throw gcnew InvalidArgumentException("eps must be non-negative.");

                List<Point2D>^ out = gcnew List<Point2D>();
                // For each p in S: all ambient points within eps must belong to S
                for each (Point2D p in S) {
                    bool isInterior = true;
                    for each (Point2D a in ambient) {
                        if (Distance2D(p, a) <= eps + Constants::Epsilon) {
                            if (!ContainsPoint(S, a)) { isInterior = false; break; }
                        }
                    }
                    if (isInterior) out->Add(p);
                }
                return UniqueArray(out->ToArray());
            }

            array<Point2D>^ TopologyTools::Boundary2D(array<Point2D>^ S, array<Point2D>^ ambient, double eps)
            {
                array<Point2D>^ cl = Closure2D(S, ambient, eps);
                array<Point2D>^ inr = Interior2D(S, ambient, eps);

                List<Point2D>^ out = gcnew List<Point2D>();
                for each (Point2D p in cl) {
                    if (!ContainsPoint(inr, p)) out->Add(p);
                }
                return UniqueArray(out->ToArray());
            }

            bool TopologyTools::IsOpen2D(array<Point2D>^ S, array<Point2D>^ ambient, double eps)
            {
                array<Point2D>^ inr = Interior2D(S, ambient, eps);
                // compare sets by membership
                List<Point2D>^ sL = gcnew List<Point2D>(S);
                List<Point2D>^ inL = gcnew List<Point2D>(inr);
                if (sL->Count != inL->Count) return false;
                for each (Point2D p in sL) if (!ContainsPoint(inr, p)) return false;
                return true;
            }

            bool TopologyTools::IsClosed2D(array<Point2D>^ S, array<Point2D>^ ambient, double eps)
            {
                array<Point2D>^ cl = Closure2D(S, ambient, eps);
                List<Point2D>^ sL = gcnew List<Point2D>(S);
                List<Point2D>^ clL = gcnew List<Point2D>(cl);
                if (sL->Count != clL->Count) return false;
                for each (Point2D p in sL) if (!ContainsPoint(cl, p)) return false;
                return true;
            }

            // -------------------------
            // Point classification
            // -------------------------
            bool TopologyTools::IsLimitPoint2D(Point2D p, array<Point2D>^ S, double eps)
            {
                if (S == nullptr) throw gcnew ArgumentNullException("S");
                if (eps < 0.0) throw gcnew InvalidArgumentException("eps must be non-negative.");
                for each (Point2D q in S) {
                    if (Utils::AreEqual(p.X, q.X) && Utils::AreEqual(p.Y, q.Y)) continue;
                    if (Distance2D(p, q) <= eps + Constants::Epsilon) return true;
                }
                return false;
            }

            bool TopologyTools::IsIsolatedPoint2D(Point2D p, array<Point2D>^ S, double eps)
            {
                if (S == nullptr) throw gcnew ArgumentNullException("S");
                if (eps < 0.0) throw gcnew InvalidArgumentException("eps must be non-negative.");
                for each (Point2D q in S) {
                    if (Utils::AreEqual(p.X, q.X) && Utils::AreEqual(p.Y, q.Y)) continue;
                    if (Distance2D(p, q) <= eps + Constants::Epsilon) return false;
                }
                return true;
            }

            // -------------------------
            // Connectivity (BFS on eps-graph)
            // -------------------------
            array<array<Point2D>^>^ TopologyTools::ConnectedComponents2D(array<Point2D>^ S, double eps)
            {
                if (S == nullptr) throw gcnew ArgumentNullException("S");
                if (eps < 0.0) throw gcnew InvalidArgumentException("eps must be non-negative.");

                int n = S->Length;
                List<array<Point2D>^>^ result = gcnew List<array<Point2D>^>();
                if (n == 0) return result->ToArray();

                array<bool>^ seen = gcnew array<bool>(n);
                for (int i = 0; i < n; ++i) seen[i] = false;

                for (int i = 0; i < n; ++i) {
                    if (seen[i]) continue;
                    List<int>^ q = gcnew List<int>();
                    q->Add(i); seen[i] = true;
                    List<Point2D>^ comp = gcnew List<Point2D>();
                    for (int qi = 0; qi < q->Count; ++qi) {
                        int idx = q[qi];
                        comp->Add(S[idx]);
                        for (int j = 0; j < n; ++j) {
                            if (seen[j]) continue;
                            if (Distance2D(S[idx], S[j]) <= eps + Constants::Epsilon) {
                                seen[j] = true;
                                q->Add(j);
                            }
                        }
                    }
                    result->Add(comp->ToArray());
                }

                return result->ToArray();
            }

            bool TopologyTools::IsConnected2D(array<Point2D>^ S, double eps)
            {
                array<array<Point2D>^>^ comps = ConnectedComponents2D(S, eps);
                return (comps->Length <= 1);
            }

            // -------------------------
            // Hausdorff distance
            // -------------------------
            static double onewayHausdorff2D(array<Point2D>^ A, array<Point2D>^ B)
            {
                double best = 0.0;
                if (A == nullptr || B == nullptr) return System::Double::PositiveInfinity;
                if (A->Length == 0) return System::Double::PositiveInfinity;
                for each (Point2D a in A) {
                    double minD = System::Double::PositiveInfinity;
                    for each (Point2D b in B) {
                        double d = Distance2D(a, b);
                        if (d < minD) minD = d;
                    }
                    if (minD > best) best = minD;
                }
                return best;
            }

            double TopologyTools::HausdorffDistance2D(array<Point2D>^ A, array<Point2D>^ B)
            {
                if (A == nullptr) throw gcnew ArgumentNullException("A");
                if (B == nullptr) throw gcnew ArgumentNullException("B");
                if (A->Length == 0 && B->Length == 0) return 0.0;
                if (A->Length == 0 || B->Length == 0) return System::Double::PositiveInfinity;

                double ab = onewayHausdorff2D(A, B);
                double ba = onewayHausdorff2D(B, A);
                return System::Math::Max(ab, ba);
            }

            // -------------------------
            // Rips 0D persistence (simple)
            // -------------------------
            ref class UF {
            public:
                array<int>^ parent;
                array<int>^ rank;
                UF(int n) {
                    parent = gcnew array<int>(n);
                    rank = gcnew array<int>(n);
                    for (int i = 0; i < n; ++i) { parent[i] = i; rank[i] = 0; }
                }
                int Find(int x) {
                    if (parent[x] != x) parent[x] = Find(parent[x]);
                    return parent[x];
                }
                bool Union(int x, int y) {
                    int rx = Find(x), ry = Find(y);
                    if (rx == ry) return false;
                    if (rank[rx] < rank[ry]) parent[rx] = ry;
                    else if (rank[ry] < rank[rx]) parent[ry] = rx;
                    else { parent[ry] = rx; rank[rx]++; }
                    return true;
                }
            };

            array<KeyValuePair<double, double>>^ TopologyTools::RipsPersistence0D(array<Point2D>^ points)
            {
                if (points == nullptr) throw gcnew ArgumentNullException("points");
                int n = points->Length;
                List<KeyValuePair<double, double>>^ intervals = gcnew List<KeyValuePair<double, double>>();

                if (n == 0) return intervals->ToArray();
                if (n == 1) {
                    intervals->Add(KeyValuePair<double, double>(0.0, System::Double::PositiveInfinity));
                    return intervals->ToArray();
                }

                // build edges list: (i,j,dist) with i<j
                List<Tuple<int, int, double>^>^ edges = gcnew List<Tuple<int, int, double>^>();
                for (int i = 0; i < n; ++i) {
                    for (int j = i + 1; j < n; ++j) {
                        double d = Distance2D(points[i], points[j]);
                        edges->Add(gcnew Tuple<int, int, double>(i, j, d));
                    }
                }

                // sort edges by distance ascending using static comparer
                edges->Sort(gcnew Comparison<Tuple<int, int, double>^>(&TopologyTools::CompareEdgeTuples));

                // birth = 0 for all, death = +inf initially
                array<double>^ birth = gcnew array<double>(n);
                array<double>^ death = gcnew array<double>(n);
                for (int i = 0; i < n; ++i) { birth[i] = 0.0; death[i] = System::Double::PositiveInfinity; }

                UF uf(n);
                for each (Tuple<int, int, double> ^ e in edges) {
                    int i = e->Item1, j = e->Item2;
                    double d = e->Item3;
                    int ri = uf.Find(i), rj = uf.Find(j);
                    if (ri == rj) continue;
                    // choose which component 'dies' at this distance - choose higher index rep for determinism
                    int toDie = (ri > rj) ? ri : rj;
                    death[toDie] = d;
                    uf.Union(ri, rj);
                }

                for (int i = 0; i < n; ++i) intervals->Add(KeyValuePair<double, double>(birth[i], death[i]));
                return intervals->ToArray();
            }

        } // namespace Topology
    } // namespace Math
} // namespace Agrine

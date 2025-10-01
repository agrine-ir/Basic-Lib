#include "Topology.h"

using namespace System;
using namespace System::Collections::Generic;
using namespace Agrine::Math::Core;
using namespace Agrine::Math::Geometry;

namespace Agrine {
    namespace Math {
        namespace Topology {

            // -----------------------
            // Helpers (generic)
            // -----------------------
            generic<typename T>
            where T : value class
                List<T>^ TopologyTools::ToList(array<T>^ arr)
            {
                List<T>^ L = gcnew List<T>();
                if (arr == nullptr) return L;
                for each (T v in arr) L->Add(v);
                return L;
            }

            generic<typename T>
            where T : value class
                array<T>^ TopologyTools::UniqueArray(array<T>^ arr)
            {
                if (arr == nullptr) return gcnew array<T>(0);
                List<T>^ lst = gcnew List<T>(arr);
                List<T>^ uniq = gcnew List<T>();
                for each (T p in lst) {
                    bool f = false;
                    for each (T q in uniq) {
                        // Compare component-wise using Utils where possible via dynamic cast? For simplicity rely on equality
                        if (p.Equals(q)) { f = true; break; }
                    }
                    if (!f) uniq->Add(p);
                }
                return uniq->ToArray();
            }

            // -----------------------
            // 2D: Neighborhood helpers
            // -----------------------
            array<Point2D>^ TopologyTools::EpsilonNeighborhood2D(array<Point2D>^ S, array<Point2D>^ ambient, double eps)
            {
                if (S == nullptr) throw gcnew ArgumentNullException("S");
                if (ambient == nullptr) throw gcnew ArgumentNullException("ambient");
                if (eps < 0) throw gcnew InvalidArgumentException("eps must be non-negative.");

                List<Point2D>^ out = gcnew List<Point2D>();
                for each (Point2D a in ambient) {
                    for each (Point2D s in S) {
                        if (Distance(a, s) <= eps + Constants::Epsilon) { out->Add(a); break; }
                    }
                }
                return UniqueArray(out->ToArray());
            }

            array<Point2D>^ TopologyTools::Closure2D(array<Point2D>^ S, array<Point2D>^ ambient, double eps)
            {
                // closure = S union {points in ambient within eps of S}
                if (S == nullptr) throw gcnew ArgumentNullException("S");
                if (ambient == nullptr) throw gcnew ArgumentNullException("ambient");
                if (eps < 0) throw gcnew InvalidArgumentException("eps must be non-negative.");

                List<Point2D>^ out = gcnew List<Point2D>(ToList<Point2D>(S));
                for each (Point2D a in ambient) {
                    for each (Point2D s in S) {
                        if (Distance(a, s) <= eps + Constants::Epsilon) { out->Add(a); break; }
                    }
                }
                return UniqueArray(out->ToArray());
            }

            array<Point2D>^ TopologyTools::Interior2D(array<Point2D>^ S, array<Point2D>^ ambient, double eps)
            {
                // interior: points p in S such that B(p,eps) subset of S (relative to ambient)
                if (S == nullptr) throw gcnew ArgumentNullException("S");
                if (ambient == nullptr) throw gcnew ArgumentNullException("ambient");
                if (eps < 0) throw gcnew InvalidArgumentException("eps must be non-negative.");

                // build a HashSet-like by comparing Equals (use List scan)
                List<Point2D>^ setS = ToList<Point2D>(S);
                List<Point2D>^ out = gcnew List<Point2D>();

                for each (Point2D p in S) {
                    bool inside = true;
                    for each (Point2D a in ambient) {
                        if (Distance(p, a) <= eps + Constants::Epsilon) {
                            // a must be in S
                            bool inS = false;
                            for each (Point2D q in setS) { if (Utils::AreEqual(a.X, q.X) && Utils::AreEqual(a.Y, q.Y)) { inS = true; break; } }
                            if (!inS) { inside = false; break; }
                        }
                    }
                    if (inside) out->Add(p);
                }
                return UniqueArray(out->ToArray());
            }

            array<Point2D>^ TopologyTools::Boundary2D(array<Point2D>^ S, array<Point2D>^ ambient, double eps)
            {
                array<Point2D>^ cl = Closure2D(S, ambient, eps);
                array<Point2D>^ inr = Interior2D(S, ambient, eps);

                List<Point2D>^ out = gcnew List<Point2D>();
                for each (Point2D p in cl) {
                    bool inInterior = false;
                    for each (Point2D q in inr) { if (Utils::AreEqual(p.X, q.X) && Utils::AreEqual(p.Y, q.Y)) { inInterior = true; break; } }
                    if (!inInterior) out->Add(p);
                }
                return UniqueArray(out->ToArray());
            }

            bool TopologyTools::IsOpen2D(array<Point2D>^ S, array<Point2D>^ ambient, double eps)
            {
                array<Point2D>^ inr = Interior2D(S, ambient, eps);
                // compare sets (as arrays) by membership
                List<Point2D>^ inrL = ToList<Point2D>(inr);
                List<Point2D>^ sL = ToList<Point2D>(S);
                if (inrL->Count != sL->Count) return false;
                for each (Point2D p in sL) {
                    bool found = false;
                    for each (Point2D q in inrL) if (Utils::AreEqual(p.X, q.X) && Utils::AreEqual(p.Y, q.Y)) { found = true; break; }
                    if (!found) return false;
                }
                return true;
            }

            bool TopologyTools::IsClosed2D(array<Point2D>^ S, array<Point2D>^ ambient, double eps)
            {
                array<Point2D>^ cl = Closure2D(S, ambient, eps);
                List<Point2D>^ clL = ToList<Point2D>(cl);
                List<Point2D>^ sL = ToList<Point2D>(S);
                if (clL->Count != sL->Count) return false;
                for each (Point2D p in sL) {
                    bool found = false;
                    for each (Point2D q in clL) if (Utils::AreEqual(p.X, q.X) && Utils::AreEqual(p.Y, q.Y)) { found = true; break; }
                    if (!found) return false;
                }
                return true;
            }

            bool TopologyTools::IsLimitPoint2D(Point2D p, array<Point2D>^ S, double eps)
            {
                if (S == nullptr) throw gcnew ArgumentNullException("S");
                if (eps < 0) throw gcnew InvalidArgumentException("eps must be non-negative.");
                for each (Point2D q in S) {
                    if (Utils::AreEqual(p.X, q.X) && Utils::AreEqual(p.Y, q.Y)) continue;
                    if (Distance(p, q) <= eps + Constants::Epsilon) return true;
                }
                return false;
            }

            bool TopologyTools::IsIsolatedPoint2D(Point2D p, array<Point2D>^ S, double eps)
            {
                if (S == nullptr) throw gcnew ArgumentNullException("S");
                if (eps < 0) throw gcnew InvalidArgumentException("eps must be non-negative.");
                for each (Point2D q in S) {
                    if (Utils::AreEqual(p.X, q.X) && Utils::AreEqual(p.Y, q.Y)) continue;
                    if (Distance(p, q) <= eps + Constants::Epsilon) return false;
                }
                return true;
            }

            // -----------------------
            // Connected components (2D) via BFS (eps graph)
            // -----------------------
            static double dist2D_local(Point2D a, Point2D b) { return Distance(a, b); }

            array<array<Point2D>^>^ TopologyTools::ConnectedComponents2D(array<Point2D>^ S, double eps)
            {
                if (S == nullptr) throw gcnew ArgumentNullException("S");
                if (eps < 0) throw gcnew InvalidArgumentException("eps must be non-negative.");

                int n = S->Length;
                List<array<Point2D>^>^ result = gcnew List<array<Point2D>^>();
                if (n == 0) return result->ToArray();

                array<bool>^ seen = gcnew array<bool>(n);
                for (int i = 0; i < n; ++i) seen[i] = false;

                for (int i = 0; i < n; ++i) {
                    if (seen[i]) continue;
                    // BFS
                    List<int>^ queue = gcnew List<int>();
                    queue->Add(i); seen[i] = true;
                    List<Point2D>^ comp = gcnew List<Point2D>();
                    for (int qi = 0; qi < queue->Count; qi++) {
                        int idx = queue[qi];
                        comp->Add(S[idx]);
                        for (int j = 0; j < n; ++j) {
                            if (seen[j]) continue;
                            if (Distance(S[idx], S[j]) <= eps + Constants::Epsilon) {
                                seen[j] = true;
                                queue->Add(j);
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

            // -----------------------
            // Hausdorff distance 2D
            // -----------------------
            static double onewayHausdorff2D(array<Point2D>^ A, array<Point2D>^ B)
            {
                double best = 0.0;
                for each (Point2D a in A) {
                    double minD = System::Double::PositiveInfinity;
                    for each (Point2D b in B) {
                        double d = Distance(a, b);
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

            // -----------------------
            // Rips (0D) persistence - simple union-find on sorted edges
            // -----------------------
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

                // build list of edges (i,j,dist) with i<j
                List<Tuple<int, int, double>^>^ edges = gcnew List<Tuple<int, int, double>^>();
                for (int i = 0; i < n; ++i) {
                    for (int j = i + 1; j < n; ++j) {
                        double d = Distance(points[i], points[j]);
                        edges->Add(gcnew Tuple<int, int, double>(i, j, d));
                    }
                }
                // sort edges by distance ascending
                edges->Sort(gcnew Comparison<Tuple<int, int, double>^>([](Tuple<int, int, double>^ a, Tuple<int, int, double>^ b) {
                    if (a->Item3 < b->Item3) return -1;
                    if (a->Item3 > b->Item3) return 1;
                    return 0;
                    }));

                // initial: each point alive with birth 0
                array<double>^ birth = gcnew array<double>(n);
                array<double>^ death = gcnew array<double>(n);
                for (int i = 0; i < n; ++i) { birth[i] = 0.0; death[i] = System::Double::PositiveInfinity; }

                UF uf(n);
                // process edges
                for each (Tuple<int, int, double> ^ e in edges) {
                    int i = e->Item1, j = e->Item2;
                    double d = e->Item3;
                    int ri = uf.Find(i), rj = uf.Find(j);
                    if (ri == rj) continue;
                    // when merging components, arbitrarily choose one to die at distance d
                    // choose to kill component with larger representative id to keep deterministic
                    int toDie = (ri > rj) ? ri : rj;
                    death[toDie] = d;
                    uf.Union(ri, rj);
                }

                for (int i = 0; i < n; ++i) {
                    intervals->Add(KeyValuePair<double, double>(birth[i], death[i]));
                }
                return intervals->ToArray();
            }

            // -----------------------
            // 3D: reuse logic but with Point3D and Distance3D
            // -----------------------
            array<Point3D>^ TopologyTools::EpsilonNeighborhood3D(array<Point3D>^ S, array<Point3D>^ ambient, double eps)
            {
                if (S == nullptr) throw gcnew ArgumentNullException("S");
                if (ambient == nullptr) throw gcnew ArgumentNullException("ambient");
                if (eps < 0) throw gcnew InvalidArgumentException("eps must be non-negative.");

                List<Point3D>^ out = gcnew List<Point3D>();
                for each (Point3D a in ambient) {
                    for each (Point3D s in S) {
                        if (Distance3D(a, s) <= eps + Constants::Epsilon) { out->Add(a); break; }
                    }
                }
                return UniqueArray(out->ToArray());
            }

            array<Point3D>^ TopologyTools::Closure3D(array<Point3D>^ S, array<Point3D>^ ambient, double eps)
            {
                if (S == nullptr) throw gcnew ArgumentNullException("S");
                if (ambient == nullptr) throw gcnew ArgumentNullException("ambient");
                if (eps < 0) throw gcnew InvalidArgumentException("eps must be non-negative.");

                List<Point3D>^ out = gcnew List<Point3D>(ToList<Point3D>(S));
                for each (Point3D a in ambient) {
                    for each (Point3D s in S) {
                        if (Distance3D(a, s) <= eps + Constants::Epsilon) { out->Add(a); break; }
                    }
                }
                return UniqueArray(out->ToArray());
            }

            array<Point3D>^ TopologyTools::Interior3D(array<Point3D>^ S, array<Point3D>^ ambient, double eps)
            {
                if (S == nullptr) throw gcnew ArgumentNullException("S");
                if (ambient == nullptr) throw gcnew ArgumentNullException("ambient");
                if (eps < 0) throw gcnew InvalidArgumentException("eps must be non-negative.");

                List<Point3D>^ setS = ToList<Point3D>(S);
                List<Point3D>^ out = gcnew List<Point3D>();

                for each (Point3D p in S) {
                    bool inside = true;
                    for each (Point3D a in ambient) {
                        if (Distance3D(p, a) <= eps + Constants::Epsilon) {
                            bool inS = false;
                            for each (Point3D q in setS) { if (Utils::AreEqual(a.X, q.X) && Utils::AreEqual(a.Y, q.Y) && Utils::AreEqual(a.Z, q.Z)) { inS = true; break; } }
                            if (!inS) { inside = false; break; }
                        }
                    }
                    if (inside) out->Add(p);
                }
                return UniqueArray(out->ToArray());
            }

            array<Point3D>^ TopologyTools::Boundary3D(array<Point3D>^ S, array<Point3D>^ ambient, double eps)
            {
                array<Point3D>^ cl = Closure3D(S, ambient, eps);
                array<Point3D>^ inr = Interior3D(S, ambient, eps);

                List<Point3D>^ out = gcnew List<Point3D>();
                for each (Point3D p in cl) {
                    bool inInterior = false;
                    for each (Point3D q in inr) { if (Utils::AreEqual(p.X, q.X) && Utils::AreEqual(p.Y, q.Y) && Utils::AreEqual(p.Z, q.Z)) { inInterior = true; break; } }
                    if (!inInterior) out->Add(p);
                }
                return UniqueArray(out->ToArray());
            }

            array<array<Point3D>^>^ TopologyTools::ConnectedComponents3D(array<Point3D>^ S, double eps)
            {
                if (S == nullptr) throw gcnew ArgumentNullException("S");
                if (eps < 0) throw gcnew InvalidArgumentException("eps must be non-negative.");

                int n = S->Length;
                List<array<Point3D>^>^ result = gcnew List<array<Point3D>^>();
                if (n == 0) return result->ToArray();

                array<bool>^ seen = gcnew array<bool>(n);
                for (int i = 0; i < n; ++i) seen[i] = false;

                for (int i = 0; i < n; ++i) {
                    if (seen[i]) continue;
                    List<int>^ queue = gcnew List<int>();
                    queue->Add(i); seen[i] = true;
                    List<Point3D>^ comp = gcnew List<Point3D>();
                    for (int qi = 0; qi < queue->Count; qi++) {
                        int idx = queue[qi];
                        comp->Add(S[idx]);
                        for (int j = 0; j < n; ++j) {
                            if (seen[j]) continue;
                            if (Distance3D(S[idx], S[j]) <= eps + Constants::Epsilon) {
                                seen[j] = true;
                                queue->Add(j);
                            }
                        }
                    }
                    result->Add(comp->ToArray());
                }
                return result->ToArray();
            }

            bool TopologyTools::IsConnected3D(array<Point3D>^ S, double eps)
            {
                array<array<Point3D>^>^ comps = ConnectedComponents3D(S, eps);
                return (comps->Length <= 1);
            }

            static double onewayHausdorff3D(array<Point3D>^ A, array<Point3D>^ B)
            {
                double best = 0.0;
                for each (Point3D a in A) {
                    double minD = System::Double::PositiveInfinity;
                    for each (Point3D b in B) {
                        double d = Distance3D(a, b);
                        if (d < minD) minD = d;
                    }
                    if (minD > best) best = minD;
                }
                return best;
            }

            double TopologyTools::HausdorffDistance3D(array<Point3D>^ A, array<Point3D>^ B)
            {
                if (A == nullptr) throw gcnew ArgumentNullException("A");
                if (B == nullptr) throw gcnew ArgumentNullException("B");
                if (A->Length == 0 && B->Length == 0) return 0.0;
                if (A->Length == 0 || B->Length == 0) return System::Double::PositiveInfinity;

                double ab = onewayHausdorff3D(A, B);
                double ba = onewayHausdorff3D(B, A);
                return System::Math::Max(ab, ba);
            }

        } // namespace Topology
    } // namespace Math
} // namespace Agrine

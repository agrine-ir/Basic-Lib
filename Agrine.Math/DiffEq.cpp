#include "DiffEq.h"

using namespace System;
using namespace Agrine::Math::Core;
using namespace Agrine::Math::Algebra;

namespace Agrine {
    namespace Math {
        namespace Calculus {

            // ------------------------------
            // Helpers
            // ------------------------------
            double DiffEq::SafeStepCount(double a, double b, double h)
            {
                if (Utils::IsZero(h)) throw gcnew InvalidArgumentException("Step size h must be non-zero.");
                double len = b - a;
                if (Utils::IsZero(len)) return 0.0;
                return System::Math::Max(1.0, System::Math::Floor(len / h));
            }

            // ==============================
            // Scalar solvers
            // ==============================
            ODESolutionScalar^ DiffEq::EulerScalar(Func<double, double, double>^ f,
                double t0, double y0, double tf, double h)
            {
                if (f == nullptr) throw gcnew ArgumentNullException("f");
                if (h <= 0) throw gcnew InvalidArgumentException("h must be positive.");

                // compute step count (ensure include final point)
                int n = (int)(System::Math::Ceiling((tf - t0) / h));
                if (n < 0) throw gcnew InvalidArgumentException("tf must be >= t0.");

                array<double>^ T = gcnew array<double>(n + 1);
                array<double>^ Y = gcnew array<double>(n + 1);

                double t = t0;
                double y = y0;

                T[0] = t; Y[0] = y;

                for (int i = 1; i <= n; i++)
                {
                    double hstep = (i == n) ? (tf - t) : h; // last step adjust
                    double dy = f(t, y);
                    y = y + hstep * dy;
                    t = t + hstep;
                    T[i] = t; Y[i] = y;
                }

                return gcnew ODESolutionScalar(T, Y);
            }

            ODESolutionScalar^ DiffEq::HeunScalar(Func<double, double, double>^ f,
                double t0, double y0, double tf, double h)
            {
                if (f == nullptr) throw gcnew ArgumentNullException("f");
                if (h <= 0) throw gcnew InvalidArgumentException("h must be positive.");

                int n = (int)(System::Math::Ceiling((tf - t0) / h));
                if (n < 0) throw gcnew InvalidArgumentException("tf must be >= t0.");

                array<double>^ T = gcnew array<double>(n + 1);
                array<double>^ Y = gcnew array<double>(n + 1);

                double t = t0;
                double y = y0;

                T[0] = t; Y[0] = y;

                for (int i = 1; i <= n; i++)
                {
                    double hstep = (i == n) ? (tf - t) : h;
                    double k1 = f(t, y);
                    double y_pred = y + hstep * k1;
                    double k2 = f(t + hstep, y_pred);
                    y = y + (hstep * 0.5) * (k1 + k2);
                    t = t + hstep;
                    T[i] = t; Y[i] = y;
                }

                return gcnew ODESolutionScalar(T, Y);
            }

            ODESolutionScalar^ DiffEq::RK4Scalar(Func<double, double, double>^ f,
                double t0, double y0, double tf, double h)
            {
                if (f == nullptr) throw gcnew ArgumentNullException("f");
                if (h <= 0) throw gcnew InvalidArgumentException("h must be positive.");

                int n = (int)(System::Math::Ceiling((tf - t0) / h));
                if (n < 0) throw gcnew InvalidArgumentException("tf must be >= t0.");

                array<double>^ T = gcnew array<double>(n + 1);
                array<double>^ Y = gcnew array<double>(n + 1);

                double t = t0;
                double y = y0;
                T[0] = t; Y[0] = y;

                for (int i = 1; i <= n; i++)
                {
                    double hstep = (i == n) ? (tf - t) : h;

                    double k1 = f(t, y);
                    double k2 = f(t + 0.5 * hstep, y + 0.5 * hstep * k1);
                    double k3 = f(t + 0.5 * hstep, y + 0.5 * hstep * k2);
                    double k4 = f(t + hstep, y + hstep * k3);

                    y = y + (hstep / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
                    t = t + hstep;
                    T[i] = t; Y[i] = y;
                }

                return gcnew ODESolutionScalar(T, Y);
            }

            // Adaptive RKF45 coefficients (classic Fehlberg) simplified implementation
            // We implement a basic adaptive stepper: try step h, compute 4th and 5th order estimates,
            // estimate error and adapt h. This is an in-code implementation (not full Dormand-Prince table)
            // but sufficient as a robust adaptive RK stepper for many problems.

            ODESolutionScalar^ DiffEq::AdaptiveRK45Scalar(Func<double, double, double>^ f,
                double t0, double y0, double tf,
                double h0, double tol)
            {
                if (f == nullptr) throw gcnew ArgumentNullException("f");
                if (h0 <= 0) throw gcnew InvalidArgumentException("Initial step h0 must be positive.");
                if (tol <= 0) throw gcnew InvalidArgumentException("tol must be positive.");

                List<double>^ times = gcnew List<double>();
                List<double>^ values = gcnew List<double>();

                double t = t0;
                double y = y0;
                double h = h0;

                times->Add(t);
                values->Add(y);

                while (t < tf)
                {
                    if (t + h > tf) h = tf - t;

                    // classical RK4 estimate
                    double k1 = f(t, y);
                    double k2 = f(t + 0.5 * h, y + 0.5 * h * k1);
                    double k3 = f(t + 0.5 * h, y + 0.5 * h * k2);
                    double k4 = f(t + h, y + h * k3);
                    double y_rk4 = y + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);

                    // simple embedded 5th order estimate using smaller step (two half RK4 steps)
                    // do two half-steps of RK4 to get higher-resolution y2
                    double hh = h * 0.5;
                    // first half
                    double k1a = f(t, y);
                    double k2a = f(t + 0.5 * hh, y + 0.5 * hh * k1a);
                    double k3a = f(t + 0.5 * hh, y + 0.5 * hh * k2a);
                    double k4a = f(t + hh, y + hh * k3a);
                    double y_mid = y + (hh / 6.0) * (k1a + 2 * k2a + 2 * k3a + k4a);

                    // second half
                    double k1b = f(t + hh, y_mid);
                    double k2b = f(t + hh + 0.5 * hh, y_mid + 0.5 * hh * k1b);
                    double k3b = f(t + hh + 0.5 * hh, y_mid + 0.5 * hh * k2b);
                    double k4b = f(t + h, y_mid + hh * k3b);
                    double y_two_half = y_mid + (hh / 6.0) * (k1b + 2 * k2b + 2 * k3b + k4b);

                    // estimate error between single-step RK4 and two-half-step RK4
                    double err = System::Math::Abs(y_two_half - y_rk4);

                    // safety factors
                    double SAFETY = 0.9;
                    double MIN_SCALE = 0.2;
                    double MAX_SCALE = 5.0;

                    // compute new step
                    if (err <= tol)
                    {
                        // accept step: use higher-accuracy value y_two_half
                        t += h;
                        y = y_two_half;
                        times->Add(t);
                        values->Add(y);

                        // increase h
                        double scale = SAFETY * System::Math::Pow(tol / (err + 1e-16), 0.25);
                        scale = System::Math::Max(MIN_SCALE, System::Math::Min(MAX_SCALE, scale));
                        h = h * scale;
                    }
                    else
                    {
                        // reject step, reduce h
                        double scale = SAFETY * System::Math::Pow(tol / (err + 1e-16), 0.25);
                        scale = System::Math::Max(MIN_SCALE, System::Math::Min(MAX_SCALE, scale));
                        h = h * scale;
                        // try again with smaller h
                    }

                    // guard: avoid extremely small step
                    if (h < Constants::Epsilon * 10.0)
                        throw gcnew ConvergenceException("Step size underflow in AdaptiveRK45Scalar.");
                }

                // convert Lists to arrays
                array<double>^ T = times->ToArray();
                array<double>^ Y = values->ToArray();
                return gcnew ODESolutionScalar(T, Y);
            }

            // ==============================
            // System solvers (Vector)
            // ==============================
            ODESolutionSystem^ DiffEq::EulerSystem(Func<double, Vector^, Vector^>^ f,
                double t0, Vector^ y0, double tf, double h)
            {
                if (f == nullptr) throw gcnew ArgumentNullException("f");
                if (y0 == nullptr) throw gcnew ArgumentNullException("y0");
                if (h <= 0) throw gcnew InvalidArgumentException("h must be positive.");

                int n = (int)(System::Math::Ceiling((tf - t0) / h));
                if (n < 0) throw gcnew InvalidArgumentException("tf must be >= t0.");

                array<double>^ T = gcnew array<double>(n + 1);
                array<Vector^>^ S = gcnew array<Vector^>(n + 1);

                double t = t0;
                Vector^ y = y0->Clone();

                T[0] = t; S[0] = y->Clone();

                for (int i = 1; i <= n; i++)
                {
                    double hstep = (i == n) ? (tf - t) : h;
                    Vector^ dy = f(t, y);
                    Vector^ ynew = y->Add(dy->Multiply(hstep));
                    t += hstep;
                    y = ynew;
                    T[i] = t; S[i] = y->Clone();
                }

                return gcnew ODESolutionSystem(T, S);
            }

            ODESolutionSystem^ DiffEq::HeunSystem(Func<double, Vector^, Vector^>^ f,
                double t0, Vector^ y0, double tf, double h)
            {
                if (f == nullptr) throw gcnew ArgumentNullException("f");
                if (y0 == nullptr) throw gcnew ArgumentNullException("y0");
                if (h <= 0) throw gcnew InvalidArgumentException("h must be positive.");

                int n = (int)(System::Math::Ceiling((tf - t0) / h));
                if (n < 0) throw gcnew InvalidArgumentException("tf must be >= t0.");

                array<double>^ T = gcnew array<double>(n + 1);
                array<Vector^>^ S = gcnew array<Vector^>(n + 1);

                double t = t0;
                Vector^ y = y0->Clone();

                T[0] = t; S[0] = y->Clone();

                for (int i = 1; i <= n; i++)
                {
                    double hstep = (i == n) ? (tf - t) : h;
                    Vector^ k1 = f(t, y);
                    Vector^ y_pred = y->Add(k1->Multiply(hstep));
                    Vector^ k2 = f(t + hstep, y_pred);
                    Vector^ ynew = y->Add(k1->Multiply(0.5 * hstep))->Add(k2->Multiply(0.5 * hstep)); // average
                    t += hstep;
                    y = ynew;
                    T[i] = t; S[i] = y->Clone();
                }

                return gcnew ODESolutionSystem(T, S);
            }

            ODESolutionSystem^ DiffEq::RK4System(Func<double, Vector^, Vector^>^ f,
                double t0, Vector^ y0, double tf, double h)
            {
                if (f == nullptr) throw gcnew ArgumentNullException("f");
                if (y0 == nullptr) throw gcnew ArgumentNullException("y0");
                if (h <= 0) throw gcnew InvalidArgumentException("h must be positive.");

                int n = (int)(System::Math::Ceiling((tf - t0) / h));
                if (n < 0) throw gcnew InvalidArgumentException("tf must be >= t0.");

                array<double>^ T = gcnew array<double>(n + 1);
                array<Vector^>^ S = gcnew array<Vector^>(n + 1);

                double t = t0;
                Vector^ y = y0->Clone();

                T[0] = t; S[0] = y->Clone();

                for (int i = 1; i <= n; i++)
                {
                    double hstep = (i == n) ? (tf - t) : h;

                    Vector^ k1 = f(t, y);
                    Vector^ k2 = f(t + 0.5 * hstep, y->Add(k1->Multiply(0.5 * hstep)));
                    Vector^ k3 = f(t + 0.5 * hstep, y->Add(k2->Multiply(0.5 * hstep)));
                    Vector^ k4 = f(t + hstep, y->Add(k3->Multiply(hstep)));

                    // y_next = y + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
                    Vector^ sum = k1->Multiply(1.0)
                        ->Add(k2->Multiply(2.0))
                        ->Add(k3->Multiply(2.0))
                        ->Add(k4->Multiply(1.0));

                    Vector^ ynew = y->Add(sum->Multiply(hstep / 6.0));

                    t += hstep;
                    y = ynew;
                    T[i] = t; S[i] = y->Clone();
                }

                return gcnew ODESolutionSystem(T, S);
            }

            ODESolutionSystem^ DiffEq::AdaptiveRK45System(Func<double, Vector^, Vector^>^ f,
                double t0, Vector^ y0, double tf,
                double h0, double tol)
            {
                if (f == nullptr) throw gcnew ArgumentNullException("f");
                if (y0 == nullptr) throw gcnew ArgumentNullException("y0");
                if (h0 <= 0) throw gcnew InvalidArgumentException("Initial step h0 must be positive.");
                if (tol <= 0) throw gcnew InvalidArgumentException("tol must be positive.");

                List<double>^ times = gcnew List<double>();
                List<Vector^>^ states = gcnew List<Vector^>();

                double t = t0;
                Vector^ y = y0->Clone();
                double h = h0;

                times->Add(t);
                states->Add(y->Clone());

                while (t < tf)
                {
                    if (t + h > tf) h = tf - t;

                    // single rk4 step
                    Vector^ k1 = f(t, y);
                    Vector^ k2 = f(t + 0.5 * h, y->Add(k1->Multiply(0.5 * h)));
                    Vector^ k3 = f(t + 0.5 * h, y->Add(k2->Multiply(0.5 * h)));
                    Vector^ k4 = f(t + h, y->Add(k3->Multiply(h)));
                    Vector^ y_rk4 = y->Add(k1->Multiply(h / 6.0)->Add(k2->Multiply(h / 3.0))->Add(k3->Multiply(h / 3.0))->Add(k4->Multiply(h / 6.0)));

                    // two half steps
                    double hh = 0.5 * h;
                    // first half
                    Vector^ k1a = f(t, y);
                    Vector^ k2a = f(t + 0.5 * hh, y->Add(k1a->Multiply(0.5 * hh)));
                    Vector^ k3a = f(t + 0.5 * hh, y->Add(k2a->Multiply(0.5 * hh)));
                    Vector^ k4a = f(t + hh, y->Add(k3a->Multiply(hh)));
                    Vector^ y_mid = y->Add(k1a->Multiply(hh / 6.0)->Add(k2a->Multiply(hh / 3.0))->Add(k3a->Multiply(hh / 3.0))->Add(k4a->Multiply(hh / 6.0)));

                    // second half
                    Vector^ k1b = f(t + hh, y_mid);
                    Vector^ k2b = f(t + hh + 0.5 * hh, y_mid->Add(k1b->Multiply(0.5 * hh)));
                    Vector^ k3b = f(t + hh + 0.5 * hh, y_mid->Add(k2b->Multiply(0.5 * hh)));
                    Vector^ k4b = f(t + h, y_mid->Add(k3b->Multiply(hh)));
                    Vector^ y_two_half = y_mid->Add(k1b->Multiply(hh / 6.0)->Add(k2b->Multiply(hh / 3.0))->Add(k3b->Multiply(hh / 3.0))->Add(k4b->Multiply(hh / 6.0)));

                    // compute error as vector norm of difference
                    Vector^ errVec = y_two_half->Subtract(y_rk4);
                    double err = errVec->Norm();

                    // safety and scaling
                    double SAFETY = 0.9;
                    double MIN_SCALE = 0.2;
                    double MAX_SCALE = 5.0;

                    if (err <= tol)
                    {
                        // accept higher resolution solution
                        t += h;
                        y = y_two_half->Clone();
                        times->Add(t);
                        states->Add(y->Clone());

                        double scale = SAFETY * System::Math::Pow(tol / (err + 1e-16), 0.25);
                        scale = System::Math::Max(MIN_SCALE, System::Math::Min(MAX_SCALE, scale));
                        h = h * scale;
                    }
                    else
                    {
                        double scale = SAFETY * System::Math::Pow(tol / (err + 1e-16), 0.25);
                        scale = System::Math::Max(MIN_SCALE, System::Math::Min(MAX_SCALE, scale));
                        h = h * scale;
                    }

                    if (h < Constants::Epsilon * 10.0)
                        throw gcnew ConvergenceException("Step size underflow in AdaptiveRK45System.");
                }

                array<double>^ T = times->ToArray();
                array<Vector^>^ S = states->ToArray();
                return gcnew ODESolutionSystem(T, S);
            }

        } // namespace Calculus
    } // namespace Math
} // namespace Agrine

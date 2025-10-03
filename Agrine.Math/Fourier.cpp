#include "Fourier.h"
#include <cmath>

using namespace System;
using namespace System::Collections::Generic;
using namespace Agrine::Math::Advanced;
using namespace Agrine::Math::Core;

namespace Agrine {
    namespace Math {
        namespace Advanced {

            // -------------------------
            // Helpers
            // -------------------------
            int Fourier::BitReverse(int x, int bits)
            {
                int y = 0;
                for (int i = 0; i < bits; ++i) {
                    y = (y << 1) | (x & 1);
                    x >>= 1;
                }
                return y;
            }

            int Fourier::NextPowerOfTwo(int n)
            {
                if (n <= 1) return 1;
                int p = 1;
                while (p < n) p <<= 1;
                return p;
            }

            array<Complex>^ Fourier::ZeroPad(array<Complex>^ input, int newSize)
            {
                if (input == nullptr) throw gcnew ArgumentNullException("input");
                if (newSize < input->Length) throw gcnew InvalidArgumentException("newSize must be >= input length.");
                array<Complex>^ out = gcnew array<Complex>(newSize);
                for (int i = 0; i < input->Length; ++i) out[i] = input[i];
                for (int i = input->Length; i < newSize; ++i) out[i] = Complex::Zero();
                return out;
            }

            // -------------------------
            // Naive DFT / IDFT
            // -------------------------
            array<Complex>^ Fourier::DFT(array<Complex>^ input)
            {
                if (input == nullptr) throw gcnew ArgumentNullException("input");
                int N = input->Length;
                array<Complex>^ out = gcnew array<Complex>(N);
                double twoPi = 2.0 * Constants::Pi;
                for (int k = 0; k < N; ++k) {
                    Complex sum = Complex::Zero();
                    for (int n = 0; n < N; ++n) {
                        double angle = -twoPi * k * n / (double)N;
                        Complex tw = Complex::FromPolar(1.0, angle);
                        sum = sum + (input[n] * tw);
                    }
                    out[k] = sum;
                }
                return out;
            }

            array<Complex>^ Fourier::IDFT(array<Complex>^ input)
            {
                if (input == nullptr) throw gcnew ArgumentNullException("input");
                int N = input->Length;
                array<Complex>^ out = gcnew array<Complex>(N);
                double twoPi = 2.0 * Constants::Pi;
                for (int n = 0; n < N; ++n) {
                    Complex sum = Complex::Zero();
                    for (int k = 0; k < N; ++k) {
                        double angle = twoPi * k * n / (double)N; // positive sign
                        Complex tw = Complex::FromPolar(1.0, angle);
                        sum = sum + (input[k] * tw);
                    }
                    // scale by 1/N
                    out[n] = sum / (double)N;
                }
                return out;
            }

            // -------------------------
            // FFT (radix-2 iterative)
            // -------------------------
            array<Complex>^ Fourier::FFT(array<Complex>^ input)
            {
                if (input == nullptr) throw gcnew ArgumentNullException("input");
                int N = input->Length;
                // if N is not power of two, fall back to DFT (safer)
                if ((N & (N - 1)) != 0) {
                    return DFT(input);
                }

                // copy input
                array<Complex>^ a = gcnew array<Complex>(N);
                for (int i = 0; i < N; ++i) a[i] = input[i];

                // bit-reversal permutation
                int bits = 0;
                while ((1 << bits) < N) ++bits;
                for (int i = 0; i < N; ++i) {
                    int j = BitReverse(i, bits);
                    if (j > i) {
                        Complex tmp = a[i]; a[i] = a[j]; a[j] = tmp;
                    }
                }

                // Cooley-Tukey iterative
                double twoPi = 2.0 * Constants::Pi;
                for (int len = 2; len <= N; len <<= 1) {
                    int half = len >> 1;
                    // wlen = exp(-2pi i / len)
                    for (int i = 0; i < N; i += len) {
                        for (int j = 0; j < half; ++j) {
                            double angle = -twoPi * j / (double)len;
                            Complex w = Complex::FromPolar(1.0, angle);
                            Complex u = a[i + j];
                            Complex v = a[i + j + half] * w;
                            // explicit component-wise add/sub to avoid operator overload ambiguity
                            a[i + j] = Complex(u.Re + v.Re, u.Im + v.Im);
                            a[i + j + half] = Complex(u.Re - v.Re, u.Im - v.Im);

                        }
                    }
                }

                return a;
            }

            array<Complex>^ Fourier::IFFT(array<Complex>^ input)
            {
                if (input == nullptr) throw gcnew ArgumentNullException("input");
                // IFFT via conjugate trick: ifft(x) = conj( fft(conj(x)) ) / N
                int N = input->Length;
                array<Complex>^ conjIn = gcnew array<Complex>(N);
                for (int i = 0; i < N; ++i) conjIn[i] = input[i].Conjugate();
                array<Complex>^ y = FFT(conjIn);
                array<Complex>^ out = gcnew array<Complex>(N);
                for (int i = 0; i < N; ++i) out[i] = y[i].Conjugate() / (double)N;
                return out;
            }

            // -------------------------
            // Real FFT helpers
            // -------------------------
            array<Complex>^ Fourier::FFTReal(array<double>^ realInput)
            {
                if (realInput == nullptr) throw gcnew ArgumentNullException("realInput");
                int n = realInput->Length;
                int N = NextPowerOfTwo(n);
                // copy into complex array and zero-pad
                array<Complex>^ a = gcnew array<Complex>(N);
                for (int i = 0; i < n; ++i) a[i] = Complex(realInput[i], 0.0);
                for (int i = n; i < N; ++i) a[i] = Complex::Zero();
                return FFT(a);
            }

            array<double>^ Fourier::IFFTReal(array<Complex>^ spectrum, int originalLength)
            {
                if (spectrum == nullptr) throw gcnew ArgumentNullException("spectrum");
                if (originalLength < 0) throw gcnew InvalidArgumentException("originalLength must be non-negative.");
                array<Complex>^ t = IFFT(spectrum);
                int N = t->Length;
                int outLen = System::Math::Min(originalLength, N);
                array<double>^ out = gcnew array<double>(outLen);
                for (int i = 0; i < outLen; ++i) out[i] = t[i].Real;
                return out;
            }

            // -------------------------
            // Convolution
            // -------------------------
            array<Complex>^ Fourier::Convolution(array<Complex>^ a, array<Complex>^ b)
            {
                if (a == nullptr) throw gcnew ArgumentNullException("a");
                if (b == nullptr) throw gcnew ArgumentNullException("b");
                int na = a->Length, nb = b->Length;
                if (na == 0 || nb == 0) return gcnew array<Complex>(0);
                int needed = na + nb - 1;
                int N = NextPowerOfTwo(needed);
                array<Complex>^ A = ZeroPad(a, N);
                array<Complex>^ B = ZeroPad(b, N);
                array<Complex>^ FA = FFT(A);
                array<Complex>^ FB = FFT(B);
                array<Complex>^ FC = gcnew array<Complex>(N);
                for (int i = 0; i < N; ++i) FC[i] = FA[i] * FB[i];
                array<Complex>^ c = IFFT(FC);
                // truncate to needed length
                array<Complex>^ out = gcnew array<Complex>(needed);
                for (int i = 0; i < needed; ++i) out[i] = c[i];
                return out;
            }

            array<double>^ Fourier::ConvolutionReal(array<double>^ a, array<double>^ b)
            {
                if (a == nullptr) throw gcnew ArgumentNullException("a");
                if (b == nullptr) throw gcnew ArgumentNullException("b");
                int na = a->Length, nb = b->Length;
                if (na == 0 || nb == 0) return gcnew array<double>(0);
                // convert to complex arrays
                array<Complex>^ A = gcnew array<Complex>(na);
                array<Complex>^ B = gcnew array<Complex>(nb);
                for (int i = 0; i < na; ++i) A[i] = Complex(a[i], 0.0);
                for (int i = 0; i < nb; ++i) B[i] = Complex(b[i], 0.0);
                array<Complex>^ C = Convolution(A, B);
                int needed = na + nb - 1;
                array<double>^ out = gcnew array<double>(needed);
                for (int i = 0; i < needed; ++i) out[i] = C[i].Real;
                return out;
            }

            // -------------------------
            // Window functions
            // -------------------------
            array<double>^ Fourier::WindowHann(int N)
            {
                if (N <= 0) throw gcnew InvalidArgumentException("N must be positive.");
                array<double>^ w = gcnew array<double>(N);
                if (N == 1) { w[0] = 1.0; return w; }
                double twoPi = 2.0 * Constants::Pi;
                for (int n = 0; n < N; ++n) {
                    w[n] = 0.5 * (1.0 - std::cos(twoPi * n / (double)(N - 1)));
                }
                return w;
            }

            array<double>^ Fourier::WindowHamming(int N)
            {
                if (N <= 0) throw gcnew InvalidArgumentException("N must be positive.");
                array<double>^ w = gcnew array<double>(N);
                if (N == 1) { w[0] = 1.0; return w; }
                double twoPi = 2.0 * Constants::Pi;
                for (int n = 0; n < N; ++n) {
                    w[n] = 0.54 - 0.46 * std::cos(twoPi * n / (double)(N - 1));
                }
                return w;
            }

            array<double>^ Fourier::WindowBlackman(int N)
            {
                if (N <= 0) throw gcnew InvalidArgumentException("N must be positive.");
                array<double>^ w = gcnew array<double>(N);
                if (N == 1) { w[0] = 1.0; return w; }
                double twoPi = 2.0 * Constants::Pi;
                for (int n = 0; n < N; ++n) {
                    double a0 = 0.42, a1 = 0.5, a2 = 0.08;
                    w[n] = a0 - a1 * std::cos(twoPi * n / (double)(N - 1)) + a2 * std::cos(2.0 * twoPi * n / (double)(N - 1));
                }
                return w;
            }

            // -------------------------
            // Spectrum helpers
            // -------------------------
            array<double>^ Fourier::MagnitudeSpectrum(array<Complex>^ spectrum)
            {
                if (spectrum == nullptr) throw gcnew ArgumentNullException("spectrum");
                int N = spectrum->Length;
                array<double>^ out = gcnew array<double>(N);
                for (int i = 0; i < N; ++i) out[i] = spectrum[i].Abs();
                return out;
            }

            array<double>^ Fourier::PhaseSpectrum(array<Complex>^ spectrum)
            {
                if (spectrum == nullptr) throw gcnew ArgumentNullException("spectrum");
                int N = spectrum->Length;
                array<double>^ out = gcnew array<double>(N);
                for (int i = 0; i < N; ++i) out[i] = spectrum[i].Arg();
                return out;
            }

            // -------------------------
            // Fourier series coefficients from samples (exponential form)
            // -------------------------
            array<Complex>^ Fourier::FourierSeriesCoefficientsFromSamples(array<double>^ samples)
            {
                if (samples == nullptr) throw gcnew ArgumentNullException("samples");
                int N = samples->Length;
                if (N == 0) return gcnew array<Complex>(0);
                // convert samples to complex
                array<Complex>^ s = gcnew array<Complex>(N);
                for (int i = 0; i < N; ++i) s[i] = Complex(samples[i], 0.0);
                // compute DFT
                array<Complex>^ C = DFT(s);
                // scale by 1/N
                for (int k = 0; k < N; ++k) C[k] = C[k] / (double)N;
                return C;
            }

        } // namespace Advanced
    } // namespace Math
} // namespace Agrine

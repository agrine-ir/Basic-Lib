#pragma once

// Fourier.h - DFT/FFT, convolution, windows and spectral helpers
// Part of Agrine.Math (C++/CLI, .NET Framework 4.6.2)

#include "Constants.h"
#include "Utils.h"
#include "Exceptions.h"
#include "Complex.h"

using namespace System;
using namespace System::Collections::Generic;
using namespace Agrine::Math::Advanced;
using namespace Agrine::Math::Core;

namespace Agrine {
    namespace Math {
        namespace Advanced {

            /// <summary>
            /// Fourier - Discrete Fourier Transform utilities (complex and real),
            /// FFT (radix-2 iterative), convolution and window functions.
            /// </summary>
            public ref class Fourier abstract sealed
            {
            public:
                // ----- DFT / IDFT (naive O(N^2)) -----
                static array<Complex>^ DFT(array<Complex>^ input);
                static array<Complex>^ IDFT(array<Complex>^ input);

                // ----- FFT / IFFT (radix-2 iterative). If length isn't power-of-two,
                //        FFT falls back to DFT (safe but slower). -----
                static array<Complex>^ FFT(array<Complex>^ input);
                static array<Complex>^ IFFT(array<Complex>^ input);

                // ----- real signal helpers -----
                // FFTReal: returns complex spectrum (size = nextPow2(N) or N if already power-of-two)
                static array<Complex>^ FFTReal(array<double>^ realInput);
                // IFFTReal: inverse for real-spectrum (returns real time-series truncated to original length)
                static array<double>^ IFFTReal(array<Complex>^ spectrum, int originalLength);

                // ----- convolution (linear) -----
                // returns linear convolution result of length (len(a)+len(b)-1)
                static array<Complex>^ Convolution(array<Complex>^ a, array<Complex>^ b);
                static array<double>^ ConvolutionReal(array<double>^ a, array<double>^ b);

                // ----- windows -----
                static array<double>^ WindowHann(int N);
                static array<double>^ WindowHamming(int N);
                static array<double>^ WindowBlackman(int N);

                // ----- spectrum helpers -----
                static array<double>^ MagnitudeSpectrum(array<Complex>^ spectrum);
                static array<double>^ PhaseSpectrum(array<Complex>^ spectrum);

                // ----- utility -----
                static int NextPowerOfTwo(int n);
                static array<Complex>^ ZeroPad(array<Complex>^ input, int newSize);

                // ----- Fourier series (exponential coefficients) -----
                // Given sampled real function on N equally spaced points over one period,
                // returns c_k (k = 0..N-1) where c_k = (1/N) * sum_n f_n * exp(-i*2pi*k*n/N)
                static array<Complex>^ FourierSeriesCoefficientsFromSamples(array<double>^ samples);

            private:
                // bit-reversal helper for FFT
                static int BitReverse(int x, int bits);
            };

        } // namespace Advanced
    } // namespace Math
} // namespace Agrine

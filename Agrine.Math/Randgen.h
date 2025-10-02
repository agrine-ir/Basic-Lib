#pragma once

// Randgen.h - Advanced random number generators
// Part of Agrine.Math (C++/CLI, .NET Framework 4.6.2)

#include "Constants.h"
#include "Exceptions.h"

using namespace System;

namespace Agrine {
    namespace Math {
        namespace Probability {

            /// <summary>
            /// Randgen - provides advanced random number generators (LCG, MersenneTwister).
            /// </summary>
            public ref class Randgen abstract sealed
            {
            public:
                // -------- Linear Congruential Generator (LCG) --------
                ref class LCG
                {
                private:
                    unsigned long long state;
                    static const unsigned long long a = 6364136223846793005ULL;
                    static const unsigned long long c = 1ULL;
                    static const unsigned long long m = (1ULL << 63);

                public:
                    LCG(unsigned long long seed);
                    void Reseed(unsigned long long seed);
                    double NextDouble();      // uniform [0,1)
                    int NextInt(int max);     // [0, max)
                    int NextRange(int min, int max);
                };

                // -------- Mersenne Twister (MT19937 simplified) --------
                ref class MersenneTwister
                {
                private:
                    literal int N = 624;
                    literal int M = 397;
                    literal unsigned int MATRIX_A = 0x9908b0dfUL;
                    literal unsigned int UPPER_MASK = 0x80000000UL;
                    literal unsigned int LOWER_MASK = 0x7fffffffUL;

                    array<unsigned int>^ mt;
                    int mti;

                public:
                    MersenneTwister(unsigned int seed);
                    void Reseed(unsigned int seed);
                    unsigned int NextUInt32();
                    double NextDouble();      // uniform [0,1)
                    int NextInt(int max);     // [0, max)
                    int NextRange(int min, int max);
                    double NextNormal(double mu, double sigma); // Gaussian via Box-Muller
                };
            };

        } // namespace Probability
    } // namespace Math
} // namespace Agrine

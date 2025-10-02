#include "Randgen.h"
#include <cmath>

using namespace System;
using namespace Agrine::Math::Core;

namespace Agrine {
    namespace Math {
        namespace Probability {

            // ================= LCG =================
            Randgen::LCG::LCG(unsigned long long seed)
            {
                if (seed == 0) seed = 88172645463325252ULL; // default nonzero seed
                state = seed;
            }

            void Randgen::LCG::Reseed(unsigned long long seed)
            {
                if (seed == 0) seed = 88172645463325252ULL;
                state = seed;
            }

            double Randgen::LCG::NextDouble()
            {
                state = (a * state + c) % m;
                return (double)state / (double)m;
            }

            int Randgen::LCG::NextInt(int max)
            {
                if (max <= 0) throw gcnew InvalidArgumentException("LCG.NextInt: max must be positive.");
                return (int)(NextDouble() * max);
            }

            int Randgen::LCG::NextRange(int min, int max)
            {
                if (min >= max) throw gcnew InvalidArgumentException("LCG.NextRange: min must be < max.");
                return min + (int)(NextDouble() * (max - min));
            }

            // ================= Mersenne Twister =================
            Randgen::MersenneTwister::MersenneTwister(unsigned int seed)
            {
                mt = gcnew array<unsigned int>(N);
                Reseed(seed);
            }

            void Randgen::MersenneTwister::Reseed(unsigned int seed)
            {
                mt[0] = seed;
                for (mti = 1; mti < N; ++mti)
                    mt[mti] = (1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
            }

            unsigned int Randgen::MersenneTwister::NextUInt32()
            {
                unsigned int y;
                static unsigned int mag01[2] = { 0x0UL, MATRIX_A };

                if (mti >= N) {
                    int kk;
                    for (kk = 0; kk < N - M; kk++) {
                        y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
                        mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1UL];
                    }
                    for (; kk < N - 1; kk++) {
                        y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
                        mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
                    }
                    y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
                    mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];

                    mti = 0;
                }

                y = mt[mti++];

                // Tempering
                y ^= (y >> 11);
                y ^= (y << 7) & 0x9d2c5680UL;
                y ^= (y << 15) & 0xefc60000UL;
                y ^= (y >> 18);

                return y;
            }

            double Randgen::MersenneTwister::NextDouble()
            {
                return (double)NextUInt32() * (1.0 / 4294967296.0); // divide by 2^32
            }

            int Randgen::MersenneTwister::NextInt(int max)
            {
                if (max <= 0) throw gcnew InvalidArgumentException("MT.NextInt: max must be positive.");
                return (int)(NextDouble() * max);
            }

            int Randgen::MersenneTwister::NextRange(int min, int max)
            {
                if (min >= max) throw gcnew InvalidArgumentException("MT.NextRange: min must be < max.");
                return min + (int)(NextDouble() * (max - min));
            }

            double Randgen::MersenneTwister::NextNormal(double mu, double sigma)
            {
                if (sigma <= 0.0) throw gcnew InvalidArgumentException("MT.NextNormal: sigma must be > 0.");
                // Box-Muller transform
                double u1 = NextDouble();
                double u2 = NextDouble();
                if (u1 < 1e-300) u1 = 1e-300;
                double r = std::sqrt(-2.0 * std::log(u1));
                double theta = 2.0 * Agrine::Math::Core::Constants::Pi * u2;
                double z = r * std::cos(theta);
                return mu + sigma * z;
            }

        } // namespace Probability
    } // namespace Math
} // namespace Agrine

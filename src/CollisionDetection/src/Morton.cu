//
// Created by David Matthews on 5/21/20.
//

#include <bitset>

#include "../include/Morton.cuh"

__host__ __device__
void init_morton_func::operator()(unsigned int idx) {
    mIds[idx] = idx;
    mPoints[idx] = mortonEncode64(xRanks[idx], yRanks[idx], zRanks[idx]);
}

__host__ __device__
void init_morton_func_fast::operator()(unsigned int idx) {
    mIds[idx] = idx;
    mPoints[idx] = mortonEncode64((unsigned long long int) min(max(((x[idx] - xmin) / xrange) * expansionRange, 0.0f),
            expansionRange - 1.f),
            (unsigned long long int) min(max(((y[idx] - ymin) / yrange) * expansionRange, 0.0f), expansionRange - 1.f),
            (unsigned long long int) min(max(((z[idx] - zmin) / zrange) * expansionRange, 0.0f), expansionRange - 1.f));
}


std::ostream &operator<<(std::ostream &out, const Morton &m) {
    return out << std::bitset<64>(m.morton_code);
}

__device__ __host__ void morton_func::operator()(int idx) {
    mortonB[idx].morton_code = mortonA[idx];
}


// https://www.forceflow.be/2013/10/07/morton-encodingdecoding-through-bit-interleaving-implementations/
// method to seperate bits from a given integer 3 positions apart
__host__ __device__
unsigned long long int expandBits64(unsigned long long int a) {
    unsigned long long int x;
    // if a uses more than 21 bits, shift to only use the 21 most significant bits.
#if defined(__CUDA_ARCH__)
    x = a >> (43 - min(43, __clzll(a)));
#else
    if (a == 0) {
        x = a;
    } else {
        x = a>>(43 - min(43, __builtin_clzll(a)));
    }
#endif


    x = (x | x << 32) & 0x1f00000000ffff;
    x = (x | x << 16) & 0x1f0000ff0000ff;
    x = (x | x << 8) & 0x100f00f00f00f00f;
    x = (x | x << 4) & 0x10c30c30c30c30c3;
    x = (x | x << 2) & 0x1249249249249249;
    return x;
}

__host__ __device__
unsigned long long int mortonEncode64(unsigned long long int x, unsigned long long int y, unsigned long long int z) {
    unsigned long long int answer = 0;
    answer |= expandBits64(x) | expandBits64(y) << 1 | expandBits64(z) << 2;
    return answer;
}

__host__ __device__
unsigned int expandBits32(unsigned int v) {
    v = (v * 0x00010001u) & 0xFF0000FFu;
    v = (v * 0x00000101u) & 0x0F00F00Fu;
    v = (v * 0x00000011u) & 0xC30C30C3u;
    v = (v * 0x00000005u) & 0x49249249u;
    return v;
}

__host__ __device__
unsigned int mortonEncode32(unsigned int x, unsigned int y, unsigned int z) {
    unsigned int xx = expandBits32(x);
    unsigned int yy = expandBits32(y);
    unsigned int zz = expandBits32(z);
    return xx * 4 + yy * 2 + zz;
}
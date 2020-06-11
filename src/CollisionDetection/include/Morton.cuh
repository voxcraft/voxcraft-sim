//
// Created by David Matthews on 5/21/20.
//

#ifndef COLLISIONDETECTION_MORTON_CUH
#define COLLISIONDETECTION_MORTON_CUH

#include <iostream>

/**
 * Functor which computes morton numbers from the x, y, z ranks of objects.
 */
struct init_morton_func {
    unsigned int *xRanks;
    unsigned int *yRanks;
    unsigned int *zRanks;
    unsigned long long int *mPoints;
    unsigned int *mIds;

    __device__ __host__
    init_morton_func() :
            init_morton_func(nullptr, nullptr, nullptr, nullptr, nullptr) {}

    __device__ __host__
    init_morton_func(unsigned int *xRanks, unsigned int *yRanks, unsigned int *zRanks, unsigned long long int *mPts,
                     unsigned int *mIds) :
            xRanks(xRanks), yRanks(yRanks), zRanks(zRanks), mPoints(mPts), mIds(mIds) {}

    /**
     * Computes morton number for object at index idx.
     * @param idx
     */
    __host__ __device__
    void operator()(unsigned int idx);
};

/**
 * Functor which computes morton numbers from the x, y, z ranks of objects.
 */
struct init_morton_func_fast {
    float xmin, ymin, zmin, xrange, yrange, zrange, expansionRange, *x, *y, *z;

    unsigned long long int *mPoints;
    unsigned int *mIds;

    __device__ __host__
    init_morton_func_fast() :
            init_morton_func_fast({0.f, 0.f}, {0.f, 0.f}, {0.f, 0.f}, nullptr, nullptr, nullptr, nullptr, nullptr) {}

    __device__ __host__
    init_morton_func_fast(float2 xlims, float2 ylims, float2 zlims, float *x, float *y, float *z,
                          unsigned long long int *mPts, unsigned int *mIds) :
            x(x), y(y), z(z), mPoints(mPts), mIds(mIds) {
        expansionRange = 2097152.f; // 1<<21

        xmin = xlims.x;
        xrange = xlims.y - xlims.x;

        ymin = ylims.x;
        yrange = ylims.y - ylims.x;

        zmin = zlims.x;
        zrange = zlims.y - zlims.x;
    }

    /**
     * Computes morton number (based on Karras method) for object at index idx.
     * @param idx
     */
    __host__ __device__
    void operator()(unsigned int idx);
};

/**
 * The Morton struct is primarly used for printing to std::cout via thrust::copy
 */
struct Morton {
    unsigned long long int morton_code;

    __host__ __device__ Morton() :
            Morton(0) {}

    __host__ __device__ Morton(unsigned long long int morton_code) :
            morton_code(morton_code) {}
};

/**
 * Prints a binary representation of a morton number (e.g. an unsigned int)
 *
 * @param out stream to write to
 * @param m morton number to write to to the outstream
 * @return the out stream
 */
std::ostream &operator<<(std::ostream &out, const Morton &m);

/**
 * Functor which transforms an array of mortons stored as unsigned ints to an array of Morton Structs for easier printing.
 */
struct morton_func {
    unsigned long long int *mortonA;
    Morton *mortonB;

    morton_func(unsigned long long int *mortonA, Morton *mortonB) :
            mortonA(mortonA), mortonB(mortonB) {}

    /**
     * Transforms unsigned int morton array to an array of morton structs for printing.
     * @param idx Array index to transform
     */
    __device__ __host__
    void operator()(int idx);

};

// Helper functions to map ints to morton numbers.


__host__ __device__
unsigned long long int expandBits64(unsigned long long int a);

__host__ __device__
unsigned long long int mortonEncode64(unsigned long long int x, unsigned long long int y, unsigned long long int z);

__host__ __device__
unsigned int expandBits32(unsigned int v);

__host__ __device__
unsigned int mortonEncode32(unsigned int x, unsigned int y, unsigned int z);

#endif //COLLISIONDETECTION_MORTON_CUH

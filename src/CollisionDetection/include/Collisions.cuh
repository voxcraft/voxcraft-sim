//
// Created by David Matthews on 5/14/20.
//

#ifndef COLLISIONDETECTION_COLLISIONS_CUH
#define COLLISIONDETECTION_COLLISIONS_CUH

#include <thrust/functional.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/random.h>
#include <stdio.h>

#include "BoundingBox.cuh"
#include "Morton.cuh"

/**
 * Collision struct used for sorting collions found on GPU / CPU for unit testing.
 */
struct Collision {
    unsigned int a, b;

    __host__ __device__ Collision() :
            Collision(0, 0) {}

    __host__ __device__ Collision(unsigned int a, unsigned int b) :
            a(a), b(b) {}
};


/**
 * Given an internal node (as an index), computes the range of leaf nodes which are covered by this node.
 * First step in building the BVH tree (next is finding split).
 *
 * @param sortedMortonCodes Sorted list of morton codes
 * @param numObjects Number of objects in the system
 * @param idx Internal node index.
 * @return tuple of (first, last) index of the range that this node covers.
 */
__device__
int2 determineRange(unsigned long long int *sortedMortonCodes, unsigned int numObjects, unsigned int idx);

/**
 * Given a range of leaf node indexes, and the corresponding morton numbers, computes the split index where
 * the most significant bit of the morton numbers changes.
 *
 * @param sortedMortonCodes The sorted morton codes.
 * @param first start index of the morton code range of interest
 * @param last end index of the morton code range of interest
 * @return the split index.
 */
__device__
int findSplit(unsigned long long int *sortedMortonCodes, unsigned int first, unsigned int last);

/**
 * Functor which maps over all internal nodes in the BVH tree, computes their ranges and splits, saving details about
 * the parent / child relations between each node.
 *
 * Root node is represented with 0xFFFFFFFF id.
 */
struct build_bvh_tree_func {
    unsigned int N;
    unsigned long long int *sortedMortonCodes;
    unsigned int *leaf_parents, *internal_parents, *internal_childA, *internal_childB;

    build_bvh_tree_func() :
            build_bvh_tree_func(0, nullptr, nullptr, nullptr, nullptr, nullptr) {}

    build_bvh_tree_func(unsigned int N, unsigned long long int *sortedMortonCodes, unsigned int *leaf_parents,
                        unsigned int *internal_parents, unsigned int *internal_childA, unsigned int *internal_childB) :
            N(N), sortedMortonCodes(sortedMortonCodes), leaf_parents(leaf_parents), internal_parents(internal_parents),
            internal_childA(internal_childA), internal_childB(internal_childB) {}

    __device__
    void operator()(unsigned int idx);

};

/**
 * Functor which traverses the BVH tree from leaf nodes up updating the bounding volumes that each node covers.
 */
struct fill_bvh_tree_with_bounding_boxes_func {
    unsigned int N;
    BoundingBox *bboxTree;
    float *xPts, *yPts, *zPts, *radius;
    unsigned int *sortedMortonIds;
    unsigned int *leaf_parents, *internal_parents, *internal_childA, *internal_childB, *bbox_complete_flag;

    fill_bvh_tree_with_bounding_boxes_func() :
            fill_bvh_tree_with_bounding_boxes_func(0,
                    nullptr,
                    nullptr,
                    nullptr,
                    nullptr,
                    nullptr,
                    nullptr,
                    nullptr,
                    nullptr,
                    nullptr,
                    nullptr,
                    nullptr) {}

    fill_bvh_tree_with_bounding_boxes_func(unsigned int N, BoundingBox *bboxTree, float *xPts, float *yPts, float *zPts,
                                           float *radius, unsigned int *sortedMortonIds, unsigned int *leaf_parents,
                                           unsigned int *internal_parents, unsigned int *internal_childA,
                                           unsigned int *internal_childB, unsigned int *bbox_complete_flag) :
            N(N), bboxTree(bboxTree), xPts(xPts), yPts(yPts), zPts(zPts), radius(radius),
            sortedMortonIds(sortedMortonIds), leaf_parents(leaf_parents), internal_parents(internal_parents),
            internal_childA(internal_childA), internal_childB(internal_childB),
            bbox_complete_flag(bbox_complete_flag) {}

    __device__
    unsigned int operator()(int leafIdx);
};

/**
 * Functor which traverses the BVH tree from each leaf node to identify potential collisions where
 * the bounding volumes of pairs of leaf nodes overlap
 */
struct find_potential_collisions_func {
    unsigned int N, NUM_INTERNAL, MAX_NUM_COLLISIONS;
    BoundingBox *boundingBoxTree;
    unsigned int *sortedMortonIds, *internalChildrenA, *internalChildrenB, *potentialCollisionIdx;
    Collision *potentialCollisionList;
    float *x, *y, *z, *r;

    find_potential_collisions_func() :
            find_potential_collisions_func(0,
                    0,
                    nullptr,
                    nullptr,
                    nullptr,
                    nullptr,
                    nullptr,
                    nullptr,
                    nullptr,
                    nullptr,
                    nullptr,
                    nullptr) {}

    find_potential_collisions_func(unsigned int N, unsigned int max_cols, unsigned int *sortedMortonIds,
                                   BoundingBox *boundingBoxTree, unsigned int *internalChildrenA,
                                   unsigned int *internalChildrenB, unsigned int *potentialCollisionIdx,
                                   Collision *potentialCollisionList, float *x, float *y, float *z, float *r) :
            N(N), sortedMortonIds(sortedMortonIds), boundingBoxTree(boundingBoxTree),
            internalChildrenA(internalChildrenA), internalChildrenB(internalChildrenB),
            potentialCollisionIdx(potentialCollisionIdx), potentialCollisionList(potentialCollisionList), x(x), y(y),
            z(z), r(r) {
        NUM_INTERNAL = N - 1;
        MAX_NUM_COLLISIONS = max_cols;
    }

    /**
     * finds and saves all potential collisions of sphere at leaf node idx with all other nodes.
     * @param idx leaf node to check
     */
    __device__
    void operator()(unsigned int idx);
};

/**
 * Functor which checks if each potential collision is actually a collision
 */
struct check_potential_collisions_func {
    float *x, *y, *z, *r;

    check_potential_collisions_func() :
            check_potential_collisions_func(nullptr, nullptr, nullptr, nullptr) {}

    check_potential_collisions_func(float *x, float *y, float *z, float *r) :
            x(x), y(y), z(z), r(r) {}

    __device__
    bool operator()(const Collision &c);
};

struct check_potential_collisions_N2_func : public thrust::binary_function<unsigned int, unsigned int, bool> {
    unsigned int MAX_COLS, *cIdx;
    float *x, *y, *z, *r;
    Collision *c;

    __host__ __device__
    check_potential_collisions_N2_func() :
            check_potential_collisions_N2_func(0, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr) {}

    __host__ __device__
    check_potential_collisions_N2_func(unsigned int MAX_COLS, unsigned int *cIdx, float *x, float *y, float *z,
                                       float *r, Collision *c) :
            MAX_COLS(MAX_COLS), cIdx(cIdx), x(x), y(y), z(z), r(r), c(c) {}

    template<typename T>
    __device__
    bool operator()(T &t) {
        unsigned int colAId = thrust::get<0>(t);
        unsigned int colBId = thrust::get<1>(t);
        if (colAId >= colBId) { // only consider collisions where a < b to avoid duplicates.
            return false;
        }
        if ((x[colAId] - x[colBId]) * (x[colAId] - x[colBId]) + (y[colAId] - y[colBId]) * (y[colAId] - y[colBId]) +
            (z[colAId] - z[colBId]) * (z[colAId] - z[colBId]) < (r[colAId] + r[colBId]) * (r[colAId] + r[colBId])) {
            unsigned int currCollisionIdx = atomicInc(cIdx, 0xFFFFFFFF);
            if (currCollisionIdx >= MAX_COLS) {
                return false;
            }
            c[currCollisionIdx].a = colAId;
            c[currCollisionIdx].b = colBId;
            return true;
        }
        return false;
    }
};


__host__ __device__
bool operator<(const Collision &lhs, const Collision &rhs);

__host__ __device__
bool operator<=(const Collision &lhs, const Collision &rhs);

__host__ __device__
bool operator>(const Collision &lhs, const Collision &rhs);

__host__ __device__
bool operator>=(const Collision &lhs, const Collision &rhs);

__host__ __device__
bool operator==(const Collision &lhs, const Collision &rhs);

__host__
std::ostream &operator<<(std::ostream &out, const Collision &c);

/**
 * functor used for merging two arrays where each index represents a pair of objects colliding.
 */
struct zip_collisions_func {
    Collision *c;
    unsigned int *a, *b;

    zip_collisions_func(Collision *c, unsigned int *a, unsigned int *b) :
            a(a), b(b), c(c) {}

    __host__ __device__
    void operator()(unsigned int idx);
};

unsigned int find_collisions_cpu(unsigned int N, float *x, float *y, float *z, float *r, Collision *cols);

#endif //COLLISIONDETECTION_COLLISIONS_CUH

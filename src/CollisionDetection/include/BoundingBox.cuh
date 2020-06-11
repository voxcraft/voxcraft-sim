//
// Created by David Matthews on 5/21/20.
//

#ifndef COLLISIONDETECTION_BOUNDINGBOX_CUH
#define COLLISIONDETECTION_BOUNDINGBOX_CUH

#include <iostream>

struct BoundingBox {
    float x_min, x_max, y_min, y_max, z_min, z_max;

    __host__ __device__ BoundingBox() :
            BoundingBox(0, 0, 0, 0, 0, 0) {}

    __host__ __device__ BoundingBox(float x_min, float x_max, float y_min, float y_max, float z_min, float z_max) :
            x_min(x_min), x_max(x_max), y_min(y_min), y_max(y_max), z_min(z_min), z_max(z_max) {}

};

std::ostream &operator<<(std::ostream &out, const BoundingBox &bb);

/**
 * Check if bounding boxes a, b overlap. To do this, check if they do not overlap in all x, y, z
 * @param a Bounding box a
 * @param b Bounding box b
 * @return True if Bounding boxes overlap, False otherwise.
 */
__device__ __host__
bool overlap(BoundingBox *a, BoundingBox *b);

#endif //COLLISIONDETECTION_BOUNDINGBOX_CUH

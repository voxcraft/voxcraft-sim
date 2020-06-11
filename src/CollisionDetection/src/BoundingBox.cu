//
// Created by David Matthews on 5/21/20.
//

#include "../include/BoundingBox.cuh"
#include <iomanip>

std::ostream &operator<<(std::ostream &out, const BoundingBox &bb) {
    return out << "((" << std::setprecision(3) << std::setw(8) << bb.x_min << " <-> " << std::setprecision(3)
               << std::setw(8) << bb.x_max << "), (" << std::setprecision(3) << std::setw(8) << bb.y_min << " <-> "
               << std::setprecision(3) << std::setw(8) << bb.y_max << "), (" << std::setprecision(3) << std::setw(8)
               << bb.z_min << " <-> " << std::setprecision(3) << std::setw(8) << bb.z_max << "))";
}

__device__ __host__
bool overlap(BoundingBox *a, BoundingBox *b) {

    return !((a->x_min > b->x_max || b->x_min > a->x_max) || // do a and b not overlap in X
             (a->y_min > b->y_max || b->y_min > a->y_max) || // do a and b not overlap in y
             (a->z_min > b->z_max || b->z_min > a->z_max)); // do a and b not overlap in z
}

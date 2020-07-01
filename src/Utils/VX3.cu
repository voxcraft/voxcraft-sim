#include "VX3.cuh"

/**
 * Convert a 3D position to 1D offset
 * @param pos A 3D position
 * @param dim The dimension of the buffer
 * @return 1D offset for the buffer
 */
__device__ int to1D(VX3_Vec3D<int> pos, VX3_Vec3D<int>dim) {
    if (pos.x < 0 || pos.y < 0 || pos.z < 0 || pos.x >= dim.x || pos.y >= dim.y || pos.z >= dim.z)
        return -1; // Overflow
    return dim.y * dim.z * pos.x + dim.z * pos.y + pos.z;
}
/**
 * Convert a 1D offset to 3D position
 * @param offset A 1D offset
 * @param dim The dimension of the buffer
 * @return A 3D position
 */
__device__ VX3_Vec3D<int> to3D(int offset, VX3_Vec3D<int>dim) {
    int x, y, z;
    z = offset % dim.z;
    int residual;
    residual = (int)((offset - z) / dim.z);
    y = residual % dim.y;
    residual = (int)((residual - y) / dim.y);
    x = residual;
    return VX3_Vec3D<int>(x, y, z);
}
#include "VX3_Quat3D.h"
#include "math_constants.h"

__device__ VX3_Vec3D<double> linkDirectionVec3D(linkDirection linkdir) {
    switch (linkdir) {
    case X_POS:
        return VX3_Vec3D<>(1, 0, 0);
    case X_NEG:
        return VX3_Vec3D<>(-1, 0, 0);
    case Y_POS:
        return VX3_Vec3D<>(0, 1, 0);
    case Y_NEG:
        return VX3_Vec3D<>(0, -1, 0);
    case Z_POS:
        return VX3_Vec3D<>(0, 0, 1);
    case Z_NEG:
        return VX3_Vec3D<>(0, 0, -1);
    default:
        printf("ERROR.\n");
    }
}

__global__ void kernel() {
    VX3_Vec3D<> pos = VX3_Vec3D<>(0,1,0);
    VX3_Quat3D<> q;
    q.FromAngleToPosX(pos);
    q.debug();

    VX3_Vec3D<> target = VX3_Vec3D<>(1,1,0);
    target = q.RotateVec3D(target);
    target.debug();

    printf("\n");
}

int main() {
    kernel<<<1, 1>>>();
    cudaDeviceSynchronize();
}
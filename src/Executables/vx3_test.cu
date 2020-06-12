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

    VX3_Quat3D<double> orientation, o1;
    VX3_Vec3D<> pos, pos1, pos2, pos3;
    double n = sqrt(1-0.98*0.98-0.15*0.15);
    orientation = VX3_Quat3D<double>(0.98, 0, 0.15, n);
    pos = linkDirectionVec3D(Z_NEG);
    pos.debug();
    
    pos1 = orientation.RotateVec3D(pos);
    pos1.debug();

    pos2 =orientation.RotateVec3DInv(pos1);
    pos2.debug();

    printf("\n");
}

int main() {
    kernel<<<1, 1>>>();
    cudaDeviceSynchronize();
}
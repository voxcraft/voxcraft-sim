#include "VX3_SignalDiffusion.h"

__device__ void VX3_SignalDiffusion::deviceInit(VX3_Vec3D<int> dim) {
    dimension = dim;
    map = (double *)malloc(dimension.x * dimension.y * dimension.z * sizeof(double));
    map_nextstep = (double *)malloc(dimension.x * dimension.y * dimension.z * sizeof(double));
    if (map==NULL || map_nextstep==NULL) {
        printf("Not enough memory.\n");
    }
}

__global__ void gpu_diffusion(double dt_diffusivity_distance, int dim_x, int dim_y, int dim_z, double *map, double *map_nextstep) {
    // Formula, analogous to 2D Cellular Automata version
    // Approximated from Fick's Law
    // new amount = (1-D*dt)*(old amount) + sum_i6(D*dt*(neighbor[i] amount))
    VX3_Vec3D<int> dim = VX3_Vec3D<int>(dim_x, dim_y, dim_z);
    VX3_Vec3D<int> pos = VX3_Vec3D<int>(
        threadIdx.x + blockDim.x * blockIdx.x, 
        threadIdx.y + blockDim.y * blockIdx.y,
        threadIdx.z + blockDim.z * blockIdx.z);
    int offset;
    double diffuse_in = 0;
    for (int i = 0; i < 6; i++) {
        VX3_Vec3D<int> neighbor_pos = pos;
        switch ((linkDirection)i) {
        case X_POS:
            neighbor_pos.x++;
            break;
        case X_NEG:
            neighbor_pos.x--;
            break;
        case Y_POS:
            neighbor_pos.y++;
            break;
        case Y_NEG:
            neighbor_pos.y--;
            break;
        case Z_POS:
            neighbor_pos.z++;
            break;
        default:
            neighbor_pos.z--;
            break;
        }
        offset = to1D(neighbor_pos, dim);
        if (offset != -1) { // within the space defined. assuming things diffused out of the edges will never come back.
            diffuse_in += *(map + offset) * dt_diffusivity_distance;
        }
    }
    offset = to1D(pos, dim);
    *(map_nextstep + offset) = (1 - dt_diffusivity_distance) * (*(map + offset)) + diffuse_in / 6;
}

__device__ void VX3_SignalDiffusion::doTimestep(double dt) {
    gpu_diffusion<<<dimension.x, dimension.y, dimension.z>>>(dt * diffusivity / distance, dimension.x, dimension.y, dimension.z, map, map_nextstep);
    CUDA_CHECK_AFTER_CALL();
    VcudaDeviceSynchronize();
    // Switch the pointers to the buffers
    double *tmp = map;
    map = map_nextstep;
    map_nextstep = tmp;
}

__device__ double VX3_SignalDiffusion::quaryQuantityAtPosition(VX3_Vec3D<int> pos) { 
    int offset = to1D(pos, dimension);
    if (offset==-1) {
        return -1;
    }
    return *(map + offset);
}

__device__ void VX3_SignalDiffusion::addSolute(double amount, VX3_Vec3D<int> pos) {
    int offset = to1D(pos, dimension);
    if (offset==-1) {
        printf("VX3_SignalDiffusion::addSolute Overflow.\n");
    }
    *(map + offset) += amount;
}

__device__ void VX3_SignalDiffusion::removeSolute(double amount, VX3_Vec3D<int> pos) {
    int offset = to1D(pos, dimension);
    if (offset==-1) {
        printf("VX3_SignalDiffusion::addSolute Overflow.\n");
    }
    double *p = map + offset;
    if (*p > amount) {
        *p -= amount;
    } else {
        *p = 0;
    }
}

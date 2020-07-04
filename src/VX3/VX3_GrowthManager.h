//
// Created by Sida Liu
//  This class implements the growth of a body.
//
#if !defined(GROWTH_MANAGER_H)
#define GROWTH_MANAGER_H

#include "VX3.cuh"

class VX3_VoxelyzeKernel;
class VX3_Voxel;

class VX3_GrowthManager
{
public:
    VX3_VoxelyzeKernel* d_kernel;
    // int growMutex = 0;

    __device__ VX3_GrowthManager(VX3_VoxelyzeKernel* k);
    __device__ bool grow(VX3_Voxel* vx, int direction, VX3_Vec3D<> new_position);
};

#endif // GROWTH_MANAGER_H

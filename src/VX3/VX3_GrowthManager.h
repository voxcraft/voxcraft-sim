//
// Created by Sida Liu
//  This class implements the growth of a body.
//
#if !defined(GROWTH_MANAGER_H)
#define GROWTH_MANAGER_H

class VX3_VoxelyzeKernel;

class VX3_GrowthManager
{
public:
    VX3_VoxelyzeKernel* d_kernel;

    __device__ VX3_GrowthManager(VX3_VoxelyzeKernel* k);
    __device__ bool grow();
};

#endif // GROWTH_MANAGER_H

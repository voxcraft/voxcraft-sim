#include "VX3_GrowthManager.h"
#include "VX3_VoxelyzeKernel.cuh"

__device__ VX3_GrowthManager::VX3_GrowthManager(VX3_VoxelyzeKernel *k) { d_kernel = k; }

__device__ bool VX3_GrowthManager::grow(VX3_Voxel* v, int direction, VX3_Vec3D<> new_position) {
    VX3_Voxel *new_voxel;
    // TODO: need check surrounding of the new voxel
    // get a new voxel's pointer from pre-allocated memory
    int current_num_voxels = atomicAdd( &(d_kernel->num_d_voxels), 1);
    if (current_num_voxels - d_kernel->num_d_init_voxels >= 10000) { // memory limitation, refer to pre-allocation.
        return false;
    }
    new_voxel = &d_kernel->d_voxels[current_num_voxels];
    new_voxel->deviceInit(d_kernel);
    new_voxel->mat = &d_kernel->d_voxelMats[0];
    new_voxel->pos = new_position;
    new_voxel->orient = v->orient;
    // d_kernel->d_attach_manager->attachForNewVoxel(v, available_direction, new_voxel, oppositeDirection(available_direction));
    d_kernel->isSurfaceChanged = true;
    return true;
}

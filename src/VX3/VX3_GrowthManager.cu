#include "VX3_GrowthManager.h"
#include "VX3_VoxelyzeKernel.cuh"

__device__ VX3_GrowthManager::VX3_GrowthManager(VX3_VoxelyzeKernel *k) { d_kernel = k; }

__device__ bool VX3_GrowthManager::grow() {
    // randomly pick one surface voxel
    // check its surround
    // add a new voxel to proper position
    if (d_kernel->num_d_voxels - d_kernel->num_d_init_voxels < 1000) { // memory limitation, refer to pre-allocation.
        int r = d_kernel->randomGenerator->randint(d_kernel->num_d_surface_voxels);
        DEBUG_PRINT("r: %d, surface_voxels: %d.\n", r, d_kernel->num_d_surface_voxels);
        VX3_Voxel *v = d_kernel->d_surface_voxels[r];
        int available_direction;
        VX3_Link *n = (VX3_Link *)1;
        while (n) {
            available_direction = d_kernel->randomGenerator->randint(6);
            n = v->links[available_direction];
        }
        VX3_Vec3D<> new_position = VX3_Vec3D<>();
        switch ((linkDirection)available_direction) {
        case X_POS:
            new_position.x += d_kernel->voxSize;
            break;
        case X_NEG:
            new_position.x -= d_kernel->voxSize;
            break;
        case Y_POS:
            new_position.y += d_kernel->voxSize;
            break;
        case Y_NEG:
            new_position.y -= d_kernel->voxSize;
            break;
        case Z_POS:
            new_position.z += d_kernel->voxSize;
            break;
        case Z_NEG:
            new_position.z -= d_kernel->voxSize;
        }
        VX3_Voxel *new_voxel = &d_kernel->d_voxels[d_kernel->num_d_voxels];
        // need orientation
        new_voxel->pos = v->pos + v->orient.RotateVec3DInv(new_position);
        new_voxel->orient = v->orient;
        // printf("d_voxels[%d].pos: %e,%e,%e. === %e.\n", d_kernel->num_d_voxels, new_voxel->pos.x, new_voxel->pos.y, new_voxel->pos.z, d_kernel->voxSize / 2);
        if (new_voxel->pos.z < d_kernel->voxSize / 2) {
            // printf("no.\n");
            return false;
        }
        // printf("yes.\n");
        // need check surrounding
        new_voxel->pos = new_position;
        new_voxel->mat = &d_kernel->d_voxelMats[0];
        new_voxel->deviceInit(d_kernel);
        new_voxel->enableFloor(true);
        new_voxel->updateGroup();
        d_kernel->num_d_voxels++;
        
        d_kernel->isSurfaceChanged = true;
        return true;
    }
    return false;
}

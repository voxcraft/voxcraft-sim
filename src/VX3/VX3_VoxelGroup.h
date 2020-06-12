#if !defined(VX3_VOXELGROUP_H)
#define VX3_VOXELGROUP_H

#include "VX3_vector.cuh"
class VX3_Voxel;

class VX3_VoxelGroup {
public:
    /* data */
    int dim_x, dim_y, dim_z;
    int min_x, min_y, min_z;
    int max_x, max_y, max_z;

    VX3_Voxel **d_group_map = NULL;
    VX3_dVector<VX3_Voxel *> d_surface_voxels;
    VX3_dVector<VX3_Voxel *> d_voxels;
    bool needRebuild = true;

    __device__ VX3_VoxelGroup(/* args */);
    __device__ void updateGroup(VX3_Voxel *voxel); // Update all the group info that voxel is in. BFS.
    __device__ int to1D(int x, int y, int z);
    __device__ void to3D(int offset, int *ret_x, int *ret_y, int *ret_z);
    __device__ int getVoxelOffset(VX3_Voxel *voxel);
    __device__ bool isCompatible(VX3_Voxel *voxel_host, VX3_Voxel *voxel_remote, int* ret_linkdir_1, int* ret_linkdir_2);
    __device__ void reorient_lattice();

};

#endif // VX3_VOXELGROUP_H

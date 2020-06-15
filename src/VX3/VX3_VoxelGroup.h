#if !defined(VX3_VOXELGROUP_H)
#define VX3_VOXELGROUP_H

#include "VX3_vector.cuh"

class VX3_VoxelyzeKernel;
class VX3_Voxel;

class VX3_VoxelGroup {
public:
    /* data */
    VX3_VoxelyzeKernel *d_kernel; // saved pointer to the whole simulation
    int dim_x, dim_y, dim_z; // dimension of this group

    VX3_Voxel **d_group_map = NULL; // store 3-dimensional voxel pointers
    VX3_dVector<VX3_Voxel *> d_surface_voxels; // all surface voxels in this group
    VX3_dVector<VX3_Voxel *> d_voxels; // all voxels in this group
    bool needRebuild = true;
    int buildMutex = 0; // only allow one call to update(build) this group

    int hasNewLink = 0; // how many new links in this group.

    __device__ VX3_VoxelGroup(VX3_VoxelyzeKernel *k);
    __device__ ~VX3_VoxelGroup();
    __device__ VX3_Vec3D<int> moveGroupPosition(VX3_Vec3D<int> from, linkDirection dir, int step = 1); // return the step next position in the group
    __device__ void updateGroup(VX3_Voxel *voxel); // Update all the group info that voxel is in. BFS.
    __device__ int to1D(VX3_Vec3D<int> groupPosition); // for generating an index(offset) for `d_group_map` from groupPosition
    __device__ VX3_Vec3D<int> to3D(int offset); // get groupPosition back from an index(offset)
    __device__ bool isCompatible(VX3_Voxel *voxel_host, VX3_Voxel *voxel_remote, int *ret_linkdir_1, int *ret_linkdir_2);
};

#endif // VX3_VOXELGROUP_H

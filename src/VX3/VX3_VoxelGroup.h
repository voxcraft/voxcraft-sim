//
// Created by Sida Liu
//  VoxelGroup is a group of connected voxels, you can also call it "body".
//
#if !defined(VX3_VOXELGROUP_H)
#define VX3_VOXELGROUP_H

#include "VX3_vector.cuh"

class VX3_VoxelyzeKernel;
class VX3_Voxel;

class VX3_VoxelGroup {
public:
    /* data */
    bool removed;
    VX3_VoxelyzeKernel *d_kernel; // saved pointer to the whole simulation
    int dim_x, dim_y, dim_z; // dimension of this group

    VX3_Voxel **d_group_map; // store 3-dimensional voxel pointers
    int buffer_size_group_map;
    VX3_dVector<VX3_Voxel *> d_surface_voxels; // all surface voxels in this group
    VX3_dVector<VX3_Voxel *> d_voxels; // all voxels in this group

    int hasNewLink; // how many new links in this group.
    int needUpdate; // avoid multiple attachments before update the group

    __device__ VX3_VoxelGroup(VX3_VoxelyzeKernel *k);
    __device__ void deviceInit(VX3_VoxelyzeKernel *k);
    __device__ void switchAllVoxelsTo(VX3_VoxelGroup* group);
    __device__ VX3_Vec3D<int> moveGroupPosition(VX3_Vec3D<int> from, linkDirection dir, int step = 1); // return the step next position in the group
    
    // These two methods need to consider the racing condition:
    __device__ void updateGroup(VX3_Voxel *start_voxel); // Update all the group info of voxels that is in this group, start from d_voxels[0]. BFS.
    __device__ bool isCompatible(VX3_Voxel *voxel_host, VX3_Voxel *voxel_remote, int *ret_linkdir_1, int *ret_linkdir_2); // Check host and remote group are compatible for attachment.

    // No need to rebuild in the same timestep.
    double lastBuildTime; // time of the simulation in seconds

};

#endif // VX3_VOXELGROUP_H

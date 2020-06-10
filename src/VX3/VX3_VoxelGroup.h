#if !defined(VX3_VOXELGROUP_H)
#define VX3_VOXELGROUP_H

class VX3_Voxel;

class VX3_VoxelGroup {
public:
    /* data */
    int x, y, z;
    VX3_Voxel *d_group_map = NULL;
    bool needRebuild = true;

    __device__ VX3_VoxelGroup(/* args */);
    __device__ void updateGroup(VX3_Voxel *voxel); // Update all the group info that voxel is in. BFS.
};

#endif // VX3_VOXELGROUP_H

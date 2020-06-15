#include "VX3_VoxelGroup.h"

#include "VX3_Voxel.h"

__device__ VX3_VoxelGroup::VX3_VoxelGroup(VX3_VoxelyzeKernel* k) {
    d_kernel = k;    
}
__device__ VX3_VoxelGroup::~VX3_VoxelGroup() {
    // if someone is deleting this group, set the d_group of all voxels in this group to NULL
    for (int i=0;i<d_voxels.size();i++) {
        if (d_voxels[i]->d_group == this) {
            d_voxels[i]->d_group = NULL;
        }
    }
}

__device__ VX3_Vec3D<int> VX3_VoxelGroup::moveGroupPosition(VX3_Vec3D<int> from, linkDirection dir, int step) {
    VX3_Vec3D<int> to = from;
    switch (dir) {
    case X_POS:
        to.x += step;
        break;
    case X_NEG:
        to.x -= step;
        break;
    case Y_POS:
        to.y += step;
        break;
    case Y_NEG:
        to.y -= step;
        break;
    case Z_POS:
        to.z += step;
        break;
    case Z_NEG:
        to.z -= step;
        break;
    default:
        printf("ERROR in moveGroupPosition.\n");
    }
    return to;
}

__device__ void VX3_VoxelGroup::updateGroup(VX3_Voxel *voxel) {
    if (atomicCAS(&buildMutex, 0, 1)==0) {
        // only allow one call to update(build) this group ( might be call simultaneously from VX3_Voxels )
        if (needRebuild) {
            needRebuild = false;
            // First set *voxel as origin, and negative number is allowed.
            // After everything is mapped, change origin to (min_x, min_y, min_z).
            voxel->groupPosition = VX3_Vec3D<int>(0, 0, 0);

            // Rebuild a map of Group
            //// BF Search for all neighbors
            VX3_dDictionary<VX3_Voxel *, int> BFS_visited;
            VX3_dQueue<VX3_Voxel *> BFS_Queue;
            VX3_dVector<VX3_Voxel *> BFS_result;

            // printf("Start from (%d, %d, %d).\n", voxel->groupPosition.x, voxel->groupPosition.y, voxel->groupPosition.z);
            BFS_result.push_back(voxel);
            BFS_visited.set(voxel, 1);
            BFS_Queue.push_back(voxel);

            while (!BFS_Queue.isEmpty()) {
                VX3_Voxel *v = BFS_Queue.pop_front();
                for (int i = 0; i < 6; i++) {
                    VX3_Voxel *neighbor = v->adjacentVoxel((linkDirection)i);
                    if (neighbor) {
                        if (BFS_visited.get(neighbor) == -1) {
                            neighbor->groupPosition = moveGroupPosition(v->groupPosition, (linkDirection)i);
                            // printf("Visit (%d, %d, %d).\n", neighbor->groupPosition.x, neighbor->groupPosition.y, neighbor->groupPosition.z);
                            // Set all connected voxels' d_group to this
                            neighbor->switchGroupTo(this);

                            BFS_result.push_back(neighbor);
                            BFS_visited.set(neighbor, 1);
                            BFS_Queue.push_back(neighbor);
                        }
                    }
                }
            }

            //// Allocate memory for the map
            int min_x, min_y, min_z, max_x, max_y, max_z;
            min_x = INT_MAX;
            min_y = INT_MAX;
            min_z = INT_MAX;
            max_x = 0;
            max_y = 0;
            max_z = 0;
            for (int i = 0; i < BFS_result.size(); i++) {
                VX3_Voxel *v = BFS_result[i];
                if (v->groupPosition.x < min_x)
                    min_x = v->groupPosition.x;
                if (v->groupPosition.y < min_y)
                    min_y = v->groupPosition.y;
                if (v->groupPosition.z < min_z)
                    min_z = v->groupPosition.z;

                if (v->groupPosition.x > max_x)
                    max_x = v->groupPosition.x;
                if (v->groupPosition.y > max_y)
                    max_y = v->groupPosition.y;
                if (v->groupPosition.z > max_z)
                    max_z = v->groupPosition.z;
            }
            dim_x = max_x - min_x + 1;
            dim_y = max_y - min_y + 1;
            dim_z = max_z - min_z + 1;

            if (d_group_map) {
                free(d_group_map);
                d_group_map = NULL;
            }
            d_group_map = (VX3_Voxel **)malloc(dim_x * dim_y * dim_z * sizeof(VX3_Voxel *));
            for (int i = 0; i < dim_x * dim_y * dim_z; i++) {
                d_group_map[i] = NULL;
            }
            d_voxels.clear();
            d_surface_voxels.clear();
            for (int i = 0; i < BFS_result.size(); i++) {
                VX3_Voxel *v = BFS_result[i];
                v->groupPosition.x -= min_x;
                v->groupPosition.y -= min_y;
                v->groupPosition.z -= min_z;

                int offset = to1D(v->groupPosition);
                d_group_map[offset] = v;
                d_voxels.push_back(v);
                for (int j = 0; j < 6; j++) {
                    if (v->links[j] == NULL) {
                        d_surface_voxels.push_back(v);
                        break;
                    }
                }
            }
        }
        atomicExch(&buildMutex, 0);
    }
}

__device__ int VX3_VoxelGroup::to1D(VX3_Vec3D<int> groupPosition) {
    int x = groupPosition.x;
    int y = groupPosition.y;
    int z = groupPosition.z;
    if (x < 0 || y < 0 || z < 0 || x >= dim_x || y >= dim_y || z >= dim_z)
        return -1; // Overflow
    int offset;
    offset = dim_y * dim_z * x + dim_z * y + z;
    return offset;
}
__device__ VX3_Vec3D<int> VX3_VoxelGroup::to3D(int offset) {
    int x, y, z;
    z = offset % dim_z;
    int residual;
    residual = (int)((offset - z) / dim_z);
    y = residual % dim_y;
    residual = (int)((residual - y) / dim_y);
    x = residual;
    return VX3_Vec3D<int>(x,y,z);
}

__device__ bool VX3_VoxelGroup::isCompatible(VX3_Voxel *voxel_host, VX3_Voxel *voxel_remote, int *ret_linkdir_1, int *ret_linkdir_2) {
    if (voxel_host->d_group != this) {
        printf("This method should be call from voxel_host's d_group.\n");
        return false;
    }
    // Given two voxels, determine the best way to attach them.
    VX3_Vec3D<int> offset_of_link = VX3_Vec3D<int>(0,0,0);
    int potential_link_1, potential_link_2;
    VX3_Quat3D<double> relativeRotation = voxel_remote->orientation().Conjugate() * voxel_host->orientation();
    if (relativeRotation.w > 0.96) // within 15 degree
    {
        VX3_Vec3D<> raw_pos = voxel_remote->position() - voxel_host->position();
        VX3_Vec3D<> pos = voxel_host->orientation().RotateVec3DInv(raw_pos); // the position of remote voxel relative to host voxel.
        if (abs(pos.x) > abs(pos.y) && abs(pos.x) > abs(pos.z)) {
            potential_link_1 = (int)(pos.x > 0 ? X_POS : X_NEG);
            offset_of_link.x = pos.x > 0 ? 1 : -1;
        } else if (abs(pos.y) > abs(pos.z)) {
            potential_link_1 = (int)(pos.y > 0 ? Y_POS : Y_NEG);
            offset_of_link.y = pos.y > 0 ? 1 : -1;
        } else {
            potential_link_1 = (int)(pos.z > 0 ? Z_POS : Z_NEG);
            offset_of_link.z = pos.z > 0 ? 1 : -1;
        }
        potential_link_2 = oppositeDirection(potential_link_1); // only support oppositeDirection attachment for now. Arbitrary attachment is much more difficult.
    } else {
        // Too large of an angle
        return false;
    }
    // Start checking for compatibility.
    // e.g. a 2D example:
    //     Host: 0,0 - 1,0 - 2,0
    //                  ?
    //     Remote:     0,1
    //                  |
    //                 0,0
    // potential_link_1 = Y_NEG
    // position_of_remote_voxel_in_host_group = Voxel(in remote group).groupPosition + remote_diff;
    // remote_diff = Host.groupPosition + offset_of_link - Remote.groupPosition
    // so (0,0) in Remote group becomes (1,-2).

    VX3_Vec3D<int> remote_diff = VX3_Vec3D<int>(0,0,0);
    remote_diff = voxel_host->groupPosition + offset_of_link - voxel_remote->groupPosition;

    for (int i = 0; i < voxel_remote->d_group->d_surface_voxels.size(); i++) {
        VX3_Voxel *v = voxel_remote->d_group->d_surface_voxels[i];
        VX3_Vec3D<int> position_of_remote_voxel_in_host_group = v->groupPosition + remote_diff;
        int offset = to1D(position_of_remote_voxel_in_host_group);
        if (offset == -1) {
            // good, because out of range
        } else if (d_group_map[offset] == NULL) {
            // good, because empty position
        } else {
            // printf("Not Compatible. Offset %d\n", offset);
            return false;
        }
    }
    // printf("Compatible.\n");
    *ret_linkdir_1 = potential_link_1;
    *ret_linkdir_2 = potential_link_2;
    return true;
}

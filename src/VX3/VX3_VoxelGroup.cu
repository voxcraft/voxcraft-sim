#include "VX3_VoxelGroup.h"

#include "VX3_Voxel.h"
#include "VX3_VoxelyzeKernel.cuh"

__device__ VX3_VoxelGroup::VX3_VoxelGroup(VX3_VoxelyzeKernel *k) { d_kernel = k; }

__device__ void VX3_VoxelGroup::deviceInit() {
    d_surface_voxels.clear();
    d_voxels.clear();
}

__device__ void VX3_VoxelGroup::switchAllVoxelsTo(VX3_VoxelGroup *group) {
    if (group == this)
        return;
    for (int i = 0; i < d_voxels.size(); i++) {
        if (d_voxels[i]->d_group == NULL) {
            d_voxels[i]->d_group = group;
        } else if (d_voxels[i]->d_group != NULL && d_voxels[i]->d_group != group) {
            d_voxels[i]->d_group->removed = true;
            d_voxels[i]->d_group = group;
        } else {
            // d_group is already group.
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

__device__ void VX3_VoxelGroup::updateGroup() {
    if (removed)
        return;
    if (lastBuildTime == d_kernel->currentTime)
        return;
    lastBuildTime = d_kernel->currentTime;

    VX3_Voxel *voxel = d_voxels[0];
    int min_x, min_y, min_z, max_x, max_y, max_z;
    min_x = 0;
    min_y = 0;
    min_z = 0;
    max_x = 0;
    max_y = 0;
    max_z = 0;

    // First set *voxel as origin, and negative number is allowed.
    // After everything is mapped, change origin to (min_x, min_y, min_z).
    voxel->groupPosition = VX3_Vec3D<int>(0, 0, 0);

    // Rebuild a map of Group
    //// BF Search for all neighbors
    VX3_dDictionary<VX3_Voxel *, int> BFS_visited;
    VX3_dQueue<VX3_Voxel *> BFS_Queue;
    VX3_dVector<VX3_Voxel *> BFS_result;

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
                    // Set all connected voxels' d_group to this
                    neighbor->d_group->switchAllVoxelsTo(this);

                    BFS_result.push_back(neighbor);
                    BFS_visited.set(neighbor, 1);
                    BFS_Queue.push_back(neighbor);
                    // update the min and max for later use
                    if (neighbor->groupPosition.x < min_x) {
                        min_x = neighbor->groupPosition.x;
                    }
                    if (neighbor->groupPosition.x > max_x) {
                        max_x = neighbor->groupPosition.x;
                    }
                    if (neighbor->groupPosition.y < min_y) {
                        min_y = neighbor->groupPosition.y;
                    }
                    if (neighbor->groupPosition.y > max_y) {
                        max_y = neighbor->groupPosition.y;
                    }
                    if (neighbor->groupPosition.z < min_z) {
                        min_z = neighbor->groupPosition.z;
                    }
                    if (neighbor->groupPosition.z > max_z) {
                        max_z = neighbor->groupPosition.z;
                    }
                } else {
                    // check the group position to make sure the connection is right
                    if (neighbor->groupPosition != moveGroupPosition(v->groupPosition, (linkDirection)i)) {
                        if (v->links[i])
                            v->links[i]->detach(); // TODO: there are potential racing conditions here.
                    }
                }
            }
        }
    }

    //// Allocate memory for the map
    dim_x = max_x - min_x + 1;
    dim_y = max_y - min_y + 1;
    dim_z = max_z - min_z + 1;
    if (buffer_size_group_map < dim_x * dim_y * dim_z) { // size of group map exceed the buffer size
        if (buffer_size_group_map == 0) {
            buffer_size_group_map = (dim_x * dim_y * dim_z >= 16) ? (dim_x * dim_y * dim_z * 2) : 32; // by default, allocate 10, so no need to go through 2,4,8,16,32
        } else {
            buffer_size_group_map = dim_x * dim_y * dim_z * 2; // double the size
        }
        if (d_group_map) {
            free(d_group_map);
            d_group_map = NULL;
        }
        d_group_map = (VX3_Voxel **)malloc(buffer_size_group_map * sizeof(VX3_Voxel *));
        if (!d_group_map) {
            printf("Out of Memory: d_group_map.\n");
        }
    }
    memset(d_group_map, 0, dim_x * dim_y * dim_z * sizeof(VX3_Voxel *)); // only use this much of the buffer
    d_voxels.clear();
    d_surface_voxels.clear();
    for (int i = 0; i < BFS_result.size(); i++) {
        VX3_Voxel *v = BFS_result[i];
        v->groupPosition.x -= min_x;
        v->groupPosition.y -= min_y;
        v->groupPosition.z -= min_z;

        // sam:
        bool absorb = false;
        if (d_kernel->keepJustOneIfManyHaveSameGroupPosition) {

            for (int j = 0; j < d_voxels.size(); j++) {
                if (v->groupPosition == d_voxels[j]->groupPosition) {
                    absorb = true;

                    // Option 1: delete it
                    // v->removed = true;

                    // Option 2: just break off
                    v->d_group = new VX3_VoxelGroup(d_kernel);
                    d_kernel->d_voxel_to_update_group.push_back(v);
                    v->d_group->d_voxels.push_back(v);
                    d_kernel->d_voxelgroups.push_back(v->d_group);

                    // either way, delete all the links
                    for (int k = 0; k < 6; k++) { // check links in all direction
                        if (v->links[k]) {
                            v->links[k]->removed = true;
                            v->adjacentVoxel((linkDirection)k)->links[oppositeDirection(k)] = NULL;
                            v->links[k] = NULL;
                        }
                    }
                }
            }
        }
        if (!absorb) { // sam
            int offset = to1D(v->groupPosition, VX3_Vec3D<int>(dim_x, dim_y, dim_z));
            d_group_map[offset] = v;
            d_voxels.push_back(v);
            // If any link is NULL => is surface voxel
            if (!(v->links[0] && v->links[1] && v->links[2] && v->links[3] && v->links[4] && v->links[5])) {
                d_surface_voxels.push_back(v);
            }
        } // sam
    }
}

////////////////////////////////////////////////////////////////////
// NOTE: There is another racing condition that I have not catched yet:
//  suppose there are three groups, it is compatible for any two of them. But consider all three groups, they are not compatible.
//  e.g.:   |    |      Three verticle bars like this.
//          |  | |      If watchDistance is large enough, this could happen.
//             |
////////////////////////////////////////////////////////////////////
__device__ bool VX3_VoxelGroup::isCompatible(VX3_Voxel *voxel_host, VX3_Voxel *voxel_remote, int *ret_linkdir_1, int *ret_linkdir_2) {
    if (voxel_host->d_group != this) {
        printf("This method should be call from voxel_host's d_group.\n"); // here?
        return false;
    }
    // Given two voxels, determine the best way to attach them.
    VX3_Vec3D<int> offset_of_link = VX3_Vec3D<int>(0, 0, 0);
    int potential_link_1, potential_link_2;
    VX3_Quat3D<double> relativeRotation = voxel_remote->orientation().Conjugate() * voxel_host->orientation();
    bool hasSingleton = false;

    if (true) { // Rotate singleton and small bar to align for attachment
        bool voxel_remote_singleton, voxel_remote_smallbar;
        int voxel_remote_direction;
        voxel_remote->isSingletonOrSmallBar(&voxel_remote_singleton, &voxel_remote_smallbar, &voxel_remote_direction);
        bool voxel_host_singleton, voxel_host_smallbar;
        int voxel_host_direction;
        voxel_host->isSingletonOrSmallBar(&voxel_host_singleton, &voxel_host_smallbar, &voxel_host_direction);
        if (voxel_remote_singleton || voxel_remote_smallbar || voxel_host_singleton || voxel_host_smallbar) {
            if (atomicCAS(&d_kernel->mutexRotateSingleton, 0, 1) == 0) {
                if (voxel_remote_singleton) { // remote has no link. so rotate the remote
                    voxel_remote->changeOrientationTo(voxel_host->orient);
                } else if (voxel_host_singleton) { // voxel host has no link, so rotate the host
                    voxel_host->changeOrientationTo(voxel_remote->orient);
                } else if (voxel_remote_smallbar) { // remote is a small bar  // they all have links, detach the one with less link
                    voxel_remote->adjacentVoxel((linkDirection)voxel_remote_direction)->changeOrientationTo(voxel_host->orient);
                    voxel_remote->changeOrientationTo(voxel_host->orient);
                    voxel_remote->links[voxel_remote_direction]->detach();
                } else if (voxel_host_smallbar) { // host is a small bar
                    voxel_host->adjacentVoxel((linkDirection)voxel_host_direction)->changeOrientationTo(voxel_remote->orient);
                    voxel_host->changeOrientationTo(voxel_remote->orient);
                    voxel_host->links[voxel_host_direction]->detach();
                }
                atomicExch(&d_kernel->mutexRotateSingleton, 0);
                hasSingleton = true;
            }
        }
    }

    if (relativeRotation.w > 0.866 || hasSingleton) // within 30 degree
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

    if (hasSingleton) { // no need to use group map to check for singletons
        if (voxel_host->links[potential_link_1] == NULL and voxel_remote->links[potential_link_2] == NULL) {
            *ret_linkdir_1 = potential_link_1;
            *ret_linkdir_2 = potential_link_2;
            return true;
        } else {
            // DEBUG_PRINT("%f) BAD position for the rotated singleton. The singleton might in between two linked voxels. Skip.\n", d_kernel->currentTime);
            return false;
        }
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
    bool ret = true;
    VX3_VoxelGroup *remote_group = voxel_remote->d_group;
    VX3_Vec3D<int> remote_diff = VX3_Vec3D<int>(0, 0, 0);
    remote_diff = voxel_host->groupPosition + offset_of_link - voxel_remote->groupPosition;
    for (int i = 0; i < remote_group->d_surface_voxels.size(); i++) {
        VX3_Voxel *v = remote_group->d_surface_voxels[i];
        VX3_Vec3D<int> position_of_remote_voxel_in_host_group = v->groupPosition + remote_diff;
        int offset = to1D(position_of_remote_voxel_in_host_group, VX3_Vec3D<int>(dim_x, dim_y, dim_z));
        if (offset == -1) {
            // good, because out of range
        } else if (d_group_map[offset] == NULL) {
            // good, because empty position
        } else {
            // printf("Not Compatible. Offset %d\n", offset); // Instead of return false, absorb the voxel!
            // TODO: Sida: Is this still in-place modifications?
            VX3_Voxel *voxel_to_absorb = d_group_map[offset];
            if (voxel_to_absorb == voxel_host) {
                ret = false; // Sida: this is a weird situation, the collision happened, but before this check, another voxel has been attached to this exact position. so this collision should not cause attachment.
            } else {
                voxel_to_absorb->removed = true;
                // delete all the links as well
                for (int i = 0; i < 6; i++) {
                    if (voxel_to_absorb->links[i]) {
                        voxel_to_absorb->links[i]->removed = true;
                        voxel_to_absorb->adjacentVoxel((linkDirection)i)->links[oppositeDirection(i)] = NULL;
                        voxel_to_absorb->links[i] = NULL;
                    }
                }
            }
        }
    }
    *ret_linkdir_1 = potential_link_1;
    *ret_linkdir_2 = potential_link_2;
    return ret;
}

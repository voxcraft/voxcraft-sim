#include "VX3_VoxelGroup.h"

#include "VX3_Voxel.h"

__device__ VX3_VoxelGroup::VX3_VoxelGroup(
    /* args */) {}

__device__ void VX3_VoxelGroup::updateGroup(VX3_Voxel *voxel) {
    if (needRebuild) {
        // Rebuild a map of Group
        //// BF Search for all neighbors
        VX3_dDictionary<VX3_Voxel *, int> BFS_visited;
        VX3_dQueue<VX3_Voxel *> BFS_Queue;
        VX3_dVector<VX3_Voxel *> BFS_result;

        printf("Start from (%d, %d, %d).\n", voxel->ix, voxel->iy, voxel->iz);
        BFS_result.push_back(voxel);
        BFS_visited.set(voxel, 1);
        BFS_Queue.push_back(voxel);

        while (!BFS_Queue.isEmpty()) {
            VX3_Voxel *v = BFS_Queue.pop_front();
            for (int i = 0; i < 6; i++) {
                VX3_Voxel *neighbor = v->adjacentVoxel((linkDirection)i);
                if (neighbor) {
                    if (BFS_visited.get(neighbor) == -1) {
                        printf("Visit (%d, %d, %d).\n", neighbor->ix, neighbor->iy, neighbor->iz);
                        // Set all connected voxels' d_group to this
                        neighbor->d_group = this;

                        BFS_result.push_back(neighbor);
                        BFS_visited.set(neighbor, 1);
                        BFS_Queue.push_back(neighbor);
                    }
                }
            }
        }

        //// Allocate memory for the map
        min_x = INT_MAX;
        min_y = INT_MAX;
        min_z = INT_MAX;
        max_x = 0;
        max_y = 0;
        max_z = 0;
        for (int i = 0; i < BFS_result.size(); i++) {
            VX3_Voxel *v = BFS_result[i];
            if (v->ix < min_x)
                min_x = v->ix;
            if (v->iy < min_y)
                min_y = v->iy;
            if (v->iz < min_z)
                min_z = v->iz;

            if (v->ix > max_x)
                max_x = v->ix;
            if (v->iy > max_y)
                max_y = v->iy;
            if (v->iz > max_z)
                max_z = v->iz;
        }
        printf("(%d,%d,%d) - (%d,%d,%d)\n", min_x, min_y, min_z, max_x, max_y, max_z);
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
        for (int i = 0; i < BFS_result.size(); i++) {
            VX3_Voxel *v = BFS_result[i];

            int offset = getVoxelOffset(v);
            d_group_map[offset] = v;
            d_voxels.push_back(v);
            for (int j = 0; j < 6; j++) {
                if (v->links[j] == NULL) {
                    d_surface_voxels.push_back(v);
                    break;
                }
            }
        }

        needRebuild = false;
    }
}

__device__ int VX3_VoxelGroup::getVoxelOffset(VX3_Voxel *voxel) {
    // calculate the offset for a voxel
    return to1D(voxel->ix - min_x, voxel->iy - min_y, voxel->iz - min_z);
}
__device__ int VX3_VoxelGroup::to1D(int x, int y, int z) {
    if (x < 0 || y < 0 || z < 0 || x >= dim_x || y >= dim_y || z >= dim_z)
        return -1; // Overflow
    int offset;
    offset = dim_y * dim_z * x + dim_z * y + z;
    return offset;
}
__device__ void VX3_VoxelGroup::to3D(int offset, int *ret_x, int *ret_y, int *ret_z) {
    int x, y, z;
    z = offset % dim_z;
    int residual;
    residual = (int)((offset - z) / dim_z);
    y = residual % dim_y;
    residual = (int)((residual - y) / dim_y);
    x = residual;
    *ret_x = x;
    *ret_y = y;
    *ret_z = z;
}

__device__ bool VX3_VoxelGroup::isCompatible(VX3_Voxel *voxel_host, VX3_Voxel *voxel_remote, int *ret_linkdir_1, int *ret_linkdir_2) {
    if (voxel_host->d_group != this) {
        printf("This method should be call from voxel_host's d_group.\n");
        return false;
    }
    // Determine Angle
    VX3_Quat3D<double> q1 = voxel_host->orientation();
    VX3_Quat3D<double> q2 = voxel_remote->orientation();
    // printf("q1 w=%f,x=%f,y=%f,z=%f, ", q1.w, q1.x, q1.y, q1.z);
    // printf("q2 w=%f,x=%f,y=%f,z=%f\n", q2.w, q2.x, q2.y, q2.z);
    VX3_Quat3D<double> q2to1 = q2.Conjugate() * q1;
    // printf("q2to1 w=%f,x=%f,y=%f,z=%f\n", q2to1.w, q2to1.x, q2to1.y, q2to1.z);

    int potential_link_1 = -1, potential_link_2 = -1;
    // Quaternion = cos(degree) + sin(degree)(a i+b j+c k)
    // when degree = -135, -90, -45, 0, 45, 90, 135, 180, it is orthogonal
    // If we allow +- 8 degree from orthogonal, that is:
    // check if w is in (cos(-135-8) ~ cos(-135+8)) or (cos(-90-8) ~ cos(-90+8)) or (cos(-45-8) ~ cos(-45+8)) , etc.
    float w = abs(q2to1.w);
    bool almost_orthogonal = false;
    if ((w < 0.139173) || (w > 0.601815 && w < 0.798636) || (w > 0.990268)) {
        almost_orthogonal = true;
    }
    if (almost_orthogonal) {
        // Rotate voxel_remote in lattice to align with voxel_host
        if (true) {
            // printf("this special case.");
            voxel_remote->d_group->reorient_lattice();
        }

        // Almost in the same direction
        VX3_Vec3D<double> raw_pos = voxel_remote->position() - voxel_host->position();
        VX3_Vec3D<double> pos_1 = q1.RotateVec3DInv(raw_pos);  // make the pos relative to voxel_host
        VX3_Vec3D<double> pos_2 = q2.RotateVec3DInv(-raw_pos); // make the pos relative to voxel_remote

        // printf("relative pos_1: %f, %f, %f\n\n", pos_1.x, pos_1.y, pos_1.z);
        if (abs(pos_1.x) > abs(pos_1.y) && abs(pos_1.x) > abs(pos_1.z)) {
            // X-axis
            potential_link_1 = (int)(pos_1.x > 0 ? X_POS : X_NEG);
        } else if (abs(pos_1.y) > abs(pos_1.z)) {
            // Y-axis
            potential_link_1 = (int)(pos_1.y > 0 ? Y_POS : Y_NEG);
        } else {
            // Z-axis
            potential_link_1 = (int)(pos_1.z > 0 ? Z_POS : Z_NEG);
        }

        // printf("relative pos_2: %f, %f, %f\n\n", pos_2.x, pos_2.y, pos_2.z);
        if (abs(pos_2.x) > abs(pos_2.y) && abs(pos_2.x) > abs(pos_2.z)) {
            // X-axis
            potential_link_2 = (int)(pos_2.x > 0 ? X_POS : X_NEG);
        } else if (abs(pos_2.y) > abs(pos_2.z)) {
            // Y-axis
            potential_link_2 = (int)(pos_2.y > 0 ? Y_POS : Y_NEG);
        } else {
            // Z-axis
            potential_link_2 = (int)(pos_2.z > 0 ? Z_POS : Z_NEG);
        }
    }

    if (potential_link_1 == -1 || potential_link_2 == -1) {
        // Too large of an angle
        return false;
    }
    int remote_ix, remote_iy, remote_iz;
    // TODO: should rotate remote according to q2to1
    // Assume direction is X_NEG for now
    remote_ix = -voxel_remote->iy;
    remote_iy = voxel_remote->ix;
    remote_iz = voxel_remote->iz;

    int diff_ix = voxel_host->ix - remote_ix - min_x;
    int diff_iy = voxel_host->iy - remote_iy - min_y;
    int diff_iz = voxel_host->iz - remote_iz - min_z;
    switch ((linkDirection)potential_link_1) {
    case X_POS:
        diff_ix++;
        break;
    case X_NEG:
        diff_ix--;
        break;
    case Y_POS:
        diff_iy++;
        break;
    case Y_NEG:
        diff_iy--;
        break;
    case Z_POS:
        diff_iz++;
        break;
    case Z_NEG:
        diff_iz--;
        break;
    default:
        printf("TODO: no valid link direction.\n");
    }

    // Only need to consider surface voxel
    // for (int i = 0; i < d_surface_voxels.size(); i++) {
    //     printf("Surface 1: %d, %d, %d. \n", d_surface_voxels[i]->ix - min_x, d_surface_voxels[i]->iy - min_y, d_surface_voxels[i]->iz - min_z);
    // }

    for (int i = 0; i < voxel_remote->d_group->d_surface_voxels.size(); i++) {
        remote_ix = -voxel_remote->d_group->d_surface_voxels[i]->iy;
        remote_iy = voxel_remote->d_group->d_surface_voxels[i]->ix;
        remote_iz = voxel_remote->d_group->d_surface_voxels[i]->iz;
        int remote_x = remote_ix + diff_ix;
        int remote_y = remote_iy + diff_iy;
        int remote_z = remote_iz + diff_iz;
        // printf("Surface 2: %d, %d, %d. \n", remote_x, remote_y, remote_z);

        int offset = to1D(remote_x, remote_y, remote_z);
        if (offset == -1) {
            // good, because out of range
        } else if (d_group_map[offset] == NULL) {
            // good, because empty position
        } else {
            // printf("Not Compatible. Offset %d\n", offset);
            // return false;
        }
    }
    // printf("Compatible.\n");
    *ret_linkdir_1 = potential_link_1;
    *ret_linkdir_2 = potential_link_2;
    return true;
}

__device__ void VX3_VoxelGroup::reorient_lattice() {
    // int c = 0;
    // VX3_Link* tmp[6];
    // VX3_dDictionary<VX3_Link *, int> link_visited;
    // for (int i=0;i<d_voxels.size();i++) {
    //     // 1. change link->axis
    //     for (int j=0;j<6;j+=2) {
    //         VX3_Link* link = d_voxels[i]->links[j];
    //         if (link) {
    //             if (link_visited.get(link)==-1) {
    //                 link_visited.set(link, 1);
    //                 printf("%d) link %p.\n", c++, link);
    //                 // if (link->axis == X_AXIS) {
    //                 //     link->axis = Y_AXIS;
    //                 //     // VX3_Voxel *tmp;
    //                 //     // tmp = link->pVNeg;
    //                 //     // link->pVNeg = link->pVPos;
    //                 //     // link->pVPos = tmp;
    //                 // } else {
    //                 //     link->axis = X_AXIS;
    //                 // }
    //                 link->reset();
    //                 link->isNewLink = 10;
    //             }
    //         }
    //     }
    //     // 2. change links
    //     for (int j=0;j<6;j++) {
    //         tmp[j] = d_voxels[i]->links[j];
    //     }
    //     d_voxels[i]->links[(int)X_POS] = tmp[(int)Y_NEG];
    //     d_voxels[i]->links[(int)X_NEG] = tmp[(int)Y_POS];
    //     d_voxels[i]->links[(int)Y_POS] = tmp[(int)X_POS];
    //     d_voxels[i]->links[(int)Y_NEG] = tmp[(int)X_NEG];
    //     // 3. change orientation
    //     // d_voxels[i]->orient = d_voxels[i]->orient * q2to1
    //     double num71 = cos(3.1415926/4.0);
    //     d_voxels[i]->orient = VX3_Quat3D<double>(num71, 0, 0, num71);
    // }
}
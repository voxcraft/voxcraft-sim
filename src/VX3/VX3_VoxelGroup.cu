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

__device__ bool VX3_VoxelGroup::isCompatible(VX3_Voxel *voxel_host, VX3_Voxel *voxel_remote, int *ret_linkdir) {
    if (voxel_host->d_group != this) {
        printf("This method should be call from voxel_host's d_group.\n");
        return false;
    }
    // Determine Angle
    VX3_Quat3D<double> q1 = voxel_host->orientation();
    VX3_Quat3D<double> q2 = voxel_remote->orientation();
    printf("q1 w=%f,x=%f,y=%f,z=%f, ", q1.w, q1.x, q1.y, q1.z);
    printf("q2 w=%f,x=%f,y=%f,z=%f\n", q2.w, q2.x, q2.y, q2.z);
    VX3_Quat3D<double> q2to1 = q2.Conjugate() * q1;
    printf("q2to1 w=%f,x=%f,y=%f,z=%f\n", q2to1.w, q2to1.x, q2to1.y, q2to1.z);

    int potential_link = -1;
    // TODO: need consider all Almost orthogonal situations! refer to Quat
    if (q2to1.w > 0.99f) {
        // Almost in the same direction
        VX3_Vec3D<double> pos = voxel_remote->position() - voxel_host->position();
        pos = q1.RotateVec3DInv(pos); // make the pos relative to voxel_host

        printf("relative pos: %f, %f, %f\n\n", pos.x, pos.y, pos.z);
        if (abs(pos.x) > abs(pos.y) && abs(pos.x) > abs(pos.z)) {
            // X-axis
            potential_link = (int)pos.x > 0 ? X_POS : X_NEG;
        } else if (abs(pos.y) > abs(pos.z)) {
            // Y-axis
            potential_link = (int)pos.y > 0 ? Y_POS : Y_NEG;
        } else {
            // Z-axis
            potential_link = (int)pos.z > 0 ? Z_POS : Z_NEG;
        }
    }

    if (potential_link == -1) {
        // Too large of an angle
        return false;
    }
    // Assume direction is X_NEG for now
    int diff_ix = voxel_host->ix - voxel_remote->ix - min_x;
    int diff_iy = voxel_host->iy - voxel_remote->iy - min_y;
    int diff_iz = voxel_host->iz - voxel_remote->iz - min_z;
    switch ((linkDirection)potential_link) {
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
    for (int i = 0; i < d_surface_voxels.size(); i++) {
        printf("Surface 1: %d, %d, %d. \n", d_surface_voxels[i]->ix - min_x, d_surface_voxels[i]->iy - min_y, d_surface_voxels[i]->iz - min_z);
    }

    for (int i = 0; i < voxel_remote->d_group->d_surface_voxels.size(); i++) {
        int remote_x = voxel_remote->d_group->d_surface_voxels[i]->ix + diff_ix;
        int remote_y = voxel_remote->d_group->d_surface_voxels[i]->iy + diff_iy;
        int remote_z = voxel_remote->d_group->d_surface_voxels[i]->iz + diff_iz;
        printf("Surface 2: %d, %d, %d. \n", remote_x, remote_y, remote_z);

        int offset = to1D(remote_x, remote_y, remote_z);
        if (offset == -1) {
            // good, because out of range
        } else if (d_group_map[offset] == NULL) {
            // good, because empty position
        } else {
            printf("Not Compatible. Offset %d\n", offset);
            return false;
        }
    }
    printf("Compatible.\n");
    *ret_linkdir = potential_link;
    return true;
}

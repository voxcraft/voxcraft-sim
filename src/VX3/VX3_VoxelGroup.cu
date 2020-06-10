#include "VX3_VoxelGroup.h"

#include "VX3_Voxel.h"
#include "VX3_dictionary.cuh"
#include "VX3_queue.cuh"

__device__ VX3_VoxelGroup::VX3_VoxelGroup(/* args */) {}

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
                        printf("Visit (%d, %d, %d).\n", neighbor->ix,
                               neighbor->iy, neighbor->iz);
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
        int min_x = INT_MAX, min_y = INT_MAX, min_z = INT_MAX, max_x = 0,
            max_y = 0, max_z = 0;
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
        printf("(%d,%d,%d) - (%d,%d,%d)\n", min_x,min_y,min_z,max_x,max_y,max_z);
        if (d_group_map) {
            free(d_group_map);
            d_group_map = NULL;
        }

        needRebuild = false;
    }
    // Update Position
}
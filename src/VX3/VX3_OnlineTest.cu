#include "VX3_OnlineTest.h"
#include "VX3_VoxelyzeKernel.cuh"

// check all the voxels, links and voxelgroups for validation.
__device__ bool VX3_OnlineTest::ThoroughTest(VX3_VoxelyzeKernel *k) {
    if (k->CurStepCount < k->ThoroughTestStartAt || k->CurStepCount % k->ThoroughTestStepSize > 0) {
        // skip most of the steps.
        return true;
    }
    printf("GPU (%d) Online Testing at step %lu ...\n", k->GPU_id, k->CurStepCount);
    bool ret = true;

    // Testing all voxels.
    for (int i = 0; i < k->num_d_voxels; i++) {
        VX3_Voxel *v = &k->d_voxels[i];
        if (v->removed)
            continue;
        if (v->d_kernel != k) {
            printf("ERROR: voxel d_kernel. A voxel should keep a ptr of kernel in it.\n");
            ret = false;
        }
        if (v->d_group == NULL) {
            printf("ERROR: voxel d_group. A voxel cannot belong to a NULL group.\n");
            ret = false;
        }
        if (v->d_group->removed) {
            printf("ERROR: voxel d_group. A voxel (%p) cannot belong to a removed group: %p.\n", v, v->d_group);
            ret = false;
        }
        
        bool hit = false;
        for (int j = 0; j < v->d_group->d_voxels.size(); j++) {
            if (v->d_group->d_voxels[j] == v) {
                hit = true;
            }
        }
        if (hit == false) {
            printf("ERROR: voxel d_group voxel. A voxel (%p) should not belong to a group (%p) that doesn't contain itself.\n", v, v->d_group);
            ret = false;
        }

        for (int j = 0; j < 6; j++) {
            if (v->links[j]) {
                if (v->links[j]->removed) {
                    printf("ERROR: links->removed. Set the link ptr to NULL if the link was destoried.\n");
                    ret = false;
                }
                if (v->adjacentVoxel((linkDirection)j)->removed) {
                    printf("ERROR: adjacent voxel removed. A voxel (%p) cannot link (%p) to a removed voxel (%p).\n", v, v->links[j], v->adjacentVoxel((linkDirection)j));
                    ret = false;
                }
                if (v->adjacentVoxel((linkDirection)j)->adjacentVoxel((linkDirection)oppositeDirection(j)) != v) {
                    printf("ERROR: adjacent voxel inconsistent. A voxel's neighbor should point back to itself.\n");
                    ret = false;
                }
            }
        }
    }

    // Testing all voxelgroups
    VX3_dDictionary<int, VX3_Voxel *> dic;
    for (int i = 0; i < k->d_voxelgroups.size(); i++) {
        VX3_VoxelGroup *g = k->d_voxelgroups[i];
        if (g->removed)
            continue;
        VX3_Vec3D<int> dim = VX3_Vec3D<int>(g->dim_x, g->dim_y, g->dim_z);
        dic.clear();
        for (int j = 0; j < g->d_voxels.size(); j++) {
            VX3_Voxel *v = g->d_voxels[j];
            if (v->removed) {
                printf("ERROR: group voxel removed. A group (%p) should not contain a removed voxel (%p).\n", g, v);
                ret = false;
            }

            if (v->d_group != g) {
                printf("ERROR: group voxel inconsistent. A group (%p) should not contain a voxel (%p) that is not belong to this group (%p).\n", g, v, v->d_group);
                ret = false;
            }

            if (dic.get(to1D(v->groupPosition, dim)) != (VX3_Voxel *)-1) {
                printf("ERROR: group position is the same. In one voxel group (%p), there cannot be two voxels (%p, %p) with the same group position.\n", v->d_group, v, dic.get(to1D(v->groupPosition, dim)));
                ret = false;
            }
            dic.set(to1D(v->groupPosition, dim), v);
        }
    }

    // Testing all links
    for (int i = 0; i < k->num_d_links; i++) {
        VX3_Link *l = &k->d_links[i];
        if (l->removed)
            continue;

        if (l->pVPos->removed) {
            printf("ERROR: link pVPos removed. A link should not attach to a removed voxel.\n");
            ret = false;
        }
        if (l->pVNeg->removed) {
            printf("ERROR: link pVNeg removed. A link should not attach to a removed voxel.\n");
            ret = false;
        }
        if (l->pVPos->links[l->linkdirPos] != l) {
            printf("ERROR: link inconsistent. A link's positive end should point to itself.\n");
            ret = false;
        }
        if (l->pVNeg->links[l->linkdirNeg] != l) {
            printf("ERROR: link inconsistent. A link's positive end should point to itself.\n");
            ret = false;
        }
    }
    return ret;
}
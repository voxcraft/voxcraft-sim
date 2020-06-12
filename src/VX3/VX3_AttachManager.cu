#include "VX3_AttachManager.h"

#include "VX3_Voxel.h"
#include "VX3_VoxelyzeKernel.cuh"

__device__ VX3_AttachManager::VX3_AttachManager(/* args */) {}

__device__ int VX3_AttachManager::oppositeDir(int linkdir) {
    // X_NEG for X_POS, etc.
    int residual = linkdir % 2;
    return linkdir - residual + !(linkdir % 2);
}
__device__ bool VX3_AttachManager::tryAttach(VX3_Voxel *voxel1, VX3_Voxel *voxel2) {
    printf("try attach: voxel1 (%d,%d,%d), voxel2 (%d,%d,%d)\n", voxel1->ix, voxel1->iy, voxel1->iz, voxel2->ix, voxel2->iy, voxel2->iz);

    // Simple preliminary filter
    //// determined by formula
    if (!voxel1->enableAttach || !voxel2->enableAttach)
        return false;
    //// fixed voxels, no need to look further for attachment
    if (voxel1->mat->fixed || voxel2->mat->fixed)
        return false;
    //// different material, no need to attach
    if (voxel1->mat != voxel2->mat)
        return false;
    //// sticky or not
    if (!voxel1->mat->sticky)
        return false;

    // Check VoxelGroup map for compatible
    if (false) {
        //// consider rotation later
        VX3_Quat3D<float> q1 = voxel1->orientation();
        VX3_Quat3D<float> q2 = voxel2->orientation();
        // printf("q1 w=%f,x=%f,y=%f,z=%f, ", q1.w, q1.x, q1.y, q1.z);
        // printf("q2 w=%f,x=%f,y=%f,z=%f\n", q2.w, q2.x, q2.y, q2.z);
        VX3_Quat3D<float> q12 = q1.Conjugate() * q2;
        // printf("q12 w=%f,x=%f,y=%f,z=%f\n", q12.w, q12.x, q12.y, q12.z);
    }
    int linkdir_1, linkdir_2;
    if (voxel1->d_group->isCompatible(voxel1, voxel2, &linkdir_1, &linkdir_2)) {
        // linkdir_2 = oppositeDir(linkdir_1);
        // printf("Compatible.");
    } else {
        return false;
    }

    // try only form one link
    if (debug)
        return false;

    bool ret = false;
    // Start Attach!
    // Other potential attachments are ignored, no wait.
    if (atomicCAS(&mutex, 0, 1) == 0) {
        // Entering Critical Area
        if (voxel1->links[linkdir_1] == NULL && voxel2->links[linkdir_2] == NULL) {
            VX3_Link *pL;
            pL = new VX3_Link(voxel2, (linkDirection)linkdir_2, voxel1, (linkDirection)linkdir_1, k); // make the new link (change to both materials, etc.
            if (!pL) {
                printf(COLORCODE_BOLD_RED "ERROR: Out of memory. Link not created.\n");
            } else {
                pL->isNewLink = k->SafetyGuard;
                k->d_v_links.push_back(pL); // add to the list

                k->isSurfaceChanged = true;
                printf("New Link formed.\n");
                k->EnableCilia = false; // for debug: no cilia after attachment.
                ret = true;
                debug++;
            }
        }
        atomicExch(&mutex, 0);
    }

    return ret;
}

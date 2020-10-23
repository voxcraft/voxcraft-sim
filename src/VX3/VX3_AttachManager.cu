#include "VX3_AttachManager.h"

#include "VX3_Voxel.h"
#include "VX3_VoxelyzeKernel.cuh"

__device__ VX3_AttachManager::VX3_AttachManager(VX3_VoxelyzeKernel *k) { d_kernel = k; }

__device__ bool VX3_AttachManager::attachWhileCollide(VX3_Voxel *voxel1, VX3_Voxel *voxel2) {
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

    // sam:
    if (voxel1->nonStickTimer > d_kernel->currentTime || voxel2->nonStickTimer > d_kernel->currentTime)
        return false;

    // Check VoxelGroup map for compatible
    int linkdir_1, linkdir_2;

    if (voxel1->d_group == voxel2->d_group) {
        // Same group, but please attach neighbors in the same group.
        VX3_Vec3D<int> diff = voxel1->groupPosition - voxel2->groupPosition;
        if (diff.Length2() > 1)
            return false;
        linkdir_1 = (diff.x == 1 ? 1 : 0) + (diff.y == 1 ? 3 : 0) + (diff.y == -1 ? 2 : 0) + (diff.z == 1 ? 5 : 0) + (diff.z == -1 ? 4 : 0);
        linkdir_2 = oppositeDirection(linkdir_1);

    } else {
        // Different groups, check for compatibility
        if (!voxel1->d_group->isCompatible(voxel1, voxel2, &linkdir_1, &linkdir_2))
            return false;
    }

    // try only form one link
    if (OnlyFormOneLink) {
        if (totalLinksFormed >= 1)
            return false;
    }
    return tryToAttach(voxel1, linkdir_1, voxel2, linkdir_2);
}


__device__ bool VX3_AttachManager::attachForNewVoxel(VX3_Voxel *voxel1, int linkdir_1, VX3_Voxel *voxel2, int linkdir_2) {
    while(1) {
        if (tryToAttach(voxel1, linkdir_1, voxel2, linkdir_2)) {
            break;
        }
    }
    return true;
}

__device__ bool VX3_AttachManager::tryToAttach(VX3_Voxel *voxel1, int linkdir_1, VX3_Voxel *voxel2, int linkdir_2) {

    bool ret = false;
    // Start Attach!
    // Only once attchment at a time, other potential attachments should be ignored, no wait.
    if (atomicCAS(&attachmentMutex, 0, 1) == 0) {
        // Entering Critical Area
        if (!voxel1->d_group->needUpdate && !voxel2->d_group->needUpdate) { // to avoid two voxels attach at the same position
            if (voxel1->links[linkdir_1] == NULL && voxel2->links[linkdir_2] == NULL) {
                VX3_Link *pL;
                pL = new VX3_Link(voxel1, (linkDirection)linkdir_1, voxel2, (linkDirection)linkdir_2, d_kernel); // make the new link (change to both materials, etc.
                if (!pL) {
                    printf(COLORCODE_BOLD_RED "ERROR: Out of memory. Link not created.\n");
                } else {
                    // update voxel1's group, it will set all connected voxels' group to voxel1's
                    voxel1->d_group->hasNewLink++;
                    voxel1->d_group->needUpdate = 1;
                    d_kernel->d_voxel_to_update_group.push_back(voxel1);

                    // if (voxel1->d_group != voxel2->d_group) {
                    //     voxel2->d_group->switchAllVoxelsTo(voxel1->d_group);
                    //     // voxel1->d_group->updateGroup(voxel1);
                    // }
                    pL->isNewLink = d_kernel->SafetyGuard;
                    d_kernel->d_v_links.push_back(pL); // add to the list

                    d_kernel->isSurfaceChanged = true;
                    // DEBUG_PRINT("%f) New Link formed.\n", d_kernel->currentTime);
                    // d_kernel->EnableCilia = false; // for debug: no cilia after attachment.
                    ret = true;
                    totalLinksFormed++;
                    
                    // sam:
                    voxel1->targetPos.clear();
                    voxel2->targetPos.clear();

                }
            }
        }
        atomicExch(&attachmentMutex, 0);
    }

    return ret;
}
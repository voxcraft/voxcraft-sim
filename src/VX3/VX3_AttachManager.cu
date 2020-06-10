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
        printf("q1 w=%f,x=%f,y=%f,z=%f, ", q1.w, q1.x, q1.y, q1.z);
        printf("q2 w=%f,x=%f,y=%f,z=%f\n", q2.w, q2.x, q2.y, q2.z);
        VX3_Quat3D<float> q12 = q1.Conjugate() * q2;
        printf("q12 w=%f,x=%f,y=%f,z=%f\n", q12.w, q12.x, q12.y, q12.z);
    }
    int linkdir_1, linkdir_2;
    if (voxel1->d_group->isCompatible(voxel1, voxel2, &linkdir_1)) {
        // printf("Compatible.");
        linkdir_2 = oppositeDir(linkdir_1);
    } else {
        return false;
    }

    bool ret = false;
    // Start Attach!
    // Other potential attachments are ignored, no wait.
    if (atomicCAS(&mutex, 0, 1) == 0) {
        // Entering Critical Area
        // determine relative position
        linkDirection link_dir_1, link_dir_2;
        linkAxis link_axis;
        auto a = voxel1->orientation();
        auto b = voxel2->orientation();
        auto c = voxel1->position();
        auto d = voxel2->position();
        auto e = c - d;
        auto ea = a.RotateVec3DInv(-e);
        auto eb = b.RotateVec3DInv(e);

        // first find which is the dominant axis, then determine which one is
        // neg which one is pos.
        VX3_Vec3D<double> f;
        bool reverseOrder = false;
        f = ea.Abs();
        if (f.x >= f.y && f.x >= f.z) { // X_AXIS
            link_axis = X_AXIS;
            if (ea.x < 0) {
                link_dir_1 = X_NEG;
                link_dir_2 = X_POS;
                reverseOrder = true;
            } else {
                link_dir_1 = X_POS;
                link_dir_2 = X_NEG;
            }
        } else if (f.y >= f.x && f.y >= f.z) { // Y_AXIS
            link_axis = Y_AXIS;
            if (ea.y < 0) {
                link_dir_1 = Y_NEG;
                link_dir_2 = Y_POS;
                reverseOrder = true;
            } else {
                link_dir_1 = Y_POS;
                link_dir_2 = Y_NEG;
            }
        } else { // Z_AXIS
            link_axis = Z_AXIS;
            if (ea.z < 0) { // voxel1 is on top
                link_dir_1 = Z_NEG;
                link_dir_2 = Z_POS;
                reverseOrder = true;
            } else {
                link_dir_1 = Z_POS;
                link_dir_2 = Z_NEG;
            }
        }

        // TODO: need to solve this. Create only when there's a right place to
        // attach
        printf("link_dir_1 %d %d, link_dir_2 %d %d\n", (int)link_dir_1, linkdir_1, (int)link_dir_2, linkdir_2);
        if (voxel1->links[link_dir_1] == NULL && voxel2->links[link_dir_2] == NULL) {
            VX3_Link *pL;
            if (reverseOrder) {
                pL = new VX3_Link(voxel1, link_dir_1, voxel2, link_dir_2, link_axis,
                                  k); // make the new link (change to both materials, etc.
            } else {
                pL = new VX3_Link(voxel2, link_dir_2, voxel1, link_dir_1, link_axis,
                                  k); // make the new link (change to both materials, etc.
            }
            if (!pL) {
                printf(COLORCODE_BOLD_RED "ERROR: Out of memory. Link not created.\n");
            } else {
                pL->isNewLink = k->SafetyGuard;
                k->d_v_links.push_back(pL); // add to the list

                k->isSurfaceChanged = true;
                printf("New Link formed.\n");
                ret = true;
            }
        }
        atomicExch(&mutex, 0);
    }

    return ret;
}

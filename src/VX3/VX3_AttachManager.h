//
// Created by Sida Liu
//  This class manages the attachment.
//
#if !defined(VX3_ATTACHMANAGER_H)
#define VX3_ATTACHMANAGER_H

class VX3_Voxel;
class VX3_VoxelyzeKernel;

class VX3_AttachManager {
public:
    /* data */
    VX3_VoxelyzeKernel *d_kernel; // saved pointer to the whole simulation
    int attachmentMutex = 0; // only one attachment can happen at a time

    /* method */
    __device__ VX3_AttachManager(VX3_VoxelyzeKernel *k);
    __device__ bool tryAttach(VX3_Voxel *voxel1, VX3_Voxel *voxel2);
    __device__ bool doAttach(VX3_Voxel *voxel1, int linkdir_1, VX3_Voxel *voxel2, int linkdir_2);


    int OnlyFormOneLink = false; // for debug use
    int totalLinksFormed = 0;
};

#endif // VX3_ATTACHMANAGER_H

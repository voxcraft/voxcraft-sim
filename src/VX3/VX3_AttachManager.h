#if !defined(VX3_ATTACHMANAGER_H)
#define VX3_ATTACHMANAGER_H

class VX3_Voxel;
class VX3_VoxelyzeKernel;

class VX3_AttachManager
{
public:
    /* data */
    VX3_VoxelyzeKernel* k;
    int mutex = 0;

    /* method */
    __device__ VX3_AttachManager(/* args */);
    __device__ bool tryAttach(VX3_Voxel* voxel1, VX3_Voxel* voxel2);
    __device__ int oppositeDir(int linkdir);


    int debug=0;
};



#endif // VX3_ATTACHMANAGER_H

#if !defined(VX3_ONLINE_TEST_H)
#define VX3_ONLINE_TEST_H

class VX3_VoxelyzeKernel;

class VX3_OnlineTest
{
private:
    /* data */
public:
    __device__ bool ThoroughTest(VX3_VoxelyzeKernel *k);
};

#endif // VX3_ONLINE_TEST_H

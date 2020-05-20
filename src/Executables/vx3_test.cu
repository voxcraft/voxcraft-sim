#include "VX3_VoxelyzeKernel.cuh"
#include "Voxelyze.h"

__global__ void kernel(VX3_VoxelyzeKernel* d_voxelyze_3) {
    int count = 0;
    for (int i=0;i<d_voxelyze_3->num_d_voxels;i++) {
        for (int j=0;j<i;j++) {
            auto pl = new VX3_Link(&d_voxelyze_3->d_voxels[i], Y_POS, &d_voxelyze_3->d_voxels[j], Y_NEG, Y_AXIS, d_voxelyze_3);
            printf("%d ,  pl %p.\n", count, pl);
            if (!pl) break;
            count ++;
        }
    }
}

void enlarge() {
    size_t HeapSize = 1;
    double ratio = 0.5; // make 10% of the total GPU memory to be heap memory
    size_t free, total;
    VcudaMemGetInfo(&free, &total);
    printf("Total GPU memory %ld bytes.\n", total);
    for (int i = 0; i < 100; i++) {
        if (HeapSize >= total * ratio)
            break;
        HeapSize *= 2;
    }
    HeapSize += 1024; // add some additional size
    printf("Set GPU heap size to be %ld bytes.\n", HeapSize);
    VcudaDeviceSetLimit(
        cudaLimitMallocHeapSize,
        HeapSize); // Set Heap Memory to 1G, instead of merely 8M.

}

int main() {
    // enlarge();
    
    CVoxelyze Vx(0.005); // 5mm voxels
    CVX_Material *pMaterial = Vx.addMaterial(
        1000000,
        1000); // A material with stiffness E=1MPa and density 1000Kg/m^3
    for (int i = 0; i < 1000; i++) {
        Vx.setVoxel(pMaterial, i * 2, 0, 0); // Voxel at index x=0, y=0. z=0
    }
    CVX_Environment MainEnv;
    CVX_Sim MainSim;
    CVX_Object MainObj;
    MainEnv.pObj = &MainObj; // connect environment to object
    MainSim.pEnv = &MainEnv; // connect Simulation to envirnment
    MainSim.Vx = Vx;

    VX3_VoxelyzeKernel h_d_tmp(&MainSim);
    VX3_VoxelyzeKernel *d_voxelyze_3;

    VcudaMalloc(&d_voxelyze_3, sizeof(VX3_VoxelyzeKernel));
    VcudaMemcpy(d_voxelyze_3, &h_d_tmp, sizeof(VX3_VoxelyzeKernel),
                cudaMemcpyHostToDevice);
    kernel<<<1, 1>>>(d_voxelyze_3);
    cudaDeviceSynchronize();
    return 0;

    // CVoxelyze Vx(0.005); // 5mm voxels
    // CVX_Material *pMaterial = Vx.addMaterial(
    //     1000000,
    //     1000); // A material with stiffness E=1MPa and density 1000Kg/m^3
    // CVX_Voxel *Voxel1 =
    //     Vx.setVoxel(pMaterial, 0, 0, 0); // Voxel at index x=0, y=0. z=0
    // CVX_Voxel *Voxel2 =
    //     Vx.setVoxel(pMaterial, 1, 0, 0); // Voxel at index x=0, y=0. z=0
    // Voxel2->pos.x -= 0.0001;
    // for (int i = 0; i < 100; i++)
    //     Vx.doTimeStep(); // simulates 100 timesteps
}
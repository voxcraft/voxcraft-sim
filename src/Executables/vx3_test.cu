#include "VX3_SignalDiffusion.h"
// GPU Heap is for in-kernel malloc(). Refer to
// https://stackoverflow.com/a/34795830/7001199
void enlargeGPUHeapSize() {
    size_t HeapSizeInBytes;
    size_t free, total;
    VcudaMemGetInfo(&free, &total);
    printf("Total GPU memory %ld bytes.\n", total);
    HeapSizeInBytes = 0.5 * total; // add some additional size
    printf("Set GPU heap size to be %ld bytes.\n", HeapSizeInBytes);
    VcudaDeviceSetLimit(cudaLimitMallocHeapSize,
                        HeapSizeInBytes); // Set Heap Memory to 1G, instead of merely 8M.

    // if "Lane User Stack Overflow" ocurs, maybe Stack Size too small, can try this:
    // VcudaDeviceSetLimit(cudaLimitStackSize, 2048);
}

__global__ void kernel() {
    VX3_SignalDiffusion d_solute;
    d_solute.deviceInit(VX3_Vec3D<int>(10,1,1));
    VX3_Vec3D<int> o = VX3_Vec3D<int>(3,0,0);
    d_solute.addSolute(100.258, o);
    for (int j=0;j<1000;j++) {
        d_solute.doTimestep(0.001);
        for (int i=0;i<10;i++) {
            printf("%f, ", *(d_solute.map+i));
        }
        printf("\n");
    }
    double v = d_solute.quaryQuantityAtPosition(o);
    printf("\nv: %f\n", v);
}

int main() {
    enlargeGPUHeapSize();
    kernel<<<1, 1>>>();
    cudaDeviceSynchronize();
}
#ifdef _0
#include <stdio.h>
#include <iostream>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <vector>

#include "TI_Object.h"
#include "VX3.cuh"


__global__ void kernel(float *A, unsigned num_A) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<num_A) {
        atomicAdd(&A[0],1.0f);
    }
}
__global__ void kernel2(Data *data, unsigned num) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<num) {
        data[0]._i ++;
        atomicAdd(&data[1]._i, 1);
    }
}
CTI_Object::CTI_Object() {
    
}

Data ** CTI_Object::Try4() {
    std::vector<Data *> h_vec_data;
    h_vec_data.push_back(new Data(3,4.0));
    h_vec_data.push_back(new Data(8,4.0));
    printf("length: %d\n", (int)h_vec_data.size());

    Data ** h_data;
    unsigned num = h_vec_data.size();
    h_data = (Data **)malloc(sizeof(Data *)*num);
    for (unsigned i=0;i<h_vec_data.size();i++) {
        h_data[i] = h_vec_data[i];
    }
    //std::copy(h_vec_data.begin(), h_vec_data.end(), &h_data);

    return h_data;

}

void CTI_Object::Try3() {
    Data * h_data;
    Data * d_data;
    int num = 10;
    int size = num*sizeof(Data);
    h_data = (Data *) malloc(size);
    for (unsigned i=0;i<num;i++) {
        h_data[i]._i = i;
        h_data[i]._d = i+1.5;
    }
    VcudaMalloc(&d_data, size);
    VcudaMemcpy(d_data, h_data, size, VcudaMemcpyHostToDevice);
    kernel2<<<1,1024>>>(d_data, num);
    VcudaMemcpy(h_data, d_data, size, VcudaMemcpyDeviceToHost);

    for (unsigned i=0;i<num;i++) {
        printf("%d ", h_data[i]._i);
    }
}

void CTI_Object::Try() {
    float* d_A;
    float* h_A;
    unsigned num_A = 1024;
    unsigned size = num_A*num_A * sizeof(float);
    VcudaMalloc(&d_A, size);
    kernel<<<1024,1024>>>(d_A, num_A*num_A);
    h_A = (float *)malloc( size );
    VcudaMemcpy(h_A, d_A, sizeof(float), VcudaMemcpyDeviceToHost );
    VcudaFree(d_A);
    printf("%f, ", h_A[0]);
    printf("\n");
    delete h_A;
}

void CTI_Object::Try2() {
    int major = THRUST_MAJOR_VERSION;
    int minor = THRUST_MINOR_VERSION;
  
    std::cout << "Thrust v" << major << "." << minor << std::endl;

    thrust::host_vector<int> H(4);
    for (auto h:H) {
        std::cout << "> " << h << std::endl;
    }
    thrust::device_vector<int> d_H(4);
    d_H[0] = 1;
    for (auto h:d_H) {
        std::cout << ">> " << h << std::endl;
    }

}
#endif
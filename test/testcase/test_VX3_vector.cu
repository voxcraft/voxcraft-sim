#include "gtest/gtest.h"

#include "Utils/VX3_vector.cuh"
#include <cuda_runtime.h>


//Basic_Usage: create in host, pass to device, copy to another memory in device, return to host
__global__ void kernel_copy_to_output(VX3_hdVector<double> d_vec, int length, double *output) {
    for (int i=0;i<length;i++) {
        output[i] = d_vec[i];
    }
}
TEST(VX3_hdVector_Test, Basic_Usage) {
    std::vector<double> h_vec;
    h_vec.push_back(3.0f);
    double *d_output;
    cudaMalloc((void **)&d_output, h_vec.size()*sizeof(double));

    VX3_hdVector<double> d_vec(h_vec);
    kernel_copy_to_output<<<1,1>>>(d_vec, h_vec.size(), d_output);
    double *h_output;
    cudaHostAlloc((void **)&h_output, h_vec.size()*sizeof(double), cudaHostAllocWriteCombined);
    cudaMemcpy(h_output, d_output, h_vec.size()*sizeof(double), cudaMemcpyDeviceToHost);
    EXPECT_EQ(h_output[0],3.0f);
    EXPECT_NE(h_output[0],0.0f);

    cudaFree(h_output);
    cudaFree(d_output);
}

//Basic_Usage_Larger_Size: create in host, pass to device, copy to another memory in device, return to host
TEST(VX3_hdVector_Test, Basic_Usage_Larger_Size) {
    std::vector<double> h_vec;
    for (int i=0; i<100000; i++) {
        h_vec.push_back(1/((double)i));
    }
    double *d_output;
    cudaMalloc((void **)&d_output, h_vec.size()*sizeof(double));

    VX3_hdVector<double> d_vec(h_vec);
    kernel_copy_to_output<<<1,1>>>(d_vec, h_vec.size(), d_output);
    double *h_output;
    cudaHostAlloc((void **)&h_output, h_vec.size()*sizeof(double), cudaHostAllocWriteCombined);
    cudaMemcpy(h_output, d_output, h_vec.size()*sizeof(double), cudaMemcpyDeviceToHost);
    
    for (int i=0;i<100000;i++) {
        EXPECT_EQ(h_output[i], 1/((double)i));
    }

    cudaFree(h_output);
    cudaFree(d_output);
}

//pure device vector Basic usage
__global__ void kernel_create_and_copy_to_output(int num, double* output){
    VX3_dVector<double> vec;
    for (int i=0; i<num; i++) {
        vec.push_back(1/((double)i));
    }
    for (int i=0; i<num; i++) {
        output[i] = vec[i];
    }
}
TEST(VX3_dVector_Test, Basic_Usage) {
    double *d_output;
    int test_num = 10;
    cudaMalloc((void **)&d_output, test_num * sizeof(double));
    kernel_create_and_copy_to_output<<<1,1>>>(test_num, d_output);
    double *h_output;
    cudaHostAlloc((void **)&h_output, test_num * sizeof(double), cudaHostAllocWriteCombined);
    cudaMemcpy(h_output, d_output, test_num*sizeof(double), cudaMemcpyDeviceToHost);
    for (int i=0;i<test_num;i++) {
        EXPECT_EQ(h_output[i], 1/((double)i));
    }
    cudaFree(h_output);
    cudaFree(d_output);
}
TEST(VX3_dVector_Test, Basic_Usage_Mid_Size) {
    double *d_output;
    int test_num = 100000;
    cudaMalloc((void **)&d_output, test_num * sizeof(double));
    kernel_create_and_copy_to_output<<<1,1>>>(test_num, d_output);
    double *h_output;
    cudaHostAlloc((void **)&h_output, test_num * sizeof(double), cudaHostAllocWriteCombined);
    cudaMemcpy(h_output, d_output, test_num*sizeof(double), cudaMemcpyDeviceToHost);
    for (int i=0;i<test_num;i++) {
        EXPECT_EQ(h_output[i], 1/((double)i));
    }
    cudaFree(h_output);
    cudaFree(d_output);
}

// //Larger version with larger default memory
// __global__ void kernel_create_and_copy_to_output_larger(int num, double* output){
//     VX3_dVector_Larger<double> vec;
//     for (int i=0; i<num; i++) {
//         if (!vec.push_back(1/((double)i))) break;
//     }
//     printf("vec.size() %d\n", vec.size());
//     for (int i=0; i<vec.size(); i++) {
//         output[i] = vec[i];
//     }
// }
// TEST(VX3_dVector_Test, Basic_Usage_Larger_version) {
//     double *d_output;
//     int test_num = 10;
//     cudaMalloc((void **)&d_output, test_num * sizeof(double));
//     kernel_create_and_copy_to_output_larger<<<1,1>>>(test_num, d_output);
//     double *h_output;
//     cudaHostAlloc((void **)&h_output, test_num * sizeof(double), cudaHostAllocWriteCombined);
//     cudaMemcpy(h_output, d_output, test_num*sizeof(double), cudaMemcpyDeviceToHost);
//     for (int i=0;i<test_num;i++) {
//         EXPECT_EQ(h_output[i], 1/((double)i));
//     }
//     cudaFree(h_output);
//     cudaFree(d_output);
// }
// TEST(VX3_dVector_Test, Basic_Usage_Mid_Size_Larger_version) {
//     double *d_output;
//     int test_num = 100000;
//     cudaMalloc((void **)&d_output, test_num * sizeof(double));
//     kernel_create_and_copy_to_output_larger<<<1,1>>>(test_num, d_output);
//     double *h_output;
//     cudaHostAlloc((void **)&h_output, test_num * sizeof(double), cudaHostAllocWriteCombined);
//     cudaMemcpy(h_output, d_output, test_num*sizeof(double), cudaMemcpyDeviceToHost);
//     for (int i=0;i<test_num;i++) {
//         EXPECT_EQ(h_output[i], 1/((double)i));
//     }
//     cudaFree(h_output);
//     cudaFree(d_output);
// }


//test clear 
__global__ void kernel_test_clear(int test_num, double *output) {
    VX3_dVector<double> vec;
    for (int i=0;i<test_num/2;i++) {
        vec.push_back(0.1f);
    }
    vec.clear();
    printf("after clear, vec.size() %d\n", vec.size());
    for (int i=0; i<test_num; i++) {
        if (!vec.push_back(10/((double)i))) break;
    }
    printf("after push, vec.size() %d\n", vec.size());
    for (int i=0; i<vec.size(); i++) {
        output[i] = vec[i];
    }
}
TEST(VX3_dVector_Test, Clear) {
    double *d_output;
    int test_num = 100000;
    double *h_output;
    cudaMalloc((void **)&d_output, test_num * sizeof(double));
    kernel_test_clear<<<1,1>>>(test_num, d_output);
    cudaHostAlloc((void **)&h_output, test_num * sizeof(double), cudaHostAllocWriteCombined);
    cudaMemcpy(h_output, d_output, test_num*sizeof(double), cudaMemcpyDeviceToHost);
    for (int i=0;i<test_num;i++) {
        EXPECT_EQ(h_output[i], 10/((double)i));
    }
    cudaFree(h_output);
    cudaFree(d_output);
}
// // Please use this function with nvidia-smi watching for memory leaking.
// TEST(VX3_dVector_Test, Clear_test_for_leak) { 
//     double *d_output;
//     int test_num = 100000;
//     double *h_output;
//     for (int i=0;i<100;i++) {
//         cudaMalloc((void **)&d_output, test_num * sizeof(double));
//         kernel_test_clear<<<1,1>>>(test_num, d_output);
//         cudaHostAlloc((void **)&h_output, test_num * sizeof(double), cudaHostAllocWriteCombined);
//         cudaMemcpy(h_output, d_output, test_num*sizeof(double), cudaMemcpyDeviceToHost);
//         for (int i=0;i<test_num;i++) {
//             EXPECT_EQ(h_output[i], 10/((double)i));
//         }
//         cudaFree(h_output);
//         cudaFree(d_output);
//     }
// }


//test has
__global__ void kernel_test_has(int total_num, int target, int * output) {
    VX3_dVector<double> d_vec;
    for (int i=0;i<total_num;i++) {
        d_vec.push_back( 1/ ((double) i));
    }
    *output = d_vec.has( 1/((double) target) );
}
TEST(VX3_dVector_Test, Has) {
    int output;
    int * d_output;
    cudaMalloc((void **)&d_output, sizeof(int));

    kernel_test_has<<<1,1>>>(3, 2, d_output);
    cudaMemcpy(&output, d_output, sizeof(int), cudaMemcpyDeviceToHost);
    EXPECT_EQ(output, true);

    kernel_test_has<<<1,1>>>(100, 20, d_output);
    cudaMemcpy(&output, d_output, sizeof(int), cudaMemcpyDeviceToHost);
    EXPECT_EQ(output, true);

    kernel_test_has<<<1,1>>>(100, 120, d_output);
    cudaMemcpy(&output, d_output, sizeof(int), cudaMemcpyDeviceToHost);
    EXPECT_EQ(output, false);

    kernel_test_has<<<1,1>>>(3, 20, d_output);
    cudaMemcpy(&output, d_output, sizeof(int), cudaMemcpyDeviceToHost);
    EXPECT_EQ(output, false);
}

//test default in device construction
struct Test_TI_Material {
    __device__ Test_TI_Material() {
        d_strainData.clear();
        d_strainData.push_back(0.0f);
    }
    __device__ Test_TI_Material(Test_TI_Material& In) {
        d_strainData = In.d_strainData;
    }
    Test_TI_Material(std::vector<double> &h_vec) :
    hd_strainData(h_vec) {
        
    }
    __device__ void sync_hdVector_to_dVector() {
        d_strainData.clear();
        d_strainData.push_back(0.0f);

        for (unsigned i=0;i<hd_strainData.size();i++) {
            d_strainData.push_back(hd_strainData[i]);
        }
    }
    VX3_hdVector<double> hd_strainData;
    VX3_dVector<double> d_strainData;
};
__global__ void kernel_test_construction( double* output ) {
    Test_TI_Material a;
    a.d_strainData.push_back(0.2f);
    output[0] = (double)a.d_strainData.size();
    output[1] = a.d_strainData[0];
    output[2] = a.d_strainData[1];
}
TEST(VX3_dVector_Test, Default_Construction) {
    double* d_output;
    double h_output[3];
    cudaMalloc((void **)&d_output, 3*sizeof(double));
    kernel_test_construction<<<1,1>>>(d_output);
    cudaMemcpy(&h_output[0], d_output, 3*sizeof(double), cudaMemcpyDeviceToHost);
    EXPECT_EQ(h_output[0], 2.0f);
    EXPECT_EQ(h_output[1], 0.0f);
    EXPECT_EQ(h_output[2], 0.2f);
}

//test construct in host, pass to device and push
__global__ void kernel_get_and_use( Test_TI_Material* a, double* output ) {
    a->sync_hdVector_to_dVector();
    for (int i=0;i<100;i++) {
        a->d_strainData.push_back(0.32f);
    }
    output[0] = (double)a->d_strainData.size();
    output[1] = a->d_strainData[1];
    output[2] = a->d_strainData[102];
}
TEST(VX3_dVector_Test, Construction_in_host_then_pass_to_device_and_use) {
    double* d_output;
    double h_output[3];
    cudaMalloc((void **)&d_output, 3*sizeof(double));
    std::vector<double> tmp_v;
    for (int i=0;i<100;i++) {
        tmp_v.push_back(1.29f);
    }

    Test_TI_Material tmp_a(tmp_v);

    Test_TI_Material* d_a;
    cudaMalloc((void **)&d_a, sizeof(Test_TI_Material));
    cudaMemcpy(d_a, &tmp_a, sizeof(Test_TI_Material), cudaMemcpyHostToDevice);
    kernel_get_and_use<<<1,1>>>(d_a, d_output);
    cudaMemcpy(&h_output[0], d_output, 3*sizeof(double), cudaMemcpyDeviceToHost);
    EXPECT_EQ(h_output[0], 201.0f);
    EXPECT_EQ(h_output[1], 1.29f);
    EXPECT_EQ(h_output[2], 0.32f);
}



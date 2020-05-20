#include "gtest/gtest.h"

#include "VX3_queue.cuh"

__global__ void dqueue_basic(int* error, int method=0) {
    VX3_dQueue<int> q;
    int test_length = 1000;
    if (method==1) {
        q.push_back(219);
        if (q.front()!=219) {
            *error=1;
            printf(COLORCODE_BOLD_RED "ERROR: front() is not right.\n");
        }
        if (q.pop_front()!=219) {
            *error=1;
            printf(COLORCODE_BOLD_RED "ERROR: pop_front() is not right.\n");
        }        
    }
    for (int i=0;i<test_length;i++) {
        q.push_back(i);
        if (i!=q.back()) {
            *error=1;
            printf(COLORCODE_BOLD_RED "ERROR: back() is not right.\n");
        }
    }
    if (q.size()!=test_length) {
        *error=1;
        printf(COLORCODE_BOLD_RED "ERROR: size\n");
    }
    for (int i=0;i<test_length;i++) {
        if (q.isEmpty()) {
            *error = 1;
            printf(COLORCODE_BOLD_RED "ERROR: is empty.\n");
        }
        auto j = q.pop_front();
        if (i!=j) {
            *error=1;
            printf(COLORCODE_BOLD_RED "ERROR: %d != %d\n", i,j);
        }
    }
    if (q.size()!=0 || !q.isEmpty()) {
        *error=1;
        printf(COLORCODE_BOLD_RED "ERROR: size at end.\n");
    }
}

TEST(VX3_dQueue_Test, Basic_Usage) {
    int *error;
    cudaMallocManaged((void**)&error, sizeof(int));
    *error = 0;
    dqueue_basic<<<1,1>>>(error);
    cudaDeviceSynchronize();
    EXPECT_NE(*error,1);

}

TEST(VX3_dQueue_Test, Basic_Usage_1) {
    int *error;
    cudaMallocManaged((void**)&error, sizeof(int));

    *error = 0;
    dqueue_basic<<<1,1>>>(error,1);
    cudaDeviceSynchronize();
    EXPECT_NE(*error,1);

}
__global__ void dqueue_init(VX3_dQueue<int>* dq) {
    dq->clear();
}
__global__ void dqueue_race(VX3_dQueue<int>* dq) {
    dq->push_back(threadIdx.x);
}
__device__ void swap(int *xp, int *yp) 
{ 
    int temp = *xp; 
    *xp = *yp; 
    *yp = temp; 
} 
__global__ void dqueue_race_check(int* error, VX3_dQueue<int>* dq) {
    int n=1000;
    int arr[1000];
    for (int i=0;i<1000;i++) {
        arr[i] = dq->pop_front();
    }
    int i, j; 
    for (i = 0; i < n-1; i++)       
        for (j = 0; j < n-i-1; j++)  
            if (arr[j] > arr[j+1]) 
                swap(&arr[j], &arr[j+1]);
    
    for (int i=0;i<1000;i++) {
        if (arr[i] != i) {
            *error = 1;
            printf(COLORCODE_BOLD_RED "ERROR: difference after push, pop, and bubble sort.\n");
        }
    }
}

TEST(VX3_dQueue_Test, Race_Condition) {
    int *error;
    cudaMallocManaged((void**)&error, sizeof(int));

    VX3_dQueue<int>* dq;
    cudaMallocManaged((void**)&dq, sizeof(VX3_dQueue<int>));
    dqueue_init<<<1,1>>>(dq);
    cudaDeviceSynchronize();
    
    dqueue_race<<<1,1000>>>(dq);
    cudaDeviceSynchronize();

    *error = 0;
    dqueue_race_check<<<1,1>>>(error,dq);
    cudaDeviceSynchronize();
    EXPECT_NE(*error,1);

}

__global__ void dqueue_push_pop(int * error) {
    VX3_dQueue<int> q;
    int test_length = 1000;
    for (int i=0;i<test_length;i++) {
        q.push_back(i);
        int j = q.pop_front();
        if (i!=j) {
            *error = 1;
            printf(COLORCODE_BOLD_RED "ERROR: empty queue with pop_front after push_back.\n");
            break;
        }
    }
}

TEST(VX3_dQueue_Test, push_pop_function) {
    int *error;
    cudaMallocManaged((void**)&error, sizeof(int));

    *error = 0;
    dqueue_push_pop<<<1,1>>>(error);
    cudaDeviceSynchronize();
    EXPECT_NE(*error,1);
}

__global__ void dqueue_back(int * error) {
    VX3_dQueue<int> q;
    int test_length = 1000;
    for (int i=0;i<test_length;i++) {
        q.push_back(i);
        int m = q.back();
        if (i!=m) {
            *error = 1;
            printf(COLORCODE_BOLD_RED "ERROR: back after push_back.\n");
            break;
        }
        int j = q.pop_front();
        if (i!=j) {
            *error = 1;
            printf(COLORCODE_BOLD_RED "ERROR: empty queue with pop_front after push_back in back.\n");
            break;
        }
    }
}

TEST(VX3_dQueue_Test, back_function) {
    int *error;
    cudaMallocManaged((void**)&error, sizeof(int));

    *error = 0;
    dqueue_back<<<1,1>>>(error);
    cudaDeviceSynchronize();
    EXPECT_NE(*error,1);
}
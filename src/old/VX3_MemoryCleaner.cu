#include "VX3_MemoryCleaner.h"
// #include "VX3_TaskManager.h"

bool VX3_MemoryCleaner_running=true;
boost::mutex MemoryCleaner_mutex;
std::vector<void *> MemoryCleaner_toBeFreedCUDAPointer;

void VX3_MemoryCleaner::operator()() {
    printf("Start. %d\n", VX3_MemoryCleaner_running);
    // while(VX3_MemoryCleaner_running) {
    //     // printf("Thread count: %ld\n", TaskManager_all_threads.size());
    //     if (TaskManager_all_threads.size()>0) {
    //          // try_join_for
    //         for (int i=0;i<TaskManager_all_threads.size();i++) {
    //             if (TaskManager_all_threads[i].try_join_for( boost::chrono::nanoseconds(1) )) {
    //                 // printf("clean one thread.\n");
    //                 TaskManager_all_threads_mutex.lock();
    //                 TaskManager_all_threads.erase(TaskManager_all_threads.begin()+i);
    //                 TaskManager_all_threads_mutex.unlock();
    //             }
    //         }
    //     }
    //     if (TaskManager_all_threads.size()==0) {
    //         if (MemoryCleaner_toBeFreedCUDAPointer.size()>0) {
    //             MemoryCleaner_mutex.lock();
    //             MemoryCleaner_toBeFreedCUDAPointer.clear();
    //             MemoryCleaner_mutex.unlock();
    //             printf("cudaDeviceReset() called! should release all previously used memory.\n");
    //             CUDA_ERROR_CHECK(cudaDeviceReset());
    //         }
    //         // while(MemoryCleaner_toBeFreedCUDAPointer.size()>0) {
    //         //     void* pointer;
    //         //     MemoryCleaner_mutex.lock();
    //         //     pointer = MemoryCleaner_toBeFreedCUDAPointer.back();
    //         //     MemoryCleaner_toBeFreedCUDAPointer.pop_back();
    //         //     MemoryCleaner_mutex.unlock();
    //         //     // printf("Cleaning %p (%ld remain)...\n", pointer, MemoryCleaner_toBeFreedCUDAPointer.size());
    //         //     VcudaFree(pointer);
    //         //     std::this_thread::sleep_for(std::chrono::milliseconds(1));
    //         // }
    //     }
    //     std::this_thread::sleep_for(std::chrono::seconds(1));
    // }

}
#if !defined(VX3_MEMORY_CLEANER_H)
#define VX3_MEMORY_CLEANER_H

#include <stdio.h>
#include <thread>
#include <vector>
#include <boost/thread/mutex.hpp>
#include <boost/chrono.hpp>
#include "VX3.cuh"

extern bool VX3_MemoryCleaner_running;
extern boost::mutex MemoryCleaner_mutex;
extern std::vector<void *> MemoryCleaner_toBeFreedCUDAPointer;

inline void MycudaFree(void * ptr) {
    MemoryCleaner_mutex.lock();
    MemoryCleaner_toBeFreedCUDAPointer.push_back(ptr);
    MemoryCleaner_mutex.unlock();
}
class VX3_MemoryCleaner
{    
public:
    VX3_MemoryCleaner() = default;
  
    void operator()();
};

#endif // VX3_MEMORY_CLEANER_H

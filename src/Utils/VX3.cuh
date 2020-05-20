#if !defined(VX3_H)
#define VX3_H

#include <string>
#include <stdexcept>
#include <chrono>
#include <iomanip>
#include <ctime>
#include <sstream>

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#include <boost/algorithm/string.hpp>

inline std::string u_format_now(std::string format) {
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);

    std::stringstream folderName;
    folderName << std::put_time(std::localtime(&in_time_t), format.c_str());
    return folderName.str();
}

inline bool u_with_ext(fs::path file, std::string ext) {

    std::string ext_file = file.filename().extension().string();
    boost::to_upper(ext);
    boost::to_upper(ext_file);

    return ext==ext_file;
}

#include <curand.h>
#include <curand_kernel.h>
__device__ __inline__ int random(int max, int random_seed=0) {
    curandState_t state;
    curand_init(random_seed, 0, 0, &state);
    return curand(&state) % max;
}

#define COLORCODE_RED "\033[0;31m" 
#define COLORCODE_BOLD_RED "\033[1;31m\n" 
#define COLORCODE_GREEN "\033[0;32m" 
#define COLORCODE_BLUE "\033[0;34m" 
#define COLORCODE_RESET "\033[0m" 

//#define DEBUG_LINE printf("%s(%d): %s\n", __FILE__, __LINE__, u_format_now("at %M:%S").c_str());
#define CUDA_DEBUG_LINE(str) {printf("%s(%d): %s\n", __FILE__, __LINE__, str);}
#define DEBUG_LINE_

#ifndef CUDA_ERROR_CHECK
    __device__ __host__ inline void CUDA_ERROR_CHECK_OUTPUT(cudaError_t code, const char *file, int line, bool abort=false) {
        if (code != cudaSuccess) {
            printf(COLORCODE_BOLD_RED "%s(%d): CUDA Function Error: %s \n" COLORCODE_RESET, file, line, cudaGetErrorString(code));
            if (abort) assert(0);
        }
    }
    #define CUDA_ERROR_CHECK(ans) { CUDA_ERROR_CHECK_OUTPUT((ans), __FILE__, __LINE__); }
#endif
#define VcudaMemGetInfo(a,b) {CUDA_ERROR_CHECK(cudaMemGetInfo(a,b))}
#define VcudaDeviceSetLimit(a,b) {CUDA_ERROR_CHECK(cudaDeviceSetLimit(a,b))}
#define VcudaSetDevice(a) {CUDA_ERROR_CHECK(cudaSetDevice(a))}
#define VcudaGetDeviceCount(a) {CUDA_ERROR_CHECK(cudaGetDeviceCount(a))}
#define VcudaMemcpy(a,b,c,d)  {CUDA_ERROR_CHECK(cudaMemcpy(a,b,c,d))}
#define VcudaMemcpyAsync(a,b,c,d,e)  {CUDA_ERROR_CHECK(cudaMemcpyAsync(a,b,c,d,e))}
#define VcudaMalloc(a,b) {CUDA_ERROR_CHECK(cudaMalloc(a,b))}
#define VcudaFree(a) {CUDA_ERROR_CHECK(cudaFree(a))}
#define VcudaGetLastError() {CUDA_ERROR_CHECK(cudaGetLastError())}
#define VcudaDeviceSynchronize() {CUDA_ERROR_CHECK(cudaDeviceSynchronize())}
#define VcudaMemcpyHostToDevice cudaMemcpyHostToDevice
#define VcudaMemcpyDeviceToHost cudaMemcpyDeviceToHost
#define CUDA_CHECK_AFTER_CALL() {CUDA_ERROR_CHECK(cudaGetLastError());}
#define DEBUG_CUDA_ERROR_CHECK_STATUS() {printf("DEBUG_CUDA_ERROR_CHECK_STATUS ON.\n");}

#include "Utils/VX3_Vec3D.h"
#include "Utils/VX3_Quat3D.h"
#include "Utils/VX3_vector.cuh"
#include "Utils/VX3_dictionary.cuh"
#include "VX3/VX3_SimulationResult.h"

#endif // VX3_H

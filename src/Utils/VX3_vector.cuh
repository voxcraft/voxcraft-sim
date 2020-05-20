#if !defined(VX3_VECTOR_H)
#define VX3_VECTOR_H
#include <cuda_runtime.h>

#include "Utils/VX3.cuh"
#include <vector>

template <typename T>
class VX3_hdVector; // A vector that can be initialized in host and passed to device. It should be freed after return from device.
template <typename T>
class VX3_dVector; // A vector that can be initialized in device, and do push_back() and get(), and it should be freed after use.
template <typename T>
class VX3_dVector_Larger; // If there will be more than 128 elements in the vector, we can use this vector class instead of VX3_dVector.

template <typename T> class VX3_hdVector {
    /* It should be initialized in host, and pass to kernel directly, and it should be freed after kernel returns. */
  public:
    VX3_hdVector<T>() = default;
    VX3_hdVector<T>(const std::vector<T> &p) {
        num_main = p.size();
        VcudaMalloc(&main, num_main * sizeof(T));
        T *temp = (T *)malloc(num_main * sizeof(T));
        // VcudaHostAlloc((void **)&temp, num_main*sizeof(T), cudaHostAllocWriteCombined);
        for (unsigned i = 0; i < num_main; i++) {
            temp[i] = p[i];
        }
        VcudaMemcpy(main, temp, num_main * sizeof(T), VcudaMemcpyHostToDevice);
        // cudaFreeHost(temp);
        delete temp;
    }
    void free() {
        // cannot use ~ method, becuase passing parameter to kernel triggers ~ method.
        VcudaFree(main);
    }
    __device__ inline T &operator[](unsigned index) { return main[index]; }
    __device__ inline T get(unsigned index) { return main[index]; }
    __device__ inline unsigned size() { return num_main; }
    unsigned num_main;
    T *main;
};

#if !defined(DEFAULT_CHUNK_SIZE)
#define DEFAULT_CHUNK_SIZE 64
#endif // Shared with VX3_queue

template <typename T> class VX3_dVector {
    /* It should be initialized in dev, and do push_back() and get(), and it should be freed after use.
       if size<Default, use GPU local memory memory[Default], otherwise use malloc and main and store things in GPU global memory. (malloc
       is slow in GPU)*/
  public:
    __device__ VX3_dVector<T>() { clear(); } // Never called, since we copy the mem to GPU.
    __device__ ~VX3_dVector<T>() {
        if (main)
            delete main;
    }
    __device__ bool inline push_back(T t) {
        // Critical area Start
        bool leave = true;
        while (leave) {
            if (atomicCAS(mutex, 0, 1) == 0) {
                // Critical area Body

                if (num_main < sizeof_chunk) {
                    if (main) { // use main
                        main[num_main] = t;
                    } else { // use memory
                        default_memory[num_main] = t;
                    }
                } else { // need allocation
                    T *new_main;
                    new_main = (T *)malloc(sizeof_chunk * 2 * sizeof(T));
                    if (new_main == NULL) {
                        printf("Out of memory when alloc %ld bytes.\n", sizeof_chunk * 2 * sizeof(T));
                        return false;
                    }
                    if (main) {
                        memcpy(new_main, main, sizeof_chunk * sizeof(T));
                        delete main;
                    } else {
                        memcpy(new_main, default_memory, sizeof_chunk * sizeof(T));
                    }
                    main = new_main;
                    main[num_main] = t;
                    sizeof_chunk *= 2;
                }
                num_main += 1;

                // Critical area Body End
                leave = false;
                atomicExch(mutex, 0);
            }
        }
        // Critical area End
        return true;
    }
    // __device__ T inline pop_back() {
    //     T ret;
    //     num_main -= 1;
    //     if (main)
    //         ret = main[num_main];
    //     else
    //         ret = default_memory[num_main];
    //
    //     return ret;
    // }
    __device__ inline VX3_dVector<T> &operator=(const VX3_dVector<T> &vIn) { // TODO: Unit Test and copy to large
        if (vIn.main) {
            sizeof_chunk = vIn.sizeof_chunk;
            num_main = vIn.num_main;
            main = (T *)malloc(sizeof_chunk * sizeof(T));
            for (unsigned i = 0; i < num_main; i++) {
                main[i] = vIn.main[i];
            }
        } else {
            sizeof_chunk = vIn.sizeof_chunk;
            num_main = vIn.num_main;
            for (unsigned i = 0; i < num_main; i++) {
                default_memory[i] = vIn.default_memory[i];
            }
        }
        return *this;
    }
    __device__ inline unsigned size() { return num_main; };
    __device__ inline T &operator[](unsigned index) { return get(index); };
    __device__ inline T &get(unsigned index) {
        if (main)
            return main[index];
        else
            return default_memory[index];
    }
    __device__ void clear() {
        if (!mutex) {
            mutex = (int *)malloc(sizeof(int));
            *mutex = 0;
        }
        num_main = 0;
        if (main) {
            delete main;
            main = NULL;
        }
        sizeof_chunk = DEFAULT_CHUNK_SIZE;
    }
    __device__ inline bool has(T value) {
        if (main) {
            for (unsigned i = 0; i < num_main; i++) {
                if (main[i] == value) {
                    return true;
                }
            }
            return false;
        } else {
            for (unsigned i = 0; i < num_main; i++) {
                if (default_memory[i] == value) {
                    return true;
                }
            }
            return false;
        }
    }
    T *main = NULL;
    T default_memory[DEFAULT_CHUNK_SIZE];
    unsigned sizeof_chunk;
    unsigned num_main;
    int *mutex = NULL;
};

// // Just copy VX_dVector and make it larger
// #define DEFAULT_CHUNK_SIZE_LARGER 2048
// template <typename T> class VX3_dVector_Larger {
//     /* It should be initialized in dev, and do push_back() and get(), and it should be freed after use.
//        if size<Default, use GPU local memory memory[Default], otherwise use malloc and main and store things in GPU global memory.
//        (malloc is slow in GPU)*/
//   public:
//     __device__ VX3_dVector_Larger<T>() { clear(); }
//     __device__ ~VX3_dVector_Larger<T>() {
//         if (main)
//             delete main;
//     }
//     __device__ bool inline push_back(T t) {
//         if (num_main < sizeof_chunk) {
//             if (main) { // use main
//                 main[num_main] = t;
//             } else { // use memory
//                 default_memory[num_main] = t;
//             }
//         } else { // need allocation
//             T *new_main;
//             new_main = (T *)malloc(sizeof_chunk * 2 * sizeof(T));
//             if (new_main == NULL) {
//                 printf("Out of memory when alloc %ld bytes.\n", sizeof_chunk * 2 * sizeof(T));
//                 return false;
//             }
//             if (main) {
//                 memcpy(new_main, main, sizeof_chunk * sizeof(T));
//                 delete main;
//             } else {
//                 memcpy(new_main, default_memory, sizeof_chunk * sizeof(T));
//             }
//             main = new_main;
//             main[num_main] = t;
//             sizeof_chunk *= 2;
//         }
//         num_main += 1;
//         return true;
//     }
//     __device__ inline unsigned size() { return num_main; };
//     __device__ inline T &operator[](unsigned index) { return get(index); };
//     __device__ inline T &get(unsigned index) {
//         if (main)
//             return main[index];
//         else
//             return default_memory[index];
//     }
//     __device__ void clear() {
//         num_main = 0;
//         if (main) {
//             delete main;
//             main = NULL;
//         }
//         sizeof_chunk = DEFAULT_CHUNK_SIZE_LARGER;
//     }
//     __device__ inline bool has(T value) {
//         if (main) {
//             for (unsigned i = 0; i < num_main; i++) {
//                 if (main[i] == value) {
//                     return true;
//                 }
//             }
//             return false;
//         } else {
//             for (unsigned i = 0; i < num_main; i++) {
//                 if (default_memory[i] == value) {
//                     return true;
//                 }
//             }
//             return false;
//         }
//     }
//     T *main = NULL;
//     T default_memory[DEFAULT_CHUNK_SIZE_LARGER];
//     unsigned sizeof_chunk;
//     unsigned num_main;
// };

#endif // VX3_VECTOR_H

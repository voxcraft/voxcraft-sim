#if !defined(VX3_QUEUE)
#define VX3_QUEUE

#include <cuda_runtime.h>

#include "Utils/VX3.cuh"
#include <vector>

template <typename T>
class VX3_dQueue; // A vector that can be initialized in device, and do push_back() and get(), and it should be freed after use.

#if !defined(DEFAULT_CHUNK_SIZE)
#define DEFAULT_CHUNK_SIZE 64
#endif // Shared with VX3_vector

template <typename T> class VX3_dQueue {
  public:
    __device__ VX3_dQueue<T>() { clear(); } // Never called, since we copy the mem to GPU.
    __device__ ~VX3_dQueue<T>() {
        if (main)
            delete main;
    }

    __device__ void clear() {
        if (!mutex) {
            mutex = (int *)malloc(sizeof(int));
            *mutex = 0;
        }
        cursor_front = 0;
        cursor_back = 0;
        if (main) {
            delete main;
            main = NULL;
        }
        sizeof_chunk = DEFAULT_CHUNK_SIZE;
    }

    __device__ inline unsigned size() {
        if (cursor_back >= cursor_front)
            return cursor_back - cursor_front;
        else
            return cursor_back - cursor_front + sizeof_chunk;
    }
    __device__ inline int step(int current) {
        int next = current + 1;
        if (next >= sizeof_chunk)
            next -= sizeof_chunk;
        return next;
    }
    __device__ bool inline push_back(T t) {
        // Critical area Start
        bool leave = true;
        while (leave) {
            if (atomicCAS(mutex, 0, 1) == 0) {
                // Critical area Body

                if (cursor_back != cursor_front - 1 && cursor_back != cursor_front + sizeof_chunk - 1) {
                    if (main) { // use main
                        main[cursor_back] = t;
                    } else { // use memory
                        default_memory[cursor_back] = t;
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
                        memcpy(new_main + sizeof_chunk * sizeof(T), main, sizeof_chunk * sizeof(T));
                        delete main;
                    } else {
                        memcpy(new_main, default_memory, sizeof_chunk * sizeof(T));
                        memcpy(new_main + sizeof_chunk * sizeof(T), default_memory, sizeof_chunk * sizeof(T));
                    }
                    main = new_main;

                    if (cursor_back < cursor_front)
                        cursor_back += sizeof_chunk;
                    main[cursor_back] = t;
                    sizeof_chunk *= 2;
                }
                cursor_back = step(cursor_back);

                // Critical area Body End
                leave = false;
                atomicExch(mutex, 0);
            }
        }
        // Critical area End
        return true;
    }

    __device__ inline T pop_front() {
        T ret;
        // Critical area Start
        bool leave = true;
        while (leave) {
            if (atomicCAS(mutex, 0, 1) == 0) {
                // Critical area Body

                ret = front();
                cursor_front = step(cursor_front);

                // Critical area Body End
                leave = false;
                atomicExch(mutex, 0);
            }
        }
        // Critical area End
        return ret;
    }
    __device__ bool isEmpty() { return size() == 0; }
    __device__ inline T front() {
        if (main) {
            return main[cursor_front];
        } else {
            return default_memory[cursor_front];
        }
    }
    __device__ inline T back() {
        if (isEmpty())
            return (T)0;
        int cursor = cursor_back - 1;
        if (cursor < 0)
            cursor += sizeof_chunk;  
        if (main) {
            return main[cursor];
        } else {
            return default_memory[cursor];
        }
    }
    T *main = NULL;
    T default_memory[DEFAULT_CHUNK_SIZE];
    unsigned sizeof_chunk;
    unsigned cursor_back, cursor_front;
    int *mutex = NULL;
};

#endif // VX3_QUEUE

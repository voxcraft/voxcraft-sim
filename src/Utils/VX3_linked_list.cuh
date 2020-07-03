//
// Created by Sida Liu, not finished.
//  This class implements a linked list, with push_back() and size() and remove().
//
#if !defined(VX3_LINKED_LIST_H)
#define VX3_LINKED_LIST_H

#if !defined(DEFAULT_CHUNK_SIZE)
#define DEFAULT_CHUNK_SIZE 64
#endif // Shared with VX3_queue

template <typename T> class VX3_linkedList {
public:
    T *main = NULL;
    T *main_ptr = NULL;
    T default_memory[DEFAULT_CHUNK_SIZE];
    t default_memory_ptr[DEFAULT_CHUNK_SIZE];
    unsigned sizeof_chunk;
    unsigned num_main;
    int mutex = 0;

public:
    __device__ VX3_linkedList<T>() { clear(); };
    __device__ ~VX3_linkedList<T>() {
        if (main)
            free(main);
        if (main_ptr)
            free(main_ptr);
    }
    __device__ void push_back(T t) {
        // Critical area Start
        bool leave = true;
        while (leave) {
            if (atomicCAS(&mutex, 0, 1) == 0) {
                // Critical area Body



                // Critical area Body End
                leave = false;
                atomicExch(&mutex, 0);
            }
            
    }
    __device__ void clear(){
        mutex = 0;
        num_main = 0;
        if (main) {
            free(main);
            main = NULL;
        }
        if (main_ptr) {
            free(main_ptr);
            main_ptr = NULL;
        }
        sizeof_chunk = DEFAULT_CHUNK_SIZE;
    };
};

#endif // VX3_LINKED_LIST_H

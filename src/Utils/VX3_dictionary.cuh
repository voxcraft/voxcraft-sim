// Implementation of double hashing dictionary. Refer to https://en.wikipedia.org/wiki/Double_hashing and https://ocw.mit.edu/courses/electrical-engineering-and-computer-science/6-006-introduction-to-algorithms-fall-2011/lecture-videos/lecture-10-open-addressing-cryptographic-hashing/
#if !defined(VX3_DICTIONARY_H)
#define VX3_DICTIONARY_H
#include <assert.h>
#include <curand.h>
#include <curand_kernel.h>
enum VX3_dictionary_flag {
    NOT_USED = 0,
    USED = 1,
    DELETED = 2,
};
#define DEFAULT_MEMORY_SIZE 16
#define LG_DEFAULT_MEMORY_SIZE 4
#define BIG_PRIME_NUMBER 2147483647 // which is 0x7FFFFFFF, I don't change prime number everytime, hope it will be ok.
#define ALPHA_RATIO 0.6 // Prof. Srini Devadas says alpha is dangous when >0.5, but I think GPU memory is also very precious.

template <typename T, typename U>
class VX3_dDictionary {
public:
    // Initialization. seed is for testing purpose.
    __device__ __host__ __inline__ VX3_dDictionary(){}
    __device__ __inline__ void clear(int seed=0) {
        random_seed = seed;
        num_elements = 0;
        num_slots = DEFAULT_MEMORY_SIZE;
        lg_of_slots = LG_DEFAULT_MEMORY_SIZE;
        num_max_elements_of_this_slots = (int) num_slots * ALPHA_RATIO;
        
        random_parameter[0] = _random(BIG_PRIME_NUMBER);
        random_parameter[1] = _random(BIG_PRIME_NUMBER);
        random_parameter[2] = _random(BIG_PRIME_NUMBER);
        random_parameter[3] = _random(BIG_PRIME_NUMBER);
        for (int i=0;i<num_slots;i++) {
            main_flag[i] = NOT_USED;
        }
        extended = false;
        if (extend_flag) delete extend_flag;
        extend_flag = NULL;
        if (extend_key) delete extend_key;
        extend_key = NULL;
        if (extend_value) delete extend_value;
        extend_value = NULL;
    }
    // Key Operation 1: set
    __device__ __inline__ void set(T t, U u) {
        expandMemoryIfNeeded();
        num_elements ++;
        if (!extended) {
            return _set(t, u, main_flag, main_key, main_value);
        } else {
            return _set(t, u, extend_flag, extend_key, extend_value);
        }
    }
    // Key Operation 2: get
    __device__ __inline__ U get(T t) {
        if (!extended) {
            return _get(t, main_flag, main_key, main_value);
        } else {
            return _get(t, extend_flag, extend_key, extend_value);
        }
    }
    __device__ __inline__ void _set(T t, U u, short* flag, T* key, U* value) {
        int num_trail = 0;
        int i;
        for (i=0;i<num_slots;i++) {
            size_t h = double_hashing(t, i);
            if (flag[h] == NOT_USED || flag[h] == DELETED) {
                key[h] = t;
                value[h] = u;
                flag[h] = USED;
                break;
            }
            num_trail ++;
        }
        if (i==num_slots) { printf("ERROR: VX3_dictionary: memory overflow. should expand earlier.\n"); }
        // printf("insert on %d trail.\n", num_trail);
    }
    __device__ __inline__ U _get(T t, short* flag, T* key, U* value) {
        int i;
        for (i=0;i<num_slots;i++) {
            size_t h = double_hashing(t, i);
            switch(flag[h]) {
                case NOT_USED:
                    return (U)-1;
                case USED:
                    if (key[h]==t) {
                        return value[h];
                    }
                    continue;
                case DELETED:
                    continue;
            }
        }
        return (U)-1;
    }
    // KEY PART: Double Hashing. (MIT OCW 6.006 Lecture 10)
    __device__ __inline__ int double_hashing(T t, int trail) {
        return (hash1(t) + trail*hash2(t)) % num_slots;
    }
    // Hash2 should never yield an index of 0, should be relatively prime to num_slots.
    __device__ __inline__ int hash2(T t) {
        return 2*hash(t, 1)+1;
    }
    // Hash1 is an ordinary Universal Hash with certain parameters
    __device__ __inline__ int hash1(T t) {
        return hash(t, 0);
    }
    // Ordinary Universal Hash
    __device__ __inline__ int hash(T t, int param) {
        unsigned long k = (unsigned long) t;
        k = (random_parameter[param*2]*k + random_parameter[param*2+1]);
        k = k % BIG_PRIME_NUMBER;
        return  k % num_slots;
    }
    // Expand the memory using Double Tabling (MIT OCW 6.006 Lecture 9)
    __device__ __inline__ void expandMemoryIfNeeded() {
        if (num_elements < num_max_elements_of_this_slots) return;
        // Saving current parameters for temperary use.
        size_t old_num_slots = num_slots;
        short* old_flag;
        T* old_key;
        U* old_value;
        if (extended) {
            old_flag = extend_flag;
            old_key = extend_key;
            old_value = extend_value;
        } else {
            old_flag = main_flag;
            old_key = main_key;
            old_value = main_value;
        }
        // Allocate fresh new memory.
        short* new_flag = (short*)malloc(2*num_slots*sizeof(short));
        T* new_key = (T*)malloc(2*num_slots*sizeof(T));
        U* new_value = (U*)malloc(2*num_slots*sizeof(U));
        if (new_flag==NULL || new_key==NULL || new_value==NULL) {
            printf("ERROR: Out of memory.\n");
            return;
        }
        // Change parameters.
        extended = true;
        lg_of_slots++;
        num_slots=1;
        for (int i=0;i<lg_of_slots;i++)
            num_slots *= 2;
        num_max_elements_of_this_slots = (int) num_slots * ALPHA_RATIO;
        for (int i=0;i<num_slots;i++) {
            new_flag[i] = NOT_USED;
        }
        // Copy key value pairs from old memory to new memory. Need to do hash again since num_slots changed.
        for (int i=0;i<old_num_slots;i++) {
            if (old_flag[i]==USED) {
                _set(old_key[i], old_value[i], new_flag, new_key, new_value);
            }
        }
        // Redirect pointer
        if (extend_flag) delete extend_flag;
        extend_flag = new_flag;
        if (extend_key) delete extend_key;
        extend_key = new_key;
        if (extend_value) delete extend_value;
        extend_value = new_value;
    }
    // print all elements for debugging purpose
    __device__ void print() {
        for (int i=0;i<num_slots;i++) {
            if (!extended) {
                if (main_flag[i]==USED) {
                    printf("%d) %d; \n", i, (int)main_value[i]);
                }
            } else {
                if (extend_flag[i]==USED) {
                    printf("%d) %d; \n", i, (int)extend_value[i]);
                }
            }
        }
    }
    // use cuda Math API to generate a random number
    __device__ __inline__ int _random(int max) {
        curandState_t state;
        curand_init(random_seed, 0, 0, &state);
        return curand(&state) % max;
    }
    long random_parameter[4];
    short main_flag[DEFAULT_MEMORY_SIZE]; //-1: Deleted; 0: Not used; 1: Used.
    T main_key[DEFAULT_MEMORY_SIZE];
    U main_value[DEFAULT_MEMORY_SIZE];
    short* extend_flag = NULL;
    T* extend_key = NULL;
    U* extend_value = NULL;
    size_t num_slots;
    size_t lg_of_slots;
    size_t num_max_elements_of_this_slots;
    size_t num_elements;
    int random_seed;
    bool extended;

};

#endif // VX3_DICTIONARY_H

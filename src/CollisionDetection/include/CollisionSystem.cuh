//
// Created by David Matthews on 5/21/20.
//

#ifndef COLLISIONDETECTION_COLLISIONSYSTEM_CUH
#define COLLISIONDETECTION_COLLISIONSYSTEM_CUH


#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/logical.h>
#include <thrust/random.h>
#include <thrust/sequence.h>

#include "cub/cub.cuh"

#include "BoundingBox.cuh"
#include "Collisions.cuh"
#include "Morton.cuh"


class CollisionSystem {

public:
    CollisionSystem() :
            CollisionSystem(2, 1, true) {}

    CollisionSystem(size_t n, size_t maxCollisionsPerObject,  bool toAllocHostVecs);

    ~CollisionSystem();

    unsigned int get_num_masses() { return N; }

    unsigned int get_max_collisions_per_mass() { return MAX_COLLISIONS_PER_OBJECT; }

    /**
     * Resize the number of masses in the system
     *
     * NOTE: if we did not allocate host vectors, we will not resize them. we will assume that you have resized them for us!.
     * @param n new size
     */
    void set_num_masses(size_t n);

    /**
     * reserve space for new masses in the system so that future calls to set_num_masses are faster.
     * Does not shrink reservations if we already have memory allocated.
     *
     * NOTE: if we did not allocate host vectors, we will not resize them. we will assume that you have resized them for us!.
     * @param n new reserved size.
     */
    void set_reserved_num_masses(size_t n);

    /**
     * Change how many collisions per mass are supported by the system.
     * @param m new max number of collisions per mass.
     */
    void set_max_num_cols_per_mass(size_t m);

    void update_device_pointers_and_functors();

    /**
     * Get a pointer to the x pos host array.
     * @return pointer to first element in x_pos host array.
     */
    float *get_x_pos_host() { return x_pos_h; }

    void set_x_pos_host(float *x_pos) { x_pos_h = x_pos; }

    /**
     * Get a pointer to the y pos host array.
     * @return pointer to first element in y_pos host array.
     */
    float *get_y_pos_host() { return y_pos_h; }

    void set_y_pos_host(float *y_pos) { y_pos_h = y_pos; }

    /**
     * Get a pointer to the z pos host array.
     * @return pointer to the first element in z_pos host array.
     */
    float *get_z_pos_host() { return z_pos_h; }

    void set_z_pos_host(float *z_pos) { z_pos_h = z_pos; }

    /**
     * Get a pointer to the radius host array.
     * @return pointer to the radius host array.
     */
    float *get_radius_host() { return radius_h; }

    void set_radius_host(float *r_h) { radius_h = r_h; }

    void set_collisions_host(Collision *c_ha) { host_collisions_a = c_ha; }

    void sync_collisions_to_host() {
        thrust::copy(collisions.begin(), collisions.end(), host_collisions_a);
    }

    /**
     * Get a pointer to the x pos device array.
     * @return pointer to first element in x_pos device array.
     */
    float *get_x_pos_device() { return thrust::raw_pointer_cast(x_pos_d.data()); }

    /**
     * Get a pointer to the y pos device array.
     * @return pointer to first element in y_pos device array.
     */
    float *get_y_pos_device() { return thrust::raw_pointer_cast(y_pos_d.data()); }

    /**
     * Get a pointer to the z pos device array.
     * @return pointer to the first element in z_pos device array.
     */
    float *get_z_pos_device() { return thrust::raw_pointer_cast(z_pos_d.data()); }

    /**
     * Get a pointer to the radius device array.
     * @return pointer to the radius device array.
     */
    float *get_radius_device() { return thrust::raw_pointer_cast(radius_d.data()); }

    Collision *get_collisions_device_ptr() { return thrust::raw_pointer_cast(collisions.data()); }

    /**
     * Inform the CollisionSystem that the X, Y, and Z host positions have changed.
     */
    void update_all_from_host();

    /**
     * Inform the CollisionSystem that the X host positions have changed.
     */
    void update_x_pos_from_host();

    /**
     * Inform the CollisionSystem that the Y host positions have changed.
     */
    void update_y_pos_from_host();

    /**
     * Inform the CollisionSystem that the Z host positions have changed.
     */
    void update_z_pos_from_host();

    /**
     * Inform the CollisionSystem that the radiuses have changed.
     */
    void update_radius_from_host();

    /**
     * Runs a first pass through building the BVH tree.
     *
     */
    void init();

    /**
     * Update the x ranks of each object.
     */
    __host__ __device__
    void update_x_pos_ranks();

    /**
     * Update the y ranks of each object.
     */
    __host__ __device__
    void update_y_pos_ranks();

    /**
     * Update the z ranks of each object.
     */
    __host__ __device__
    void update_z_pos_ranks();

    /**
     * Update the morton numbers of each object and rebuild the BVH tree.
     */
    __host__ __device__
    void update_mortons();

    /**
     * Uses the Karras method for updating morton numbers and rebuilds the BVH tree.
     *
     * @param xlims float2 of (min x, max x) to map positions to unsigned ints for morton numbers.
     * @param ylims float2 of (min y, max y) to map positions to unsigned ints for morton numbers.
     * @param zlims float2 of (min z, max z) to map positions to unsigned ints for morton numbers.
     */
    __host__ __device__
    void update_mortons_fast(float2 xlims, float2 ylims, float2 zlims);

    /**
     * Build or rebuild the BVH tree based on current morton numbers.
     */
    __host__ __device__
    void build_tree();

    /**
     * Update the bounding box info for each node in the BVH tree.
     */
    __host__ __device__
    void update_bounding_boxes();

    /**
     * compute a list of collisions that have occured.
     * @return number of collisions found. returns -1 if we did not have enough memory to save all found collisions.
     */
    __host__
    int find_collisions();

    /**
     * to be called from device code: computes a list of collisions that have occured.
     *
     * @return number of collisions found. returns -1 if we did not have enough memory to save all collisions.
     */
    __device__
    int find_collisions_device();

    /**
     * Brute force N^2 algorithm to find collisions. This is faster for small values of N.
     * @return number of collisions that were found.
     */
     __host__
    int find_collisions_N2();

    __device__
    int find_collisions_N2_device();

    size_t N, RESERVATION_SIZE, MAX_COLLISIONS_PER_OBJECT;

    int * num_collisions_d_ptr;
    thrust::device_vector<int> num_collisions_d;

    thrust::counting_iterator<unsigned int> start;
    thrust::counting_iterator<unsigned int> end;

    // buffer for sorting with CUB
    size_t cub_sort_bytes_size;
    void *cub_sort_bytes_ptr;

    // host positions
    float *x_pos_h, *y_pos_h, *z_pos_h, *radius_h;
    bool allocatedByUs;

    // device positions
    float *x_pos_d_ptr, *y_pos_d_ptr, *z_pos_d_ptr, *tmp_pos_d_ptr, *radius_d_ptr;

    thrust::device_vector<float> x_pos_d;
    thrust::device_vector<float> y_pos_d;
    thrust::device_vector<float> z_pos_d;
    thrust::device_vector<float> tmp_pos_d;

    // radiuses of each voxel
    thrust::device_vector<float> radius_d;

    //rank and id vectors.
    unsigned int *x_rank_d_ptr, *y_rank_d_ptr, *z_rank_d_ptr, *tmp_id_a_d_ptr, *tmp_id_b_d_ptr;

    thrust::device_vector<unsigned int> x_rank_d;
    thrust::device_vector<unsigned int> y_rank_d;
    thrust::device_vector<unsigned int> z_rank_d;

    thrust::device_vector<unsigned int> tmp_id_a_d;
    thrust::device_vector<unsigned int> tmp_id_b_d;

    // Morton numbers vector.
    unsigned long long int *mortons_d_ptr, *mortons_tmp_d_ptr;
    unsigned int *mortons_id_d_ptr;

    thrust::device_vector<unsigned long long int> mortons_d;
    thrust::device_vector<unsigned long long int> mortons_tmp_d;
    thrust::device_vector<unsigned int> mortons_id_d;


    // vectors for the BVH tree
    // for the leaf nodes & internal nodes
    unsigned int *leaf_parent_d_ptr, *internal_parent_d_ptr, *internal_childA_d_ptr, *internal_childB_d_ptr;
    thrust::device_vector<unsigned int> leaf_parent_d;
    thrust::device_vector<unsigned int> internal_parent_d;
    thrust::device_vector<unsigned int> internal_childA_d;
    thrust::device_vector<unsigned int> internal_childB_d;

    // for the bounding boxes for all leaf and internal nodes.
    unsigned int *internal_node_bbox_complete_flag_d_ptr;
    BoundingBox *bounding_boxes_d_ptr;
    thrust::device_vector<unsigned int> internal_node_bbox_complete_flag;
    thrust::device_vector<BoundingBox> bounding_boxes_d;

    unsigned int *potential_collisions_idx_d_ptr;
    Collision *potential_collisions_d_ptr, *collisions_d_ptr;
    thrust::device_vector<unsigned int> potential_collisions_idx; // keep track of the idx for the collision we are at.
    thrust::device_vector<Collision> potential_collisions;
    thrust::device_vector<Collision> collisions;
    Collision* host_collisions_a;
    int* host_collisions_b;

    init_morton_func compute_morton_numbers;
    build_bvh_tree_func build_bvh_tree;
    fill_bvh_tree_with_bounding_boxes_func compute_bounding_boxes;
    find_potential_collisions_func find_potential_collisions;
    check_potential_collisions_func check_potential_collisions;
};

__global__ void find_potential_collisions_kernel(int startIdx, int num, find_potential_collisions_func functor);
__global__ void build_tree_kernel(int startIdx, int num, build_bvh_tree_func functor);
#endif //COLLISIONDETECTION_COLLISIONSYSTEM_CUH
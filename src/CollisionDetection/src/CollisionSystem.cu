//
// Created by David Matthews on 5/21/20.
//


#include "../include/CollisionSystem.cuh"
#include "../../Utils/VX3.cuh"

CollisionSystem::CollisionSystem(size_t n, size_t maxCollisionsPerObject, bool toAllocHostVecs) :
        N(n), MAX_COLLISIONS_PER_OBJECT(maxCollisionsPerObject) {
    RESERVATION_SIZE = N;
    start = thrust::counting_iterator<unsigned int>(0);
    end = start + N;

    // allocate memory for positions
    // host
    allocatedByUs = toAllocHostVecs;
    if (allocatedByUs) {
        x_pos_h = new float[RESERVATION_SIZE];
        y_pos_h = new float[RESERVATION_SIZE];
        z_pos_h = new float[RESERVATION_SIZE];
        radius_h = new float[RESERVATION_SIZE];
        host_collisions_a = new Collision[RESERVATION_SIZE * MAX_COLLISIONS_PER_OBJECT];
        host_collisions_b = nullptr;
    }

    num_collisions_d = thrust::device_vector<int>(1);

    // CUB sort buffer
    cub_sort_bytes_size = sizeof(unsigned long long int) * 4 * RESERVATION_SIZE; // Enough memory to sort the morton numbers.
    cudaMalloc(&cub_sort_bytes_ptr, cub_sort_bytes_size);

    // device
    x_pos_d = thrust::device_vector<float>(N);
    y_pos_d = thrust::device_vector<float>(N);
    z_pos_d = thrust::device_vector<float>(N);
    radius_d = thrust::device_vector<float>(N);
    tmp_pos_d = thrust::device_vector<float>(N);

    // alloc rank and id vectors.
    x_rank_d = thrust::device_vector<unsigned int>(N);
    y_rank_d = thrust::device_vector<unsigned int>(N);
    z_rank_d = thrust::device_vector<unsigned int>(N);
    tmp_id_a_d = thrust::device_vector<unsigned int>(N);
    tmp_id_b_d = thrust::device_vector<unsigned int>(N);

    // allocate Morton number vectors
    mortons_d = thrust::device_vector<unsigned long long int>(N);
    mortons_tmp_d = thrust::device_vector<unsigned long long int>(N);
    mortons_id_d = thrust::device_vector<unsigned int>(N);

    // alloc vectors for the BVH tree
    // for the leaf nodes
    leaf_parent_d = thrust::device_vector<unsigned int>(N);


    // for the internal nodes
    internal_parent_d = thrust::device_vector<unsigned int>(N - 1);
    internal_childA_d = thrust::device_vector<unsigned int>(N - 1);
    internal_childB_d = thrust::device_vector<unsigned int>(N - 1);

    // for the bounding boxes for all leaf and internal nodes.
    internal_node_bbox_complete_flag = thrust::device_vector<unsigned int>(N - 1);
    bounding_boxes_d = thrust::device_vector<BoundingBox>(2 * N - 1);
    potential_collisions_idx = thrust::device_vector<unsigned int>(1);

    if (MAX_COLLISIONS_PER_OBJECT * N == 0) {
        throw "Size of potential_collisions array must be > 0";
    }
    potential_collisions = thrust::device_vector<Collision>(MAX_COLLISIONS_PER_OBJECT * N);
    collisions = thrust::device_vector<Collision>(MAX_COLLISIONS_PER_OBJECT * N);

    update_device_pointers_and_functors();

    // init flags to zero.
    thrust::fill(thrust::device, internal_node_bbox_complete_flag_d_ptr, internal_node_bbox_complete_flag_d_ptr + N - 1, 0);

}

CollisionSystem::~CollisionSystem() {
    if (allocatedByUs) {
        delete[] x_pos_h;
        delete[] y_pos_h;
        delete[] z_pos_h;
        delete[] radius_h;
        delete[] host_collisions_a;
    }
    cudaFree(cub_sort_bytes_ptr);

}

void CollisionSystem::set_num_masses(size_t n) {
    bool resizeVectors = n > N;
    set_reserved_num_masses(n);
    N = n;
    end = start + N;

    // don't need to resize the host vectors. Just allocate memory, and assume that it is filled (by you) correctly.
    if (resizeVectors) {
        x_pos_d.resize(N);
        y_pos_d.resize(N);
        z_pos_d.resize(N);
        radius_d.resize(N);

        tmp_pos_d.resize(N);

        x_rank_d.resize(N);
        y_rank_d.resize(N);
        z_rank_d.resize(N);

        tmp_id_a_d.resize(N);
        tmp_id_b_d.resize(N);

        mortons_d.resize(N);
        mortons_tmp_d.resize(N);
        mortons_id_d.resize(N);

        leaf_parent_d.resize(N);

        internal_parent_d.resize(N);
        internal_childA_d.resize(N);
        internal_childB_d.resize(N);

        internal_node_bbox_complete_flag.resize(N, 0);

        bounding_boxes_d.resize(2 * N - 1);

        if (MAX_COLLISIONS_PER_OBJECT * N == 0) {
            throw "Size of potential_collisions array must be > 0";
        }

        potential_collisions.resize(MAX_COLLISIONS_PER_OBJECT * N);
        collisions.resize(MAX_COLLISIONS_PER_OBJECT * N);
    }

}

void CollisionSystem::set_reserved_num_masses(size_t n) {
    if (n < RESERVATION_SIZE) {
        return;
    }

    RESERVATION_SIZE = n;
    // allocate memory for positions
    // host

    // if the memory is allocated by us, re-allocate and copy.
    if (allocatedByUs) {
        float *tmp_x, *tmp_y, *tmp_z, *tmp_r;

        tmp_x = new float[RESERVATION_SIZE];
        thrust::copy(x_pos_h, x_pos_h + N, tmp_x);
        delete[] x_pos_h;
        x_pos_h = tmp_x;

        tmp_y = new float[RESERVATION_SIZE];
        thrust::copy(y_pos_h, y_pos_h + N, tmp_y);
        delete[] y_pos_h;
        y_pos_h = tmp_y;

        tmp_z = new float[RESERVATION_SIZE];
        thrust::copy(z_pos_h, z_pos_h + N, tmp_z);
        delete[] z_pos_h;
        z_pos_h = tmp_z;

        tmp_r = new float[RESERVATION_SIZE];
        thrust::copy(radius_h, radius_h + N, tmp_r);
        delete[] radius_h;
        radius_h = tmp_r;

        Collision *tmp_h = new Collision[RESERVATION_SIZE * MAX_COLLISIONS_PER_OBJECT];
        thrust::copy(host_collisions_a, host_collisions_a + N, tmp_h);
        delete[] host_collisions_a;
        host_collisions_a = tmp_h;
    }

    // CUB sort buffer
    cudaFree(cub_sort_bytes_ptr);
    cub_sort_bytes_size = sizeof(unsigned long long int) * 4 * RESERVATION_SIZE; // Enough memory to sort the morton numbers.
    cudaMalloc(&cub_sort_bytes_ptr, cub_sort_bytes_size);


    // device
    x_pos_d.reserve(RESERVATION_SIZE);
    y_pos_d.reserve(RESERVATION_SIZE);
    z_pos_d.reserve(RESERVATION_SIZE);
    radius_d.reserve(RESERVATION_SIZE);

    //tmp device positions used for calcuating ranks.
    tmp_pos_d.reserve(RESERVATION_SIZE);

    // alloc rank and id vectors.
    x_rank_d.reserve(RESERVATION_SIZE);
    y_rank_d.reserve(RESERVATION_SIZE);
    z_rank_d.reserve(RESERVATION_SIZE);

    tmp_id_a_d.reserve(RESERVATION_SIZE);
    tmp_id_b_d.reserve(RESERVATION_SIZE);

    // allocate Morton number vectors
    mortons_d.reserve(RESERVATION_SIZE);
    mortons_tmp_d.reserve(RESERVATION_SIZE);
    mortons_id_d.reserve(RESERVATION_SIZE);

    // alloc vectors for the BVH tree
    // for the leaf nodes
    leaf_parent_d.reserve(RESERVATION_SIZE);


    // for the internal nodes
    internal_parent_d.reserve(RESERVATION_SIZE);
    internal_childA_d.reserve(RESERVATION_SIZE);
    internal_childB_d.reserve(RESERVATION_SIZE);

    // for the bounding boxes for all leaf and internal nodes.
    internal_node_bbox_complete_flag.reserve(RESERVATION_SIZE);

    bounding_boxes_d.reserve(2 * RESERVATION_SIZE - 1);

    if (MAX_COLLISIONS_PER_OBJECT * RESERVATION_SIZE == 0) {
        throw "Size of potential_collisions array must be > 0";
    }

    potential_collisions.reserve(MAX_COLLISIONS_PER_OBJECT * RESERVATION_SIZE);
    collisions.reserve(MAX_COLLISIONS_PER_OBJECT * RESERVATION_SIZE);


    update_device_pointers_and_functors();
}

void CollisionSystem::set_max_num_cols_per_mass(size_t m) {
    MAX_COLLISIONS_PER_OBJECT = m;

    potential_collisions.reserve(MAX_COLLISIONS_PER_OBJECT * RESERVATION_SIZE);
    potential_collisions.resize(MAX_COLLISIONS_PER_OBJECT * N);

    collisions.reserve(MAX_COLLISIONS_PER_OBJECT * RESERVATION_SIZE);
    collisions.resize(MAX_COLLISIONS_PER_OBJECT * N);

    // update the device pointers that were affected.
    potential_collisions_d_ptr = thrust::raw_pointer_cast(potential_collisions.data());
    collisions_d_ptr = thrust::raw_pointer_cast(collisions.data());
}

void CollisionSystem::update_device_pointers_and_functors() {
    num_collisions_d_ptr = thrust::raw_pointer_cast(num_collisions_d.data());
    x_pos_d_ptr = thrust::raw_pointer_cast(x_pos_d.data());
    y_pos_d_ptr = thrust::raw_pointer_cast(y_pos_d.data());
    z_pos_d_ptr = thrust::raw_pointer_cast(z_pos_d.data());
    radius_d_ptr = thrust::raw_pointer_cast(radius_d.data());
    tmp_pos_d_ptr = thrust::raw_pointer_cast(tmp_pos_d.data());
    x_rank_d_ptr = thrust::raw_pointer_cast(x_rank_d.data());
    y_rank_d_ptr = thrust::raw_pointer_cast(y_rank_d.data());
    z_rank_d_ptr = thrust::raw_pointer_cast(z_rank_d.data());
    tmp_id_a_d_ptr = thrust::raw_pointer_cast(tmp_id_a_d.data());
    tmp_id_b_d_ptr = thrust::raw_pointer_cast(tmp_id_b_d.data());
    mortons_d_ptr = thrust::raw_pointer_cast(mortons_d.data());
    mortons_tmp_d_ptr = thrust::raw_pointer_cast(mortons_tmp_d.data());
    mortons_id_d_ptr = thrust::raw_pointer_cast(mortons_id_d.data());
    leaf_parent_d_ptr = thrust::raw_pointer_cast(leaf_parent_d.data());
    internal_parent_d_ptr = thrust::raw_pointer_cast(internal_parent_d.data());
    internal_childA_d_ptr = thrust::raw_pointer_cast(internal_childA_d.data());
    internal_childB_d_ptr = thrust::raw_pointer_cast(internal_childB_d.data());
    internal_node_bbox_complete_flag_d_ptr = thrust::raw_pointer_cast(internal_node_bbox_complete_flag.data());
    bounding_boxes_d_ptr = thrust::raw_pointer_cast(bounding_boxes_d.data());
    potential_collisions_idx_d_ptr = thrust::raw_pointer_cast(potential_collisions_idx.data());
    potential_collisions_d_ptr = thrust::raw_pointer_cast(potential_collisions.data());
    collisions_d_ptr = thrust::raw_pointer_cast(collisions.data());


    compute_morton_numbers = init_morton_func(x_rank_d_ptr,
            y_rank_d_ptr,
            z_rank_d_ptr,
            mortons_d_ptr,
            mortons_id_d_ptr);

    build_bvh_tree = build_bvh_tree_func(N,
            mortons_d_ptr,
            leaf_parent_d_ptr,
            internal_parent_d_ptr,
            internal_childA_d_ptr,
            internal_childB_d_ptr);

    compute_bounding_boxes = fill_bvh_tree_with_bounding_boxes_func(N,
            bounding_boxes_d_ptr,
            x_pos_d_ptr,
            y_pos_d_ptr,
            z_pos_d_ptr,
            radius_d_ptr,
            mortons_id_d_ptr,
            leaf_parent_d_ptr,
            internal_parent_d_ptr,
            internal_childA_d_ptr,
            internal_childB_d_ptr,
            internal_node_bbox_complete_flag_d_ptr);

    find_potential_collisions = find_potential_collisions_func(N,
            N * MAX_COLLISIONS_PER_OBJECT,
            mortons_id_d_ptr, bounding_boxes_d_ptr, internal_childA_d_ptr, internal_childB_d_ptr, potential_collisions_idx_d_ptr,
            potential_collisions_d_ptr, x_pos_d_ptr, y_pos_d_ptr, z_pos_d_ptr, radius_d_ptr);

    check_potential_collisions = check_potential_collisions_func(x_pos_d_ptr, y_pos_d_ptr, z_pos_d_ptr, radius_d_ptr);
}

void CollisionSystem::update_all_from_host() {
    update_x_pos_from_host();
    update_y_pos_from_host();
    update_z_pos_from_host();
    update_radius_from_host();
}

void CollisionSystem::update_x_pos_from_host() {
    thrust::copy(x_pos_h, x_pos_h + N, x_pos_d.begin());
}

void CollisionSystem::update_y_pos_from_host() {
    thrust::copy(y_pos_h, y_pos_h + N, y_pos_d.begin());
}

void CollisionSystem::update_z_pos_from_host() {
    thrust::copy(z_pos_h, z_pos_h + N, z_pos_d.begin());
}

void CollisionSystem::update_radius_from_host() {
    thrust::copy(radius_h, radius_h + N, radius_d.begin());
}

void CollisionSystem::init() {
    // copy from host to device
    update_all_from_host();

    // compute ranks
    update_x_pos_ranks();
    update_y_pos_ranks();
    update_z_pos_ranks();

    // build and sort mortons
    update_mortons();

    // build BVH tree
    build_tree();

    // fill BVH tree with bounding boxes
    update_bounding_boxes();
}

__host__ __device__
void CollisionSystem::update_x_pos_ranks() {
    // keep track of x object ids after sorting.
    thrust::sequence(thrust::device, tmp_id_a_d_ptr, tmp_id_a_d_ptr + N);


    // copy x positions to tmp one that will get mutated
//    thrust::copy(thrust::device, x_pos_d_ptr, x_pos_d_ptr + N , tmp_pos_d_ptr);

    // sort the positions to determine new ranks.
//    thrust::sort_by_key(thrust::device, tmp_pos_d_ptr, tmp_pos_d_ptr + N, tmp_id_a_d_ptr);

    cub::DeviceRadixSort::SortPairs(cub_sort_bytes_ptr, cub_sort_bytes_size, x_pos_d_ptr, tmp_pos_d_ptr, tmp_id_a_d_ptr, tmp_id_b_d_ptr, N);

    // save the new rank information
    thrust::scatter(thrust::device, start, end, tmp_id_b_d_ptr, x_rank_d_ptr);
}

__host__ __device__
void CollisionSystem::update_y_pos_ranks() {
    // keep track of y object ids after sorting
    thrust::sequence(thrust::device, tmp_id_a_d_ptr, tmp_id_a_d_ptr + N);

    // copy y positions to tmp one that will get mutated.
//    thrust::copy(thrust::device, y_pos_d_ptr, y_pos_d_ptr + N , tmp_pos_d_ptr);

    // sort the positions to determine rank.
//    thrust::sort_by_key(thrust::device, tmp_pos_d_ptr, tmp_pos_d_ptr + N, tmp_id_a_d_ptr);
    cub::DeviceRadixSort::SortPairs(cub_sort_bytes_ptr, cub_sort_bytes_size, y_pos_d_ptr, tmp_pos_d_ptr, tmp_id_a_d_ptr, tmp_id_b_d_ptr, N);

    // save the new rank information
    thrust::scatter(thrust::device, start, end, tmp_id_b_d_ptr, y_rank_d_ptr);
}

__host__ __device__
void CollisionSystem::update_z_pos_ranks() {
    // keep track of z object ids after sorting
    thrust::sequence(thrust::device, tmp_id_a_d_ptr, tmp_id_a_d_ptr + N);

    // copy y positions to tmp one that will get mutated.
//    thrust::copy(thrust::device, z_pos_d_ptr, z_pos_d_ptr + N , tmp_pos_d_ptr);

    // sort the positions to determine rank.
//    thrust::sort_by_key(thrust::device, tmp_pos_d_ptr, tmp_pos_d_ptr + N, tmp_id_a_d_ptr);
    cub::DeviceRadixSort::SortPairs(cub_sort_bytes_ptr, cub_sort_bytes_size, z_pos_d_ptr, tmp_pos_d_ptr, tmp_id_a_d_ptr, tmp_id_b_d_ptr, N);

    // save the new rank information
    thrust::scatter(thrust::device, start, end, tmp_id_b_d_ptr, z_rank_d_ptr);
}

__host__ __device__
void CollisionSystem::update_mortons() {
    // keep track of object ids after sorting.
    thrust::sequence(thrust::device, tmp_id_a_d_ptr, tmp_id_a_d_ptr + N);

    // build morton numbers.
    thrust::for_each(thrust::device, start, end, compute_morton_numbers);

    thrust::copy(thrust::device, mortons_d_ptr, mortons_d_ptr + N, mortons_tmp_d_ptr); // copy mortons to tmp array as source for sorting.

    // sort morton numbers
//    thrust::sort_by_key(thrust::device, mortons_d_ptr, mortons_d_ptr + N, mortons_id_d_ptr);

    cub::DeviceRadixSort::SortPairs(cub_sort_bytes_ptr, cub_sort_bytes_size, mortons_tmp_d_ptr, mortons_d_ptr, tmp_id_a_d_ptr, mortons_id_d_ptr, N);

}

__host__ __device__
void CollisionSystem::update_mortons_fast(float2 xlims, float2 ylims, float2 zlims) {
    thrust::sequence(thrust::device, mortons_id_d_ptr, mortons_id_d_ptr + N);

    // build morton numbers using the Karras method.
    // this will be faster if we are not simulating swarms of particles (e.g. if voxels are evenly distributed across range)
    // but slower if we are simulating swarms of particles that encompass large amounts of area and are not evenly distributed
    // e.g. voxels clumping to form new robots would likely be slower with this method.

    thrust::for_each(thrust::device, start,
            end,
            init_morton_func_fast(xlims,
                    ylims,
                    zlims,
                    x_pos_d_ptr, y_pos_d_ptr, z_pos_d_ptr, mortons_d_ptr, mortons_id_d_ptr));

    // sort morton numbers
    thrust::sort_by_key(thrust::device, mortons_d_ptr, mortons_d_ptr + N, mortons_id_d_ptr);
}

__host__ __device__
void CollisionSystem::build_tree() {
    int num_SM, curDeviceId, gridSize, blockSize;
    cudaGetDevice(&curDeviceId);
    cudaDeviceGetAttribute(&num_SM, cudaDevAttrMultiProcessorCount, curDeviceId);
    blockSize = (int)(N-1)/num_SM;
    if (num_SM * blockSize < (N-1)) {
        blockSize += 1;
    }
    if (blockSize > 256) {
        blockSize = 256;
        gridSize = ((int)N + 254)/256; // N - 1 + 255 leaf nodes.
    } else {
        gridSize = num_SM;
    }
    build_tree_kernel<<<gridSize, blockSize>>>(N - 1, N, build_bvh_tree);
    CUDA_CHECK_AFTER_CALL();
    VcudaDeviceSynchronize();

    // thrust::for_each(thrust::device, start, start + (N - 1), build_bvh_tree);
}

__host__ __device__
void CollisionSystem::update_bounding_boxes() {
    thrust::for_each(thrust::device, start, end, compute_bounding_boxes);
}

__host__
int CollisionSystem::find_collisions() {
    potential_collisions_idx[0] = 0;
    thrust::for_each(thrust::device, start + N - 1, start + 2 * N - 1, find_potential_collisions);

    if (potential_collisions_idx[0] > MAX_COLLISIONS_PER_OBJECT * N) {
        num_collisions_d[0] = -1;
        return -1;
    }

    unsigned int colCount = thrust::copy_if(thrust::device, potential_collisions.begin(),
            potential_collisions.begin() + potential_collisions_idx[0],
            collisions.begin(),
            check_potential_collisions) - collisions.begin();

    num_collisions_d[0]= colCount;
    return colCount;
}

__device__
int CollisionSystem::find_collisions_device() {
    potential_collisions_idx_d_ptr[0] = 0;

    int num_SM, curDeviceId, gridSize, blockSize;
    cudaGetDevice(&curDeviceId);
    cudaDeviceGetAttribute(&num_SM, cudaDevAttrMultiProcessorCount, curDeviceId);
    blockSize = (int)N/num_SM;
    if (num_SM * blockSize < N) {
        blockSize += 1;
    }
    if (blockSize > 256) {
        blockSize = 256;
        gridSize = ((int)N + 255)/256;
    } else {
        gridSize = num_SM;
    }
    find_potential_collisions_kernel<<<gridSize, blockSize>>>(N - 1, N, find_potential_collisions);
    CUDA_CHECK_AFTER_CALL();
    VcudaDeviceSynchronize();

    if (potential_collisions_idx_d_ptr[0] > MAX_COLLISIONS_PER_OBJECT * N) {
        num_collisions_d_ptr[0] = -1;
        return -1;
    }
    num_collisions_d_ptr[0] = potential_collisions_idx_d_ptr[0];
    thrust::copy(thrust::device, potential_collisions_d_ptr, potential_collisions_d_ptr + potential_collisions_idx_d_ptr[0], collisions_d_ptr);
    return potential_collisions_idx_d_ptr[0];
}

__host__
int CollisionSystem::find_collisions_N2() {
    auto keys_a_start = thrust::make_transform_iterator(start, thrust::placeholders::_1 / N);
    auto keys_b_start = thrust::make_transform_iterator(start, thrust::placeholders::_1 % N);

    auto keys_zip_start = thrust::make_zip_iterator(thrust::make_tuple(keys_a_start, keys_b_start));
    auto keys_zip_end = thrust::make_zip_iterator(thrust::make_tuple(keys_a_start + N * N, keys_b_start + N * N));

    potential_collisions_idx[0] = 0;
    thrust::for_each(thrust::device, keys_zip_start,
            keys_zip_end,
            check_potential_collisions_N2_func(N * MAX_COLLISIONS_PER_OBJECT,
                    potential_collisions_idx_d_ptr, x_pos_d_ptr, y_pos_d_ptr, z_pos_d_ptr, radius_d_ptr, collisions_d_ptr));

    if (potential_collisions_idx[0] > MAX_COLLISIONS_PER_OBJECT * N) {
        return -1;
    }

    num_collisions_d[0] = potential_collisions_idx[0];
    return potential_collisions_idx[0];
}

__device__
int CollisionSystem::find_collisions_N2_device() {
    auto keys_a_start = thrust::make_transform_iterator(start, thrust::placeholders::_1 / N);
    auto keys_b_start = thrust::make_transform_iterator(start, thrust::placeholders::_1 % N);

    auto keys_zip_start = thrust::make_zip_iterator(thrust::make_tuple(keys_a_start, keys_b_start));
    auto keys_zip_end = thrust::make_zip_iterator(thrust::make_tuple(keys_a_start + N * N, keys_b_start + N * N));

    potential_collisions_idx_d_ptr[0] = 0;
    thrust::for_each(thrust::device, keys_zip_start,
            keys_zip_end,
            check_potential_collisions_N2_func(N * MAX_COLLISIONS_PER_OBJECT,
                    potential_collisions_idx_d_ptr, x_pos_d_ptr, y_pos_d_ptr, z_pos_d_ptr, radius_d_ptr, collisions_d_ptr));

    if (potential_collisions_idx_d_ptr[0] > MAX_COLLISIONS_PER_OBJECT * N) {
        return -1;
    }

    num_collisions_d_ptr[0] = potential_collisions_idx_d_ptr[0];
    return potential_collisions_idx_d_ptr[0];
}

__global__ void find_potential_collisions_kernel(int startIdx, int num, find_potential_collisions_func functor) {
    int tid = threadIdx.x + blockDim.x * blockIdx.x;
    if (tid < num) {
        functor(tid + startIdx);
    }
}

__global__ void build_tree_kernel(int startIdx, int num, build_bvh_tree_func functor) {
    int tid = threadIdx.x + blockDim.x * blockIdx.x;
    if (tid < num) {
        functor(tid + startIdx);
    }
}
//
// Created by David Matthews on 5/14/20.
//
#include <stdio.h>
#include <stdlib.h>

#include <bitset>
#include <iomanip>

#include "../include/Collisions.cuh"


__device__
int2 determineRange(unsigned long long int *sortedMortonCodes, unsigned int numObjects, unsigned int idx) {
    unsigned long long int morton_i = sortedMortonCodes[idx];
    if (idx == 0) { // if we are on the root node
        return make_int2(0, numObjects - 1);
    }
    // internal node: need to pick direction.
    int d = min(max(__clzll(morton_i ^ sortedMortonCodes[idx + 1]) - __clzll(morton_i ^ sortedMortonCodes[idx - 1]),
            -1), 1);
    int d_min = __clzll(morton_i ^ sortedMortonCodes[idx - d]);

    // compute the upper bound of the range. ensure that the range is valid.
    int l_max = 2;
    int testIdx = idx + l_max * d;
    while (0 <= testIdx && testIdx < numObjects && __clzll(morton_i ^ sortedMortonCodes[testIdx]) > d_min) {
        l_max *= 2;
        testIdx = idx + l_max * d;
    }
//    printf("%d %d %d %d\n", idx,  l_max, __clzll(morton_i ^ sortedMortonCodes[idx + l_max/2*d]), d_min);


    // compute the lower bound of the range with binary search.
    int l = 0;
    for (int t = l_max / 2; t >= 1; t /= 2) {
        int newTestIdx = idx + (l + t) * d;
        if (0 <= newTestIdx && newTestIdx < numObjects && __clzll(morton_i ^ sortedMortonCodes[newTestIdx]) > d_min) {
            l += t;
        }
    }
//    printf("%d %d %d %d\n", idx, l, l_max, d);
    if (d == 1) {
        return make_int2(idx, idx + l * d);
    }
    return make_int2(idx + l * d, idx);
//    return make_int2(1,1);
}

__device__
int findSplit(unsigned long long int *sortedMortonCodes, unsigned int first, unsigned int last) {
    // Identical Morton codes => split the range in the middle.

    unsigned long long int firstCode = sortedMortonCodes[first];
    unsigned long long int lastCode = sortedMortonCodes[last];

    // this should never occur based on how we create the morton numbers from ranked position.
    if (firstCode == lastCode) {
        return (first + last) >> 1;
    }

    // Calculate the number of highest bits that are the same
    // for all objects, using the count-leading-zeros intrinsic.

    int commonPrefix = __clzll(firstCode ^ lastCode);

    // Use binary search to find where the next bit differs.
    // Specifically, we are looking for the highest object that
    // shares more than commonPrefix bits with the first one.

    unsigned int split = first; // initial guess
    unsigned int step = last - first;

    do {
        step = (step + 1) >> 1; // exponential decrease
        int newSplit = split + step; // proposed new position

        if (newSplit < last) {
            unsigned long long int splitCode = sortedMortonCodes[newSplit];
            int splitPrefix = __clzll(firstCode ^ splitCode);
            if (splitPrefix > commonPrefix) {
                split = newSplit;
            } // accept proposal
        }
    } while (step > 1);

    return split;
}

__device__
void build_bvh_tree_func::operator()(unsigned int idx) {
    // Find out which range of objects the node corresponds to.
    // (This is where the magic happens!)

    int2 range = determineRange(sortedMortonCodes, N, idx);
    int first = range.x;
    int last = range.y;

    // Determine where to split the range.
    int split = findSplit(sortedMortonCodes, first, last);


    // Select childA and record parent child relationships.
    if (split == first) {
        leaf_parents[split] = idx;
        internal_childA[idx] = split + N - 1; // denote leaf node by having id be greater than number of internal nodes.
    } else {
        internal_parents[split] = idx;
        internal_childA[idx] = split;
    }

    // Select childB and record parent child relationships.
    if (split + 1 == last) {
        leaf_parents[split + 1] = idx;
        internal_childB[idx] = split + N; // denote leaf node by having id be greater than number of internal nodes.
    } else {
        internal_parents[split + 1] = idx;
        internal_childB[idx] = split + 1;
    }
    if (idx == 0) { // track root node.
        internal_parents[0] = 0XFFFFFFFF;
    }
}


__device__
unsigned int fill_bvh_tree_with_bounding_boxes_func::operator()(int leafIdx) {
    unsigned int depth = 0;
//    BoundingBox *box, *boxL, *boxR; // pointers to the bounding boxes that we need to reference.
    volatile BoundingBox *box, *boxL, *boxR; // pointers to the bounding boxes that we need to reference.

    // always start at a leaf node.
    box = const_cast<volatile BoundingBox *> (bboxTree + leafIdx + N - 1);

    unsigned int objIdx = sortedMortonIds[leafIdx];

    float nodeRadius = radius[objIdx];
    box->x_max = xPts[objIdx] + nodeRadius;
    box->x_min = xPts[objIdx] - nodeRadius;
    box->y_max = yPts[objIdx] + nodeRadius;
    box->y_min = yPts[objIdx] - nodeRadius;
    box->z_max = zPts[objIdx] + nodeRadius;
    box->z_min = zPts[objIdx] - nodeRadius;

    // get parent node of the leaf node.
    unsigned int parentNodeId = leaf_parents[leafIdx];

    while (parentNodeId != 0xFFFFFFFF) {
        __threadfence();
        unsigned int flag = atomicInc(&(bbox_complete_flag[parentNodeId]), 10);
        if (flag == 0) { // we were the first to arrive at the parent node.
            return depth;
        }

        bbox_complete_flag[parentNodeId] = 0; // reset complete flag for future calls to this function.
        box = const_cast<volatile BoundingBox *> (bboxTree +parentNodeId);
        boxL = const_cast<volatile BoundingBox *> (bboxTree + internal_childA[parentNodeId]);
        boxR = const_cast<volatile BoundingBox *> (bboxTree + internal_childB[parentNodeId]);

        box->x_max = max(boxL->x_max, boxR->x_max);
        box->x_min = min(boxL->x_min, boxR->x_min);
        box->y_max = max(boxL->y_max, boxR->y_max);
        box->y_min = min(boxL->y_min, boxR->y_min);
        box->z_max = max(boxL->z_max, boxR->z_max);
        box->z_min = min(boxL->z_min, boxR->z_min);

        parentNodeId = internal_parents[parentNodeId];
        depth++;
    }
    return depth;
}

__device__
void find_potential_collisions_func::operator()(unsigned int idx) {
    // int IDX_OF_INTEREST = 205;
    if (potentialCollisionIdx[0] >= MAX_NUM_COLLISIONS) {
        return;
    }

    unsigned int thisObjectId = sortedMortonIds[idx - NUM_INTERNAL];
    unsigned int currCollisionIdx = 0;
    unsigned int stack[128]; // allocate the stack. (> max depth of tree)
    unsigned int *stackPtr = stack;
    *stackPtr++ = 0xFFFFFFFF; // push onto stack.
    unsigned int node = 0; // root of BVH tree.

    // if (idx == IDX_OF_INTEREST) {
    //     printf("\ninternalChildrenA: ");
    //     for (int i = 0; i < NUM_INTERNAL; i++) {
    //         printf("%u ", internalChildrenA[i]);
    //     }
    //     printf("\n\n");   
    //     printf("internalChildrenB: ");
    //     for (int i = 0; i < NUM_INTERNAL; i++) {
    //         printf("%u ", internalChildrenB[i]);
    //     }
    //     printf("\n");     
    // }

    do {
        unsigned int childL = internalChildrenA[node];
        unsigned int childR = internalChildrenB[node];
        bool overlapL =  overlap(boundingBoxTree + idx, boundingBoxTree + childL);
        bool overlapR =  overlap(boundingBoxTree + idx, boundingBoxTree + childR);

        unsigned int childLR = childL;
        unsigned int childRR = childR;

        while (childLR < NUM_INTERNAL) {
            childLR = internalChildrenB[childLR];
        }
        while (childRR < NUM_INTERNAL) {
            childRR = internalChildrenB[childRR];
        }

        if (childLR <= idx) {
            overlapL = false;
        }
        if (childRR <= idx) {
            overlapR = false;
        }

        // if (idx == IDX_OF_INTEREST) {
        //     printf("\n\nBEGIN: NUM_INTERNAL = %u\n", NUM_INTERNAL);
        //     printf("Idx %d\n", IDX_OF_INTEREST);
        //     printf("childL %u childR %u\n", childL, childR);
        //     printf("childLR %u childRR %u\n", childLR, childRR);
        //     printf("OverlapL: %d and OverlapR: %d\n", int(overlapL), int(overlapR));
        // }

        if (childL != idx && overlapL && childL >= NUM_INTERNAL) {
            // potential collision found. add to list.
            unsigned int childLId = sortedMortonIds[childL - NUM_INTERNAL];
            // if (idx == IDX_OF_INTEREST) {
            //     printf("Collision between %d and %d\n", thisObjectId, childLId);
            // }

            currCollisionIdx = atomicInc(potentialCollisionIdx, 0xFFFFFFFF);
            if (currCollisionIdx >= MAX_NUM_COLLISIONS) {
                return;
            }
            potentialCollisionList[currCollisionIdx].a = min(childLId, thisObjectId); //thisObjectId;
            potentialCollisionList[currCollisionIdx].b = max(childLId, thisObjectId); //childLId;
        }
        if (childR != idx && overlapR && childR >= NUM_INTERNAL) {
            // potential collision found. add to list
            unsigned int childRId = sortedMortonIds[childR - NUM_INTERNAL];
            // if (idx == IDX_OF_INTEREST) {
            //     printf("Collision between %d and %d\n", thisObjectId, childRId);
            // }
            
            currCollisionIdx = atomicInc(potentialCollisionIdx, 0xFFFFFFFF);
            if (currCollisionIdx >= MAX_NUM_COLLISIONS) {
                return;
            }
            potentialCollisionList[currCollisionIdx].a = min(childRId, thisObjectId); //thisObjectId;
            potentialCollisionList[currCollisionIdx].b = max(childRId, thisObjectId); //chilRLId;
        }
        bool traverseL = (overlapL && childL < NUM_INTERNAL);
        bool traverseR = (overlapR && childR < NUM_INTERNAL);

        if (!traverseL && !traverseR) {
            node = *--stackPtr; // pop
            // if (idx == IDX_OF_INTEREST) {
            //     printf("Idx %d popping an element (node: %u) off stack\n", IDX_OF_INTEREST, node);
            // }
        } else {
            node = (traverseL) ? childL : childR;
            // if (idx == IDX_OF_INTEREST) {
            //     printf("Idx %d will traverse to node %u next.\n", IDX_OF_INTEREST, node);
            // }
            if (traverseL && traverseR) {
                *stackPtr++ = childR; // push
                // if (idx == IDX_OF_INTEREST) {
                //     printf("Idx %d pushing an element (node: %u) on stack\n", IDX_OF_INTEREST,  childR);
                // }
            }
        }
    } while (node != 0xFFFFFFFF);
}

__device__
bool check_potential_collisions_func::operator()(const Collision &c) {
    unsigned int colAId = c.a;
    unsigned int colBId = c.b;

    return ((x[colAId] - x[colBId]) * (x[colAId] - x[colBId]) + (y[colAId] - y[colBId]) * (y[colAId] - y[colBId]) +
            (z[colAId] - z[colBId]) * (z[colAId] - z[colBId]) < (r[colAId] + r[colBId]) * (r[colAId] + r[colBId]));
}

__host__ __device__
bool operator<(const Collision &lhs, const Collision &rhs) {
    return (lhs.a < rhs.a || (lhs.a == rhs.a && lhs.b < rhs.b));
}

__host__ __device__
bool operator<=(const Collision &lhs, const Collision &rhs) {
    return (lhs.a <= rhs.a || (lhs.a == rhs.a && lhs.b <= rhs.b));
}

__host__ __device__
bool operator>(const Collision &lhs, const Collision &rhs) {
    return rhs < lhs;
}

__host__ __device__
bool operator>=(const Collision &lhs, const Collision &rhs) {
    return rhs <= lhs;
}

__host__ __device__
bool operator==(const Collision &lhs, const Collision &rhs) {
    return lhs.a == rhs.a && lhs.b == rhs.b;
}

__host__
std::ostream &operator<<(std::ostream &out, const Collision &c) {
    return out << "(" << c.a << ", " << c.b << ")";
}

__host__ __device__
void zip_collisions_func::operator()(unsigned int idx) {
    c[idx].a = min(a[idx], b[idx]);
    c[idx].b = max(a[idx], b[idx]);
}

//__host__ __device__
//void unzip_collisions_func::operator()(unsigned int idx) {
//    a[idx] = c[idx].a;
//    b[idx] = c[idx].b;
//}

unsigned int find_collisions_cpu(unsigned int N, float *x, float *y, float *z, float *r, Collision *cols) {
    int numCols = 0;
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            if ((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]) + (z[i] - z[j]) * (z[i] - z[j]) <
                (r[i] + r[j]) * (r[i] + r[j])) {
                cols[numCols].a = i;
                cols[numCols].b = j;
                numCols++;
            }
        }
    }
    return numCols;
}

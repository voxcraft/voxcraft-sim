//
// Created by David Matthews on 5/15/20.
//

#ifndef COLLISIONDETECTION_UTILS_H
#define COLLISIONDETECTION_UTILS_H

inline unsigned long int next_pow2(unsigned long int x) {
    return x == 1 ? 1 : 1 << (64 - __builtin_clzl(x - 1));
}

inline unsigned int next_pow2(unsigned int x) {
    return x == 1 ? 1 : 1 << (32 - __builtin_clzl(x - 1));
}

#endif //COLLISIONDETECTION_UTILS_H

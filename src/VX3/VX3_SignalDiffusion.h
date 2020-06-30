//
// Created by Sida Liu
//  There are several mechanisms for signal propagation in a system. This class implements a ideal space in which the signal particles diffuse.
//
#if !defined(VX3_SIGNAL_DIFFUSION_H)
#define VX3_SIGNAL_DIFFUSION_H

#include "VX3.cuh"

class VX3_SignalDiffusion
{
public:
    __device__ void doTimestep(); // diffusion
    __device__ void addPartiles(double amount, VX3_Vec3D<> pos); // add signal particles to this environment
    __device__ void removePartiles(double amount, VX3_Vec3D<> pos); // remove signal particles to this environment
};



#endif // VX3_SIGNAL_DIFFUSION_H

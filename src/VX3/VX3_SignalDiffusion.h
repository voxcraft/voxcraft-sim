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
    __device__ void deviceInit(VX3_Vec3D<int> dim=VX3_Vec3D<int>(100,100,100));
    __device__ void doTimestep(double dt); // diffusion
    __device__ double quaryQuantityAtPosition(VX3_Vec3D<int> pos);
    __device__ void addSolute(double amount, VX3_Vec3D<int> pos); // add signal particles to this environment
    __device__ void removeSolute(double amount, VX3_Vec3D<int> pos); // remove signal particles to this environment

    /* data */
    double distance = 1e-2;
    double diffusivity = 2e-9; // ionic in water. refer to: https://www.tandfonline.com/doi/pdf/10.1080/18811248.1996.9732037
    VX3_Vec3D<int> dimension; // make the space into discrete subspace
    double * map, * map_nextstep;
};



#endif // VX3_SIGNAL_DIFFUSION_H

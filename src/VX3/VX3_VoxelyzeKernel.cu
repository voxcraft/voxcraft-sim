#include "VX3_MemoryCleaner.h"
#include "VX3_VoxelyzeKernel.cuh"

/* Tools */
__device__ int bound(int x, int min, int max) {
    if (x < min)
        return min;
    if (x > max)
        return max;
    return x;
}

/* Sub GPU Threads */
__global__ void gpu_update_links(VX3_Link **links, int num);
__global__ void gpu_update_voxels(VX3_Voxel *voxels, int num, double dt, double currentTime, VX3_VoxelyzeKernel *k);
__global__ void gpu_update_temperature(VX3_Voxel *voxels, int num, double TempAmplitude, double TempPeriod, double currentTime, VX3_VoxelyzeKernel* k);
__global__ void gpu_update_attach(VX3_Voxel **surface_voxels, int num, double watchDistance, VX3_VoxelyzeKernel *k);
__global__ void gpu_update_cilia_force(VX3_Voxel **surface_voxels, int num, VX3_VoxelyzeKernel *k);
__global__ void gpu_clear_lookupgrid(VX3_dVector<VX3_Voxel *> *d_collisionLookupGrid, int num);
__global__ void gpu_insert_lookupgrid(VX3_Voxel **d_surface_voxels, int num, VX3_dVector<VX3_Voxel *> *d_collisionLookupGrid,
                                      VX3_Vec3D<> *gridLowerBound, VX3_Vec3D<> *gridDelta, int lookupGrid_n);
__global__ void gpu_collision_attachment_lookupgrid(VX3_dVector<VX3_Voxel *> *d_collisionLookupGrid, int num, double watchDistance,
                                                    VX3_VoxelyzeKernel *k);
__global__ void gpu_update_detach(VX3_Link **links, int num, VX3_VoxelyzeKernel *k);
/* Host methods */

VX3_VoxelyzeKernel::VX3_VoxelyzeKernel(CVX_Sim *In) {

    voxSize = In->Vx.voxSize;

    num_d_voxelMats = In->Vx.voxelMats.size();
    VcudaMalloc((void **)&d_voxelMats, num_d_voxelMats * sizeof(VX3_MaterialVoxel));
    {
        // push all h first, since there will be reference below
        for (auto mat : In->Vx.voxelMats) {
            h_voxelMats.push_back(mat);
        }
        int i = 0;
        for (auto mat : In->Vx.voxelMats) {
            VX3_MaterialVoxel tmp_voxelMat(mat, this);
            VcudaMemcpy(d_voxelMats + i, &tmp_voxelMat, sizeof(VX3_MaterialVoxel), VcudaMemcpyHostToDevice);
            i++;
        }
    }

    num_d_linkMats = In->Vx.linkMats.size();
    VcudaMalloc((void **)&d_linkMats, num_d_linkMats * sizeof(VX3_MaterialLink));
    {
        int i = 0;
        std::vector<VX3_MaterialLink *> tmp_v_linkMats;
        for (CVX_MaterialLink *mat : In->Vx.linkMats) {
            // printf("mat->vox1Mat %p, mat->vox2Mat %p.\n", mat->vox1Mat,
            // mat->vox2Mat);
            VX3_MaterialLink tmp_linkMat(mat, this);
            VcudaMemcpy(d_linkMats + i, &tmp_linkMat, sizeof(VX3_MaterialLink), VcudaMemcpyHostToDevice);
            tmp_v_linkMats.push_back(d_linkMats + i);
            h_linkMats.push_back(mat);
            i++;
        }
        hd_v_linkMats = VX3_hdVector<VX3_MaterialLink *>(tmp_v_linkMats);
    }

    num_d_voxels = In->Vx.voxelsList.size();
    VcudaMalloc((void **)&d_voxels, num_d_voxels * sizeof(VX3_Voxel));
    for (int i = 0; i < num_d_voxels; i++) {
        h_voxels.push_back(In->Vx.voxelsList[i]);
        h_lookup_voxels[In->Vx.voxelsList[i]] = d_voxels + i;
    }
    VcudaMalloc((void **)&d_initialPosition, num_d_voxels * sizeof(Vec3D<>));

    num_d_links = In->Vx.linksList.size();
    std::vector<VX3_Link *> tmp_v_links;
    VcudaMalloc((void **)&d_links, num_d_links * sizeof(VX3_Link));
    VX3_Link *tmp_link_cache = (VX3_Link *)malloc(num_d_links * sizeof(VX3_Link));
    for (int i = 0; i < num_d_links; i++) {
        VX3_Link tmp_link(In->Vx.linksList[i], this);
        memcpy(tmp_link_cache + i, &tmp_link, sizeof(VX3_Link));
        tmp_v_links.push_back(d_links + i); // not copied yet, but still ok to get the address
        h_links.push_back(In->Vx.linksList[i]);
    }
    VcudaMemcpy(d_links, tmp_link_cache, num_d_links * sizeof(VX3_Link), VcudaMemcpyHostToDevice);
    hd_v_links = VX3_hdVector<VX3_Link *>(tmp_v_links);
    for (int i = 0; i < num_d_links; i++) {
        h_lookup_links[In->Vx.linksList[i]] = d_links + i;
    }

    for (int i = 0; i < num_d_voxels; i++) {
        // set values for GPU memory space
        VX3_Voxel tmp_voxel(In->Vx.voxelsList[i], this);
        VcudaMemcpy(d_voxels + i, &tmp_voxel, sizeof(VX3_Voxel), VcudaMemcpyHostToDevice);
    }

    // Not all data is in Vx, here are others:
    DtFrac = In->DtFrac;
    StopConditionType = In->StopConditionType;
    StopConditionValue = In->StopConditionValue;
    TempEnabled = In->pEnv->TempEnabled;
    VaryTempEnabled = In->pEnv->VaryTempEnabled;
    TempBase = In->pEnv->TempBase;
    TempAmplitude = In->pEnv->TempAmplitude;
    TempPeriod = In->pEnv->TempPeriod;
    // currentTemperature = TempBase + TempAmplitude;

    d_surface_voxels = NULL;
}

void VX3_VoxelyzeKernel::cleanup() {
    // The reason not use ~VX3_VoxelyzeKernel is that will be automatically call
    // multiple times after we use memcpy to clone objects.
    MycudaFree(d_linkMats);
    MycudaFree(d_voxels);
    MycudaFree(d_links);
    // MycudaFree(d_collisionsStale);
    if (d_surface_voxels) {
        MycudaFree(d_surface_voxels); // can __device__ malloc pointer be freed
                                      // by cudaFree in __host__??
    }
    // MycudaFree(d_collisions);
}

/* Cuda methods : cannot use any CVX_xxx, and no std::, no boost::, and no
 * filesystem. */

__device__ void VX3_VoxelyzeKernel::syncVectors() {
    d_v_linkMats.clear();
    d_v_collisions.clear();
    d_targets.clear();
    // allocate memory for collision lookup table
    num_lookupGrids = lookupGrid_n * lookupGrid_n * lookupGrid_n;
    d_collisionLookupGrid = (VX3_dVector<VX3_Voxel *> *)malloc(num_lookupGrids * sizeof(VX3_dVector<VX3_Voxel *>));
    if (d_collisionLookupGrid == NULL) {
        printf(COLORCODE_BOLD_RED "ERROR: not enough memory.\n");
    }
    for (int i = 0; i < hd_v_linkMats.size(); i++) {
        d_v_linkMats.push_back(hd_v_linkMats[i]);
    }

    d_v_links.clear();
    for (int i = 0; i < hd_v_links.size(); i++) {
        d_v_links.push_back(hd_v_links[i]);
    }

    for (int i = 0; i < num_d_voxelMats; i++) {
        d_voxelMats[i].syncVectors();
    }

    for (int i = 0; i < num_d_linkMats; i++) {
        d_linkMats[i].syncVectors();
    }

    for (int i = 0; i < num_d_voxels; i++) {
        d_voxels[i].syncVectors();
    }
}
__device__ void VX3_VoxelyzeKernel::saveInitialPosition() {
    for (int i = 0; i < num_d_voxels; i++) {
        d_initialPosition[i] = d_voxels[i].pos;
        // Save this value to voxel, so it can be read out when collecting results in cpu.
        d_voxels[i].isMeasured = (bool) d_voxels[i].mat->isMeasured;
    }
}
__device__ bool VX3_VoxelyzeKernel::StopConditionMet(void) // have we met the stop condition yet?
{
    if (VX3_MathTree::eval(currentCenterOfMass.x, currentCenterOfMass.y, currentCenterOfMass.z, collisionCount, currentTime, recentAngle,
                           targetCloseness, numClosePairs, num_d_voxels, StopConditionFormula) > 0) {
        // double a =
        //     VX3_MathTree::eval(currentCenterOfMass.x, currentCenterOfMass.y, currentCenterOfMass.z, collisionCount, currentTime,
        //     StopConditionFormula);
        // printf("stop score: %f.\n\n", a);
        return true;
    }
    if (forceExit)
        return true;
    return false;
    // if (StopConditionType != SC_MAX_SIM_TIME) {
    //     printf(COLORCODE_BOLD_RED "StopConditionType: %d. Type of stop condition no supported for "
    //                               "now.\n" COLORCODE_RESET,
    //            StopConditionType);
    //     return true;
    // }
    // return currentTime > StopConditionValue ? true : false;
}

__device__ double VX3_VoxelyzeKernel::recommendedTimeStep() {
    // find the largest natural frequency (sqrt(k/m)) that anything in the
    // simulation will experience, then multiply by 2*pi and invert to get the
    // optimally largest timestep that should retain stability
    double MaxFreq2 = 0.0f; // maximum frequency in the simulation in rad/sec
    if (!num_d_links) {
        printf("WARNING: No links.\n");
    }
    if (!num_d_voxels) {
        printf(COLORCODE_BOLD_RED "ERROR: No voxels.\n");
    }
    for (int i = 0; i < num_d_links; i++) {
        VX3_Link *pL = d_links + i;
        // axial
        double m1 = pL->pVNeg->mat->mass(), m2 = pL->pVPos->mat->mass();
        double thisMaxFreq2 = pL->axialStiffness() / (m1 < m2 ? m1 : m2);
        if (thisMaxFreq2 > MaxFreq2)
            MaxFreq2 = thisMaxFreq2;
        // rotational will always be less than or equal
    }
    if (MaxFreq2 <= 0.0f) {                      // didn't find anything (i.e no links) check for
                                                 // individual voxelss
        for (int i = 0; i < num_d_voxels; i++) { // for each link
            double thisMaxFreq2 = d_voxels[i].mat->youngsModulus() * d_voxels[i].mat->nomSize / d_voxels[i].mat->mass();
            if (thisMaxFreq2 > MaxFreq2)
                MaxFreq2 = thisMaxFreq2;
        }
    }
    if (MaxFreq2 <= 0.0f)
        return 0.0f;
    else
        return 1.0f / (6.283185f * sqrt(MaxFreq2)); // the optimal timestep is to advance one
                                                    // radian of the highest natural frequency
}

__device__ void VX3_VoxelyzeKernel::updateTemperature() {
    // updates the temperatures For Actuation!
    // different temperatures in different objs are not support for now.
    if (VaryTempEnabled) {
        if (TempPeriod > 0) {
            int blockSize = 512;
            int minGridSize;
            cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, gpu_update_temperature, 0,
                                               num_d_voxels); // Dynamically calculate blockSize
            int gridSize_voxels = (num_d_voxels + blockSize - 1) / blockSize;
            int blockSize_voxels = num_d_voxels < blockSize ? num_d_voxels : blockSize;
            gpu_update_temperature<<<gridSize_voxels, blockSize_voxels>>>(d_voxels, num_d_voxels, TempAmplitude, TempPeriod, currentTime, this);
            CUDA_CHECK_AFTER_CALL();
            VcudaDeviceSynchronize();
        }
    }
}

__device__ bool VX3_VoxelyzeKernel::doTimeStep(float dt) {
    // clock_t time_measures[10];
    // time_measures[0] = clock();
    updateTemperature();
    CurStepCount++;
    if (dt == 0)
        return true;
    else if (dt < 0) {
        if (!OptimalDt) {
            OptimalDt = recommendedTimeStep();
        }
        if (OptimalDt < 1e-10) {
            CUDA_DEBUG_LINE("recommendedTimeStep is zero.");
            OptimalDt = 1e-10;
            // return false;
        }
        dt = DtFrac * OptimalDt;
    }
    bool Diverged = false;

    int blockSize;
    int minGridSize;
    if (d_v_links.size()) {
        cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, gpu_update_links, 0,
                                           d_v_links.size()); // Dynamically calculate blockSize
        int gridSize_links = (d_v_links.size() + blockSize - 1) / blockSize;
        int blockSize_links = d_v_links.size() < blockSize ? d_v_links.size() : blockSize;
        // if (CurStepCount % 1000 == 0 || currentTime>1.0) {
        //     printf("&d_v_links[0] %p; d_v_links.size() %d. \n", &d_v_links[0], d_v_links.size());
        // }
        gpu_update_links<<<gridSize_links, blockSize_links>>>(&d_v_links[0], d_v_links.size());
        CUDA_CHECK_AFTER_CALL();
        VcudaDeviceSynchronize();

        // checking every link for diverge is too wasteful! using random
        // sampling.
        int r = random(d_v_links.size(), clock());
        if (d_v_links[r]->axialStrain() > 100) {
            CUDA_DEBUG_LINE("Diverged.");
            Diverged = true; // catch divergent condition! (if any thread sets
                             // true we will fail, so don't need mutex...
        }
        if (Diverged)
            return false;
    }

    if (isSurfaceChanged) {
        isSurfaceChanged = false;

        regenerateSurfaceVoxels();
    }

    if (enableAttach || EnableCollision) { // either attachment and collision need measurement for pairwise distances
        updateAttach();
    }
    if (enableDetach) {
        updateDetach();
    }

    if (EnableCilia) {
        cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, gpu_update_cilia_force, 0,
                                           num_d_surface_voxels); // Dynamically calculate blockSize
        int gridSize_voxels = (num_d_surface_voxels + blockSize - 1) / blockSize;
        int blockSize_voxels = num_d_surface_voxels < blockSize ? num_d_surface_voxels : blockSize;
        gpu_update_cilia_force<<<gridSize_voxels, blockSize_voxels>>>(d_surface_voxels, num_d_surface_voxels, this);
        CUDA_CHECK_AFTER_CALL();
        VcudaDeviceSynchronize();
    }

    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, gpu_update_voxels, 0,
                                       num_d_voxels); // Dynamically calculate blockSize
    int gridSize_voxels = (num_d_voxels + blockSize - 1) / blockSize;
    int blockSize_voxels = num_d_voxels < blockSize ? num_d_voxels : blockSize;
    gpu_update_voxels<<<gridSize_voxels, blockSize_voxels>>>(d_voxels, num_d_voxels, dt, currentTime, this);
    CUDA_CHECK_AFTER_CALL();
    VcudaDeviceSynchronize();

    int CycleStep =
        int(TempPeriod / dt); // Sample at the same time point in the cycle, to avoid the impact of actuation as much as possible.
    if (CurStepCount % CycleStep == 0) {
        angleSampleTimes++;

        currentCenterOfMass_history[0] = currentCenterOfMass_history[1];
        currentCenterOfMass_history[1] = currentCenterOfMass;
        updateCurrentCenterOfMass();
        auto A = currentCenterOfMass_history[0];
        auto B = currentCenterOfMass_history[1];
        auto C = currentCenterOfMass;
        if (B == C || A == B || angleSampleTimes < 3) {
            recentAngle = 0; // avoid divide by zero, and don't include first two steps where A and B are still 0.
        } else {
            recentAngle = acos((B - A).Dot(C - B) / (B.Dist(A) * C.Dist(B)));
        }
        // printf("(%d) recentAngle = %f\n", angleSampleTimes, recentAngle);

        // Also calculate targetCloseness here.
        computeTargetCloseness();
    }

    if (SecondaryExperiment) {
        // SecondaryExperiment handle tags:
        // RemoveFromSimulationAfterThisManySeconds
        // ReinitializeInitialPositionAfterThisManySeconds
        // TurnOnThermalExpansionAfterThisManySeconds
        // TurnOnCiliaAfterThisManySeconds

        removeVoxels();
        if (InitialPositionReinitialized == false && ReinitializeInitialPositionAfterThisManySeconds < currentTime) {
            InitialPositionReinitialized = true;
            InitializeCenterOfMass();
            saveInitialPosition();
        }

    }

    currentTime += dt;
    // time_measures[1] = clock();
    // printf("running time for each step: \n");
    // for (int i=0;i<1;i++)
    //     printf("\t%d) %ld clock cycles.\n", i,
    //     time_measures[i+1]-time_measures[i]);
    return true;
}

__device__ void VX3_VoxelyzeKernel::InitializeCenterOfMass() {
    initialCenterOfMass = currentCenterOfMass;
}

__device__ void VX3_VoxelyzeKernel::removeVoxels() {
    for (int i=0;i<num_d_voxelMats;i++) {
        if (d_voxelMats[i].removed == false &&
        d_voxelMats[i].RemoveFromSimulationAfterThisManySeconds > 0 &&
        d_voxelMats[i].RemoveFromSimulationAfterThisManySeconds < currentTime ) {
            VX3_Voxel* neighbor_voxel;

            for (int j=0;j<num_d_voxels;j++) {
                if (d_voxels[j].mat == &d_voxelMats[i] && d_voxels[j].removed == false) {
                    d_voxels[j].removed = true; // mark this voxel as removed
                    for (int k=0;k<6;k++) { // check links in all direction
                        if (d_voxels[j].links[k]) {
                            d_voxels[j].links[k]->removed = true; // mark the link as removed
                            if (d_voxels[j].links[k]->pVNeg == &d_voxels[j]) { // this voxel is pVNeg
                                neighbor_voxel = d_voxels[j].links[k]->pVPos;
                            } else {
                                neighbor_voxel = d_voxels[j].links[k]->pVNeg;
                            }
                            for (int m=0;m<6;m++) {
                                if (neighbor_voxel->links[m] == d_voxels[j].links[k]) {
                                    neighbor_voxel->links[m] = NULL; // delete the neighbor's link
                                    break;
                                }
                            }
                            d_voxels[j].links[k] = NULL; // delete this voxel's link
                        }
                    }
                }
            }
            d_voxelMats[i].removed = true;
            isSurfaceChanged = true;
        }
    }

}

__device__ void VX3_VoxelyzeKernel::updateAttach() {
    // for each surface voxel pair, check distance < watchDistance, make a new
    // link between these two voxels, updateSurface().
    int blockSize;
    int minGridSize;
    if (false) {
        // the parameters of grid are set in gpu_update_voxels, so detection only useful after initialization
        if (gridLowerBound != gridUpperBound) {
            gridDelta = (gridUpperBound - gridLowerBound) / lookupGrid_n;
            if (gridDelta.x < voxSize * 2) {
                gridDelta.x = voxSize * 2;
            }
            if (gridDelta.y < voxSize * 2) {
                gridDelta.y = voxSize * 2;
            }
            if (gridDelta.z < voxSize * 2) {
                gridDelta.z = voxSize * 2;
            }
            // printf("gridLowerBound (%f,%f,%f), gridDelta (%f,%f,%f), gridUpperBound (%f,%f,%f).\n\n", gridLowerBound.x, gridLowerBound.y,
            //        gridLowerBound.z, gridDelta.x, gridDelta.y, gridDelta.z, gridUpperBound.x, gridUpperBound.y, gridUpperBound.z);
            // clear all lookupGrids
            cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, gpu_clear_lookupgrid, 0,
                                               num_lookupGrids); // Dynamically calculate blockSize
            int gridSize_voxels = (num_lookupGrids + blockSize - 1) / blockSize;
            int blockSize_voxels = num_lookupGrids < blockSize ? num_lookupGrids : blockSize;
            gpu_clear_lookupgrid<<<gridSize_voxels, blockSize_voxels>>>(d_collisionLookupGrid, num_lookupGrids);
            CUDA_CHECK_AFTER_CALL();
            VcudaDeviceSynchronize();
            // build lookupGrids: put surface voxels into grids
            cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, gpu_insert_lookupgrid, 0,
                                               num_d_surface_voxels); // Dynamically calculate blockSize
            gridSize_voxels = (num_d_surface_voxels + blockSize - 1) / blockSize;
            blockSize_voxels = num_d_surface_voxels < blockSize ? num_d_surface_voxels : blockSize;
            gpu_insert_lookupgrid<<<gridSize_voxels, blockSize_voxels>>>(d_surface_voxels, num_d_surface_voxels, d_collisionLookupGrid,
                                                                         &gridLowerBound, &gridDelta, lookupGrid_n);
            CUDA_CHECK_AFTER_CALL();
            VcudaDeviceSynchronize();
            // detect collision: voxels in each grid with voxels within this grid and its neighbors
            cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, gpu_collision_attachment_lookupgrid, 0,
                                               num_lookupGrids); // Dynamically calculate blockSize
            gridSize_voxels = (num_lookupGrids + blockSize - 1) / blockSize;
            blockSize_voxels = num_lookupGrids < blockSize ? num_lookupGrids : blockSize;
            gpu_collision_attachment_lookupgrid<<<gridSize_voxels, blockSize_voxels>>>(d_collisionLookupGrid, num_lookupGrids,
                                                                                       watchDistance, this);
            CUDA_CHECK_AFTER_CALL();
            VcudaDeviceSynchronize();
        }
    } else {
        // Pairwise detection O(n ^ 2)
        blockSize = 16;
        dim3 dimBlock(blockSize, blockSize);
        dim3 dimGrid((num_d_surface_voxels + dimBlock.x - 1) / dimBlock.x, (num_d_surface_voxels + dimBlock.y - 1) / dimBlock.y);
        // printf("num_d_surface_voxels %d\n", num_d_surface_voxels);
        gpu_update_attach<<<dimGrid, dimBlock>>>(d_surface_voxels, num_d_surface_voxels, watchDistance,
                                                 this); // invoke two dimensional gpu threads 'CUDA C++ Programming
                                                        // Guide', Nov 2019, P52.
        CUDA_CHECK_AFTER_CALL();
    }
}

__device__ void VX3_VoxelyzeKernel::updateDetach() {
    if (d_v_links.size()) {
        int minGridSize, blockSize;
        cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, gpu_update_detach, 0,
                                           d_v_links.size()); // Dynamically calculate blockSize
        int gridSize_links = (d_v_links.size() + blockSize - 1) / blockSize;
        int blockSize_links = d_v_links.size() < blockSize ? d_v_links.size() : blockSize;
        // if (CurStepCount % 1000 == 0 || currentTime>1.0) {
        //     printf("&d_v_links[0] %p; d_v_links.size() %d. \n", &d_v_links[0], d_v_links.size());
        // }
        gpu_update_detach<<<gridSize_links, blockSize_links>>>(&d_v_links[0], d_v_links.size(), this);
        CUDA_CHECK_AFTER_CALL();
        VcudaDeviceSynchronize();
    }
}

__device__ void VX3_VoxelyzeKernel::updateCurrentCenterOfMass() {
    double TotalMass = 0;
    VX3_Vec3D<> Sum(0, 0, 0);
    for (int i = 0; i < num_d_voxels; i++) {
        if (!d_voxels[i].mat->isMeasured) {
            continue;
        }
        double ThisMass = d_voxels[i].material()->mass();
        Sum += d_voxels[i].position() * ThisMass;
        TotalMass += ThisMass;
    }
    if (TotalMass==0) {
        currentCenterOfMass = VX3_Vec3D<>();
        return;
    }
    currentCenterOfMass = Sum / TotalMass;
}

__device__ void VX3_VoxelyzeKernel::regenerateSurfaceVoxels() {
    // regenerate d_surface_voxels
    if (d_surface_voxels) {
        delete d_surface_voxels;
        d_surface_voxels = NULL;
    }
    VX3_dVector<VX3_Voxel *> tmp;
    for (int i = 0; i < num_d_voxels; i++) {
        d_voxels[i].updateSurface();
        if (d_voxels[i].isSurface() && !d_voxels[i].removed) {
            tmp.push_back(&d_voxels[i]);
        }
    }
    num_d_surface_voxels = tmp.size();
    d_surface_voxels = (VX3_Voxel **)malloc(num_d_surface_voxels * sizeof(VX3_Voxel));
    for (int i = 0; i < num_d_surface_voxels; i++) {
        d_surface_voxels[i] = tmp[i];
    }
}

__device__ VX3_MaterialLink *VX3_VoxelyzeKernel::combinedMaterial(VX3_MaterialVoxel *mat1, VX3_MaterialVoxel *mat2) {
    for (int i = 0; i < d_v_linkMats.size(); i++) {
        VX3_MaterialLink *thisMat = d_v_linkMats[i];
        if ((thisMat->vox1Mat == mat1 && thisMat->vox2Mat == mat2) || (thisMat->vox1Mat == mat2 && thisMat->vox2Mat == mat1))
            return thisMat; // already exist
    }

    VX3_MaterialLink *newMat = new VX3_MaterialLink(mat1, mat2); // where to free this?
    d_v_linkMats.push_back(newMat);
    mat1->d_dependentMaterials.push_back(newMat);
    mat2->d_dependentMaterials.push_back(newMat);

    return newMat;
}

__device__ void VX3_VoxelyzeKernel::computeFitness() {
    VX3_Vec3D<> offset = currentCenterOfMass - initialCenterOfMass;
    fitness_score = VX3_MathTree::eval(offset.x, offset.y, offset.z, collisionCount, currentTime, recentAngle, targetCloseness,
                                       numClosePairs, num_d_voxels, fitness_function);
}

__device__ void VX3_VoxelyzeKernel::registerTargets() {
    for (int i = 0; i < num_d_voxels; i++) {
        auto v = &d_voxels[i];
        if (v->mat->isTarget) {
            d_targets.push_back(v);
        }
    }
}

__device__ void VX3_VoxelyzeKernel::computeTargetCloseness() {
    // this function is called periodically. not very often. once every thousands of steps.
    if (MaxDistInVoxelLengthsToCountAsPair==0)
        return;
    double R = MaxDistInVoxelLengthsToCountAsPair * voxSize;
    double ret = 0;
    numClosePairs = 0;
    for (int i = 0; i < d_targets.size(); i++) {
        for (int j = i + 1; j < d_targets.size(); j++) {
            double distance = d_targets[i]->pos.Dist(d_targets[j]->pos);
            if (distance < R) {
                numClosePairs++;
            }
            ret += 1 / distance;
        }
    }
    targetCloseness = ret;
    // printf("targetCloseness: %f\n", targetCloseness);
}

/* Sub GPU Threads */
__global__ void gpu_update_links(VX3_Link **links, int num) {
    int gindex = threadIdx.x + blockIdx.x * blockDim.x;
    if (gindex < num) {
        VX3_Link *t = links[gindex];
        if (t->removed)
            return;
        if (t->pVPos->mat->fixed && t->pVNeg->mat->fixed)
            return;
        if (t->isDetached)
            return;
        t->updateForces();
        if (t->axialStrain() > 100) {
            printf(COLORCODE_BOLD_RED "ERROR: Diverged.");
        }
    }
}
__global__ void gpu_update_voxels(VX3_Voxel *voxels, int num, double dt, double currentTime, VX3_VoxelyzeKernel *k) {
    int gindex = threadIdx.x + blockIdx.x * blockDim.x;
    if (gindex < num) {
        VX3_Voxel *t = &voxels[gindex];
        if (t->removed)
            return;
        if (t->mat->fixed)
            return; // fixed voxels, no need to update position
        t->timeStep(dt, currentTime, k);

        // update lower bound and upper bound
        if (t->pos.x < k->gridLowerBound.x) {
            k->gridLowerBound.x = t->pos.x;
        } else if (t->pos.x > k->gridUpperBound.x) {
            k->gridUpperBound.x = t->pos.x;
        }
        if (t->pos.y < k->gridLowerBound.y) {
            k->gridLowerBound.y = t->pos.y;
        } else if (t->pos.y > k->gridUpperBound.y) {
            k->gridUpperBound.y = t->pos.y;
        }
        if (t->pos.z < k->gridLowerBound.z) {
            k->gridLowerBound.z = t->pos.z;
        } else if (t->pos.z > k->gridUpperBound.z) {
            k->gridUpperBound.z = t->pos.z;
        }
        // update sticky status
        t->enableAttach = false;
        if (VX3_MathTree::eval(t->pos.x, t->pos.y, t->pos.z, k->collisionCount, currentTime, k->recentAngle, k->targetCloseness,
                               k->numClosePairs, k->num_d_voxels, k->AttachCondition[0]) > 0 &&
            VX3_MathTree::eval(t->pos.x, t->pos.y, t->pos.z, k->collisionCount, currentTime, k->recentAngle, k->targetCloseness,
                               k->numClosePairs, k->num_d_voxels, k->AttachCondition[1]) > 0 &&
            VX3_MathTree::eval(t->pos.x, t->pos.y, t->pos.z, k->collisionCount, currentTime, k->recentAngle, k->targetCloseness,
                               k->numClosePairs, k->num_d_voxels, k->AttachCondition[2]) > 0 &&
            VX3_MathTree::eval(t->pos.x, t->pos.y, t->pos.z, k->collisionCount, currentTime, k->recentAngle, k->targetCloseness,
                               k->numClosePairs, k->num_d_voxels, k->AttachCondition[3]) > 0 &&
            VX3_MathTree::eval(t->pos.x, t->pos.y, t->pos.z, k->collisionCount, currentTime, k->recentAngle, k->targetCloseness,
                               k->numClosePairs, k->num_d_voxels, k->AttachCondition[4]) > 0) {
            t->enableAttach = true;
        };
    }
}

__global__ void gpu_update_temperature(VX3_Voxel *voxels, int num, double TempAmplitude, double TempPeriod, double currentTime, VX3_VoxelyzeKernel* k) {
    int gindex = threadIdx.x + blockIdx.x * blockDim.x;
    if (gindex < num) {
        // vfloat tmp = pEnv->GetTempAmplitude() *
        // sin(2*3.1415926f*(CurTime/pEnv->GetTempPeriod() + pV->phaseOffset)) -
        // pEnv->GetTempBase();
        VX3_Voxel *t = &voxels[gindex];
        if (t->removed)
            return;
        if (t->mat->TurnOnThermalExpansionAfterThisManySeconds > currentTime)
            return;
        if (t->mat->fixed)
            return; // fixed voxels, no need to update temperature
        double currentTemperature =
            TempAmplitude * sin(2 * 3.1415926f * (currentTime / TempPeriod + t->phaseOffset)); // update the global temperature
        // TODO: if we decide not to use PhaseOffset any more, we can move this calculation outside.
        // By default we don't enable expansion. But we can enable that in VXA.
        if (!k->EnableExpansion) {
            if (currentTemperature > 0) {
                currentTemperature = 0;
            }
        }
        t->setTemperature(currentTemperature);
        // t->setTemperature(0.0f);
    }
}
__device__ bool is_neighbor(VX3_Voxel *voxel1, VX3_Voxel *voxel2, VX3_Link *incoming_link, int depth) {
    // printf("Checking (%d,%d,%d) and (%d,%d,%d) in depth %d.\n",
    //             voxel1->ix, voxel1->iy, voxel1->iz,
    //             voxel2->ix, voxel2->iy, voxel2->iz, depth);
    if (voxel1 == voxel2) {
        // printf("found.\n");
        return true;
    }
    if (depth <= 0) { // cannot find in depth
        // printf("not found.\n");
        return false;
    }
    for (int i = 0; i < 6; i++) {
        if (voxel1->links[i]) {
            if (voxel1->links[i] != incoming_link) {
                if (voxel1->links[i]->pVNeg == voxel1) {
                    if (is_neighbor(voxel1->links[i]->pVPos, voxel2, voxel1->links[i], depth - 1)) {
                        return true;
                    }
                } else {
                    if (is_neighbor(voxel1->links[i]->pVNeg, voxel2, voxel1->links[i], depth - 1)) {
                        return true;
                    }
                }
            }
        }
    }
    // printf("not found.\n");
    return false;
}

__device__ void handle_collision_attachment(VX3_Voxel *voxel1, VX3_Voxel *voxel2, double watchDistance, VX3_VoxelyzeKernel *k) {
    // if both of the voxels are fixed, no need to compute.
    if (voxel1->mat->fixed && voxel2->mat->fixed)
        return;

    VX3_Vec3D<double> diff = voxel1->pos - voxel2->pos;
    watchDistance = (voxel1->baseSizeAverage() + voxel2->baseSizeAverage()) * COLLISION_ENVELOPE_RADIUS;

    if (diff.x > watchDistance || diff.x < -watchDistance)
        return;
    if (diff.y > watchDistance || diff.y < -watchDistance)
        return;
    if (diff.z > watchDistance || diff.z < -watchDistance)
        return;

    if (diff.Length() > watchDistance)
        return;

    // to exclude voxels already have link between them. check in depth of
    // 1, direct neighbor ignore the collision
    if (is_neighbor(voxel1, voxel2, NULL, 1)) {
        return;
    }
    // calculate and store contact force, apply and clean in
    // VX3_Voxel::force()
    // if (voxel1->mat !=
    //     voxel2->mat) { // disable same material collision for now
    VX3_Vec3D<> cache_contactForce1, cache_contactForce2;
    if (k->EnableCollision) {
        VX3_Collision collision(voxel1, voxel2);
        collision.updateContactForce();
        cache_contactForce1 = collision.contactForce(voxel1);
        cache_contactForce2 = collision.contactForce(voxel2);
        voxel1->contactForce += cache_contactForce1;
        voxel2->contactForce += cache_contactForce2;
        if ((voxel1->mat->isTarget && !voxel2->mat->isTarget) || (voxel2->mat->isTarget && !voxel1->mat->isTarget)) {
            atomicAdd(&k->collisionCount, 1);
            if (k->EnableSignals) {
                if (voxel1->mat->isTarget) {
                    voxel2->receiveSignal(100, k->currentTime, true);
                } else {
                    voxel1->receiveSignal(100, k->currentTime, true);
                }
            }
        }
    }

    // determined by formula
    if (!voxel1->enableAttach || !voxel2->enableAttach)
        return;

    // fixed voxels, no need to look further for attachment
    if (voxel1->mat->fixed || voxel2->mat->fixed)
        return;
    // different material, no need to attach
    if (voxel1->mat != voxel2->mat)
        return;
    if (!voxel1->mat->sticky)
        return;

    // to exclude voxels already have link between them. check in depth 5.
    // closely connected part ignore the link creation.
    if (is_neighbor(voxel1, voxel2, NULL, 5)) {
        return;
    }

    // determine relative position
    linkDirection link_dir_1, link_dir_2;
    linkAxis link_axis;
    auto a = voxel1->orientation();
    auto b = voxel2->orientation();
    auto c = voxel1->position();
    auto d = voxel2->position();
    auto e = c - d;
    auto ea = a.RotateVec3DInv(-e);
    auto eb = b.RotateVec3DInv(e);

    // first find which is the dominant axis, then determine which one is
    // neg which one is pos.
    VX3_Vec3D<double> f;
    bool reverseOrder = false;
    f = ea.Abs();
    if (f.x >= f.y && f.x >= f.z) { // X_AXIS
        link_axis = X_AXIS;
        if (ea.x < 0) {
            link_dir_1 = X_NEG;
            link_dir_2 = X_POS;
            reverseOrder = true;
        } else {
            link_dir_1 = X_POS;
            link_dir_2 = X_NEG;
        }
    } else if (f.y >= f.x && f.y >= f.z) { // Y_AXIS
        link_axis = Y_AXIS;
        if (ea.y < 0) {
            link_dir_1 = Y_NEG;
            link_dir_2 = Y_POS;
            reverseOrder = true;
        } else {
            link_dir_1 = Y_POS;
            link_dir_2 = Y_NEG;
        }
    } else { // Z_AXIS
        link_axis = Z_AXIS;
        if (ea.z < 0) { // voxel1 is on top
            link_dir_1 = Z_NEG;
            link_dir_2 = Z_POS;
            reverseOrder = true;
        } else {
            link_dir_1 = Z_POS;
            link_dir_2 = Z_NEG;
        }
    }

    // TODO: need to solve this. Create only when there's a right place to
    // attach
    if (voxel1->links[link_dir_1] == NULL && voxel2->links[link_dir_2] == NULL) {
        VX3_Link *pL;
        if (reverseOrder) {
            pL = new VX3_Link(voxel1, link_dir_1, voxel2, link_dir_2, link_axis,
                              k); // make the new link (change to both materials, etc.
        } else {
            pL = new VX3_Link(voxel2, link_dir_2, voxel1, link_dir_1, link_axis,
                              k); // make the new link (change to both materials, etc.
        }
        if (!pL) {
            printf(COLORCODE_BOLD_RED "ERROR: Out of memory. Link not created.\n");
            return;
        }
        pL->isNewLink = k->SafetyGuard;
        k->d_v_links.push_back(pL); // add to the list

        k->isSurfaceChanged = true;

        // printf("createLink.... %p %p distance=> %f %f %f (%f), dir (%d and "
        //        "%d), watchDistance %f.\n",
        //        voxel1, voxel2, diff.x, diff.y, diff.z, diff.Length(),
        //        link_dir_1, link_dir_2, watchDistance);
        // printf("orientation (%f; %f, %f, %f) and (%f; %f, %f, %f).\n", a.w,
        //        a.x, a.y, a.z, b.w, b.x, b.y, b.z);
        // printf("ea, after inv rotate (%f, %f, %f)", ea.x, ea.y, ea.z);
        // printf("newLink: rest %f.\n", pL->currentRestLength);
        // printf("between (%d,%d,%d) and (%d,%d,%d).\n", voxel1->ix,
        //        voxel1->iy, voxel1->iz, voxel2->ix, voxel2->iy, voxel2->iz);

        // if a link is created, set contact force = 0 , for stable reason. (if they are connected, they should not collide.)
        voxel1->contactForce -= cache_contactForce1;
        voxel2->contactForce -= cache_contactForce2;
    }
}

__global__ void gpu_update_attach(VX3_Voxel **surface_voxels, int num, double watchDistance, VX3_VoxelyzeKernel *k) {
    int first = threadIdx.x + blockIdx.x * blockDim.x;
    int second = threadIdx.y + blockIdx.y * blockDim.y;
    if (first < num && second < first) {
        VX3_Voxel *voxel1 = surface_voxels[first];
        VX3_Voxel *voxel2 = surface_voxels[second];
        if (voxel1->removed || voxel2->removed)
            return;
        handle_collision_attachment(voxel1, voxel2, watchDistance, k);
    }
}

// TODO: only need to update after attachment changes.
__global__ void gpu_update_cilia_force(VX3_Voxel **surface_voxels, int num, VX3_VoxelyzeKernel *k) {
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    if (index < num) {
        if (surface_voxels[index]->removed)
            return;
        if (surface_voxels[index]->mat->Cilia == 0)
            return;
        if (surface_voxels[index]->mat->TurnOnCiliaAfterThisManySeconds > k->currentTime)
            return;
        // rotate base cilia force and update it into voxel.
        surface_voxels[index]->CiliaForce = surface_voxels[index]->orient.RotateVec3D(
            surface_voxels[index]->baseCiliaForce + surface_voxels[index]->localSignal * surface_voxels[index]->shiftCiliaForce);
    }
}

__global__ void gpu_clear_lookupgrid(VX3_dVector<VX3_Voxel *> *d_collisionLookupGrid, int num) {
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    if (index < num) {
        d_collisionLookupGrid[index].clear();
    }
}

__global__ void gpu_insert_lookupgrid(VX3_Voxel **d_surface_voxels, int num, VX3_dVector<VX3_Voxel *> *d_collisionLookupGrid,
                                      VX3_Vec3D<> *gridLowerBound, VX3_Vec3D<> *gridDelta, int lookupGrid_n) {
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    if (index < num) {
        VX3_Voxel *v = d_surface_voxels[index];
        int ix = int((v->pos.x - gridLowerBound->x) / gridDelta->x);
        int iy = int((v->pos.y - gridLowerBound->y) / gridDelta->y);
        int iz = int((v->pos.z - gridLowerBound->z) / gridDelta->z);
        bound(ix, 0, lookupGrid_n);
        bound(iy, 0, lookupGrid_n);
        bound(iz, 0, lookupGrid_n);
        d_collisionLookupGrid[ix * lookupGrid_n * lookupGrid_n + iy * lookupGrid_n + iz].push_back(v);
    }
}

__global__ void gpu_pairwise_detection(VX3_Voxel **voxel1, VX3_Voxel **voxel2, int num_v1, int num_v2, double watchDistance,
                                       VX3_VoxelyzeKernel *k) {
    int index_x = threadIdx.x + blockIdx.x * blockDim.x;
    int index_y = threadIdx.y + blockIdx.y * blockDim.y;
    if (index_x < num_v1 && index_y < num_v2) {
        if (voxel1[index_x]->removed || voxel2[index_y]->removed)
            return;
        handle_collision_attachment(voxel1[index_x], voxel2[index_y], watchDistance, k);
    }
}

__device__ int index_3d_to_1d(int x, int y, int z, int dim_len) { return x * dim_len * dim_len + y * dim_len + z; }
__device__ VX3_Vec3D<int> index_1d_to_3d(int n, int dim_len) {
    VX3_Vec3D<int> v;
    v.x = int(floor(double(n / (dim_len * dim_len)))) % dim_len;
    v.y = int(floor(double(n / dim_len))) % dim_len;
    v.z = n % dim_len;
    return v;
}

__global__ void gpu_collision_attachment_lookupgrid(VX3_dVector<VX3_Voxel *> *d_collisionLookupGrid, int num, double watchDistance,
                                                    VX3_VoxelyzeKernel *k) {
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    if (index < num) {
        int num_voxel_in_grid = d_collisionLookupGrid[index].size();
        if (num_voxel_in_grid == 0)
            return;
        // within the grid
        int dim_len = k->lookupGrid_n;
        auto index_3d = index_1d_to_3d(index, dim_len);
        int ix = index_3d.x;
        int iy = index_3d.y;
        int iz = index_3d.z;
        // printf("num_voxel_in_grid %d[%d][%d][%d]: %d\n", index, ix, iy, iz, num_voxel_in_grid);
        int blockSize = 16;
        dim3 dimBlock(blockSize, blockSize);
        dim3 dimGrid((num_voxel_in_grid + dimBlock.x - 1) / dimBlock.x, (num_voxel_in_grid + dimBlock.y - 1) / dimBlock.y);
        gpu_pairwise_detection<<<dimGrid, dimBlock>>>(&d_collisionLookupGrid[index][0], &d_collisionLookupGrid[index][0], num_voxel_in_grid,
                                                      num_voxel_in_grid, watchDistance, k);
        // invoke two dimensional gpu threads 'CUDA C++ Programming
        // Guide', Nov 2019, P52.
        CUDA_CHECK_AFTER_CALL();
        // with neighbors
        for (int dix = -1; dix <= 1; dix++) {
            for (int diy = -1; diy <= 1; diy++) {
                for (int diz = -1; diz <= 1; diz++) {
                    int index_2 = index_3d_to_1d(ix + dix, iy + diy, iz + diz, dim_len);
                    if (index_2 > index && index_2 < num) {
                        int num_voxel_in_grid_2 = d_collisionLookupGrid[index_2].size();
                        if (num_voxel_in_grid_2 > 0) {
                            gpu_pairwise_detection<<<dimGrid, dimBlock>>>(
                                &d_collisionLookupGrid[index][0],
                                &d_collisionLookupGrid[index_3d_to_1d(ix + dix, iy + diy, iz + diz, dim_len)][0], num_voxel_in_grid,
                                num_voxel_in_grid_2, watchDistance, k);
                        }
                    }
                }
            }
        }
        CUDA_CHECK_AFTER_CALL();
    }
}

__global__ void gpu_update_detach(VX3_Link **links, int num, VX3_VoxelyzeKernel* k) {
    int gindex = threadIdx.x + blockIdx.x * blockDim.x;
    if (gindex < num) {
        VX3_Link *t = links[gindex];
        if (t->removed)
            return;
        if (t->isDetached)
            return;
        // clu: vxa: MatModel=1, Fail_Stress=1e+6 => Fail_Stress => failureStress => isFailed.
        if (t->isFailed()) {
            t->isDetached = true;
            for (int i = 0; i < 6; i++) {
                if (t->pVNeg->links[i] == t) {
                    t->pVNeg->links[i] = NULL;
                }
                if (t->pVPos->links[i] == t) {
                    t->pVPos->links[i] = NULL;
                }
            }
            k->isSurfaceChanged = true;
        }
    }
}
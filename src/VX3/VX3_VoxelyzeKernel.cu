#include <stdarg.h>
#include "VX3_Collision.h"
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
__global__ void gpu_update_collision_system_pos_radius(VX3_Voxel **surface_voxels, int num, double watchDistance, VX3_VoxelyzeKernel *k);
__global__ void gpu_update_sync_collisions(VX3_Voxel **surface_voxels, int num, double watchDistance, VX3_VoxelyzeKernel *k);
__global__ void gpu_update_attach(VX3_Voxel **surface_voxels, int num, double watchDistance, VX3_VoxelyzeKernel *k);
__global__ void gpu_update_cilia_force(VX3_Voxel **surface_voxels, int num, VX3_VoxelyzeKernel *k);
__global__ void gpu_update_brownian_motion(VX3_Voxel **surface_voxels, int num, int WorldSize, double seed, double currentTime, VX3_VoxelyzeKernel *k);  // sam
__global__ void gpu_clear_lookupgrid(VX3_dVector<VX3_Voxel *> *d_collisionLookupGrid, int num);
__global__ void gpu_insert_lookupgrid(VX3_Voxel **d_surface_voxels, int num, VX3_dVector<VX3_Voxel *> *d_collisionLookupGrid,
                                      VX3_Vec3D<> *gridLowerBound, VX3_Vec3D<> *gridDelta, int lookupGrid_n);
__global__ void gpu_collision_attachment_lookupgrid(VX3_dVector<VX3_Voxel *> *d_collisionLookupGrid, int num, double watchDistance,
                                                    VX3_VoxelyzeKernel *k);
__global__ void gpu_surface_grow(VX3_Voxel ** surface_voxels, int num);
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
    num_d_init_voxels = num_d_voxels;
    VcudaMalloc((void **)&d_voxels, (num_d_voxels + MaxNewVoxelsAddedMidSim) * sizeof(VX3_Voxel)); // pre-allocate memory for new Voxels
    CUDA_CHECK_AFTER_CALL();

    for (int i = 0; i < num_d_voxels; i++) {
        h_voxels.push_back(In->Vx.voxelsList[i]);
        h_lookup_voxels[In->Vx.voxelsList[i]] = d_voxels + i;
    }
    VcudaMalloc((void **) &d_initialPosition, (num_d_voxels + MaxNewVoxelsAddedMidSim) * sizeof(Vec3D<>));
    CUDA_CHECK_AFTER_CALL();

    // Create the collison system and copy it to the device.
    h_collision_system = new CollisionSystem((num_d_voxels + MaxNewVoxelsAddedMidSim), 128, false);
    VcudaMalloc((void **) &d_collision_system, sizeof(CollisionSystem));
    CUDA_CHECK_AFTER_CALL();
    VcudaMemcpy(d_collision_system, h_collision_system, sizeof(CollisionSystem), cudaMemcpyHostToDevice);

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

__device__ void VX3_VoxelyzeKernel::deviceInit() {
    PRINT(this, "Kernel (%p) deviceInit\n", this);
    d_v_linkMats.clear();
    d_v_collisions.clear();
    d_targets.clear();
    d_voxelgroups.clear();
    d_voxel_to_update_group.clear();
    d_voxels_to_detach.clear();

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
        d_voxelMats[i].deviceInit();
    }

    for (int i = 0; i < num_d_linkMats; i++) {
        d_linkMats[i].deviceInit();
    }

    for (int i = 0; i < num_d_voxels; i++) {
        d_voxels[i].deviceInit(this);
    }
    
    for (int i = 0; i < num_d_links; i++) {
        d_links[i].deviceInit(this);
    }

    d_attach_manager = new VX3_AttachManager(this);

    d_growth_manager = new VX3_GrowthManager(this);

    staticWatchDistance = 2 * COLLISION_ENVELOPE_RADIUS * watchDistance * voxSize;
    staticWatchDistance_square = staticWatchDistance * staticWatchDistance;

    randomGenerator = new RandomGenerator();

    regenerateSurfaceVoxels();

    registerTargets();  // sam
    
    for (int i=0;i<num_d_voxels;i++) {
        d_voxels[i].d_group->updateGroup(&d_voxels[i]); //it'll automatically skip those don't need update.
    }
    
}
__device__ void VX3_VoxelyzeKernel::saveInitialPosition() {
    PRINT(this, "d_initialPosition (%p).\n", d_initialPosition);
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
    if (!SkipThoroughTest) {
        if (!ThoroughValidationCheck())
            return false;
    }

    CurStepCount++;
    updateTemperature();
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
            printf("Diverged.");
            Diverged = true; // catch divergent condition! (if any thread sets
                             // true we will fail, so don't need mutex...
        }
        if (Diverged)
            return false;
    }

    if (enableAttach || EnableCollision) { // either attachment and collision need measurement for pairwise distances
        // updateAttach(0);
        // int num_cols = tmpCollisionCount;
        updateAttach(CollisionMode);
        // if(num_cols != tmpCollisionCount) {
            // printf("ERROR!!\n N2 algorithm found a different number of collisions than the Tree algorithm did!\n N2: %d\nTree: %d\n", num_cols, tmpCollisionCount);
        // }
        // assert(tmpCollisionCount == num_cols);    
    }

    if (EnableCilia) {

        // sam:
        if (RandomizeCiliaEvery > 0) { 
            // int CiliaStep = int(RandomizeCiliaEvery / dt);
            // if (CurStepCount % CiliaStep == 0) {
            if (currentTime >= lastBrownianUpdateTime + RandomizeCiliaEvery) {
                // in parallel:
                cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, gpu_update_brownian_motion, 0, num_d_surface_voxels);
                int gridSize_voxels = (num_d_surface_voxels + blockSize - 1) / blockSize;
                int blockSize_voxels = num_d_surface_voxels < blockSize ? num_d_surface_voxels : blockSize;
                gpu_update_brownian_motion<<<gridSize_voxels, blockSize_voxels>>>(d_surface_voxels, num_d_surface_voxels, WorldSize, RandomSeed, currentTime, this);
                CUDA_CHECK_AFTER_CALL();
                VcudaDeviceSynchronize();

                // // just do a big loop:
                // updateBrownianMotion(RandomSeed, currentTime);

                lastBrownianUpdateTime = currentTime;
            }
        }


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
        computeLargestStickyGroupSize();
    }

    if (EnableSurfaceGrowth) {
        if (currentTime>=SurfaceGrowth_activeTime) {
            SurfaceGrowth_activeTime = currentTime + SurfaceGrowth_Interval;
            PRINT(this, "call surfaceGrow at currentTime: %f. with SurfaceGrowth_Rate %f.\n", currentTime, SurfaceGrowth_Rate);
            surfaceGrow();
            // SurfaceGrowth_Growed = 0;
        }
        // if (SurfaceGrowth_Growed < ceil(num_d_surface_voxels * SurfaceGrowth_Rate)) {
        //     if (d_growth_manager->grow()) {
        //         SurfaceGrowth_Growed++;
        //     };
        // }
    }

    if (isSurfaceChanged) {
        if (currentTime >= lastRegenerateSurfaceTime) {
            lastRegenerateSurfaceTime = currentTime + 0.05; // regenerate at most once per 0.1 second simulation time.
            isSurfaceChanged = false;
            regenerateSurfaceVoxels();
        }
    }

    {
        // only update Groups after all operation is done at each step
        updateGroups();
    }

    // Sida: after update Groups, information will be good to determine who's going to be removed.
    // sam:
    if (SecondaryExperiment) {
        // SecondaryExperiment handle tags:
        // RemoveFromSimulationAfterThisManySeconds
        // ReinitializeInitialPositionAfterThisManySeconds
        // TurnOnThermalExpansionAfterThisManySeconds
        // TurnOnCiliaAfterThisManySeconds

        // removeVoxels();
        // if (InitialPositionReinitialized == false && ReinitializeInitialPositionAfterThisManySeconds < currentTime) {
        //     InitialPositionReinitialized = true;
        //     InitializeCenterOfMass();
        //     saveInitialPosition();

        //     registerTargets();
        // }     

        double nextReplicationTime = lastReplicationTime + ReinitializeInitialPositionAfterThisManySeconds;
        
        if (false) {
            if (currentTime >= nextReplicationTime) {
                addVoxel(0, 0, 20, 1);
                addVoxel(2, 2, 20, 1);
                lastReplicationTime = currentTime;  // reset timer
            }
        }

        if (SelfReplication) {  // kinematic self replication experiments

            // debris treadmill
            if ( (ReplenishDebrisEvery > 0) && (currentTime >= lastReplenishDebrisTime + ReplenishDebrisEvery) ) {
                replenishMaterial(2, WorldSize, SpaceBetweenDebris+1, DebrisMat-1, DebrisHeight);  // replenish sticky building material along the surface plane
                lastReplenishDebrisTime = currentTime;  // reset timer
            }

            // Step 1: Remove small bodies and allow self-attachment with momentum before tranisiton
            if ( (InitialPositionReinitialized) && (currentTime >= nextReplicationTime) ) {
                // InitialPositionReinitialized is true at t=0

                // removeVoxels();  // remove mat 0 voxels that are flagged: removed=true
                // reInitAllGroups();  // EXTREME

                computeTargetCloseness(); // in case no bots remain and sim ends
                computeLargestStickyGroupSize();
                
                convertMatIfSmallBody(1, 0, DebrisHeight+1);   // convert small mat1 bodies>1 to mat0; flags material as not yet removed
                removeVoxels();  // remove bots and newly converted small mat 0 bodies
                InitialPositionReinitialized = false; // false = transition period
            }

            // Transition period
            if (!InitialPositionReinitialized) {
                // do stuff here before next round
                // push debris to ground?
                // vary gravity in a concave function returning to normal
            }

            // Step 2: Transform piles into organisms (mat1->mat0)
            if (currentTime >= nextReplicationTime + SettleTimeBeforeNextRoundOfReplication) {
                if (ReplenishDebrisEvery == 0) {
                    replenishMaterial(2, WorldSize, SpaceBetweenDebris+1, DebrisMat-1, DebrisHeight);  // replenish just before new filial generation
                }
                convertMatIfLargeBody(1, 0); // finally, convert large mat1 bodies>MinimumBotSize to mat0 [new organisms]
                InitializeCenterOfMass();  // in case fitness is a function of x,y,z
                saveInitialPosition();
                InitialPositionReinitialized = true;  // switch that allows settle time between removing voxels and next round
                lastReplicationTime = currentTime;  // reset timer

                // check for inconsistent voxel groups (bad attach/detach)
                if (!ThoroughValidationCheck()) {
                    convertMatIfSmallBody(1, 0, DebrisHeight+1);  // so numClosePairs will be 0
                    removeVoxels();  // failed the test so remove all the bots (causes sim to end)
                }
            }
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



// sam:
__device__ void VX3_VoxelyzeKernel::computeLargestStickyGroupSize() {
    largestStickyGroupSize = 0;
    for (int i = 0; i < d_voxelgroups.size(); i++) {
        
        int thisSize = d_voxelgroups[i]->d_voxels.size();
        bool sticky = d_voxelgroups[i]->d_voxels[0]->mat->sticky;

        if ( (sticky) && (thisSize > largestStickyGroupSize) )
            largestStickyGroupSize = thisSize;

    }
}

// sam:
__device__ void VX3_VoxelyzeKernel::computeNumRealLinks() {
    numRealLinks = 0;
    for (int i = 0; i < d_v_links.size(); i++) {
        auto l = d_v_links[i];
        if (l->removed)
            continue;
        if (l->pVNeg->mat->fixed || l->pVPos->mat->fixed)
            continue;
        if (l->isDetached)
            continue;

        numRealLinks++;
    }
}

// sam:
__device__ bool VX3_VoxelyzeKernel::EarlyStopIfNoBotsRemain() {

    if (!InitialPositionReinitialized){
        return false; // provides some time in between voxel removal and next round
    }
    
    bool earlyStoppage = false;
    for (int i = 0; i < num_d_voxelMats; i++) {
        if(d_voxelMats[i].EndSimIfCompletelyRemoved) {
            earlyStoppage = true;
        }
    }
    if (!earlyStoppage) { 
        return false;
    }

    for (int i=0;i<num_d_voxels;i++) {
        if ( (d_voxels[i].mat->EndSimIfCompletelyRemoved) && (!d_voxels[i].removed) ) {
            return false;
        }
    }
    return true;
}

// sam:
__device__ void VX3_VoxelyzeKernel::replenishMaterial(int start, int end, int step, int mat, int height) {
    // TODO: make this a material attribute, replenish at initialized position.
    for (int x = start; x < end; x+=step) {
        for (int y = start; y < end; y+=step) {
            for (int z = 0; z < height; z++) {
                addVoxel(x, y, z, mat);
            }
        }
    }
}


// sam:
__device__ void VX3_VoxelyzeKernel::updateBrownianMotion(double seed, double currentTime) {
    curandState state;
    curand_init(seed + currentTime, 0, 0, &state);  // seed, sequence number (in a loop this can be 0), offset, state
    for (int i = 0; i < num_d_surface_voxels; i++) {

        if (d_surface_voxels[i]->removed)
            return;
        if (d_surface_voxels[i]->mat->Cilia == 0)
            return;
        if (d_surface_voxels[i]->mat->TurnOnCiliaAfterThisManySeconds > currentTime)
            return;

        d_surface_voxels[i]->baseCiliaForce.x = 2*curand_uniform(&state)-1;
        d_surface_voxels[i]->baseCiliaForce.y = 2*curand_uniform(&state)-1;

    }
}

// sam:
__device__ void VX3_VoxelyzeKernel::reInitAllGroups() {

    d_voxelgroups.clear();

    for (int i = 0; i < num_d_voxels; i++) {
        d_voxels[i].deviceInit(this);
        d_voxels[i].d_group->updateGroup(&d_voxels[i]);
    }
    
    // for (int i=0;i<d_voxelgroups.size();i++) {
    //     d_voxelgroups[i]->updateGroup();
    // }
}

// sam:
__device__ void VX3_VoxelyzeKernel::convertMatIfSmallBody(int mat1, int mat2, int minSizeToConvert) {

    d_voxelMats[mat2].removed = false;  // this material will be removed next removeVoxels() call

    for (int i=0;i<num_d_voxels;i++) {
        if (d_voxels[i].removed) continue;
        if (d_voxels[i].mat == &d_voxelMats[mat1]) {
            if (d_voxels[i].d_group->removed) {
                printf("Sida: A voxel's group should not be removed. d_voxels[%d].\n", i);
            }
            if ( (d_voxels[i].d_group->d_voxels.size() >= minSizeToConvert) && (d_voxels[i].d_group->d_voxels.size() < MinimumBotSize) ) {
                d_voxels[i].mat = &d_voxelMats[mat2];
            }
            // resets debris that moved but didn't attach
            // TODO: debug initial position of newly added voxels
            // if ( (singleton) && (abs(d_initialPosition[i].x - d_voxels[i].pos.x) > voxSize*0.1) ) {
            //         d_voxels[i].mat = &d_voxelMats[mat2]; // moved out of position: convert it 
            // } 
        }
    }
}


// sam:
__device__ void VX3_VoxelyzeKernel::convertMatIfLargeBody(int mat1, int mat2) {

    for (int i=0;i<num_d_voxels;i++) {
        if ( (d_voxels[i].mat == &d_voxelMats[mat1]) && (d_voxels[i].d_group->d_voxels.size() >= MinimumBotSize) ) {
            d_voxels[i].mat = &d_voxelMats[mat2];
        }
    }
    d_voxelMats[mat2].removed = false;
}

// sam:
__device__ bool VX3_VoxelyzeKernel::addVoxel(int x, int y, int z, int mat) {

    float r = float(voxSize);

    bool voxAlreadyThere = d_collision_system->check_collisions_device(float(x)*r, float(y)*r, float(z)*r, r);

    if ( (!voxAlreadyThere) && (num_d_voxels - num_d_init_voxels < MaxNewVoxelsAddedMidSim) ) { // memory limitation, refer to pre-allocation.
        d_voxels[num_d_voxels].deviceInit(this); // do this first
        d_voxels[num_d_voxels].pos = VX3_Vec3D<>(float(x)*r, float(y)*r, float(z)*r);
        d_voxels[num_d_voxels].orient = VX3_Quat3D<>(); // default orientation
        d_voxels[num_d_voxels].mat = &d_voxelMats[mat];
        d_voxels[num_d_voxels].baseCiliaForce = VX3_Vec3D<>(0.0, -1.0, -1.0);
        atomicAdd(&num_d_voxels, 1); // safer to use atomic add.
        isSurfaceChanged = true;
        return true;
    }
    return false;
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
                    d_voxels[j].d_group->removed = true;
                    PRINT(this, "Remove voxel (%p) and group (%p) explicitly.\n", &d_voxels[j], d_voxels[j].d_group);

                    for (int k=0;k<6;k++) { // check links in all direction
                        if (d_voxels[j].links[k]) {

                            neighbor_voxel = d_voxels[j].adjacentVoxel((linkDirection)k);
                            d_voxels[j].links[k]->removed = true; // mark the link as removed
                            neighbor_voxel->links[oppositeDirection(k)] = NULL; // delete the neighbor's link
                            d_voxels[j].links[k] = NULL; // delete this voxel's link
                            PRINT(this, "Remove link (%p).\n", d_voxels[j].links[k])
                        }
                    }
                }
            }
            d_voxelMats[i].removed = true;
            PRINT(this, "Remove voxel material (%p).\n", &d_voxelMats[i]);
            isSurfaceChanged = true;
        }
    }
}

__device__ void VX3_VoxelyzeKernel::updateAttach(int mode, bool needFullRebuild) {
    // for each surface voxel pair, check distance < watchDistance, make a new
    // link between these two voxels, updateSurface().
    int blockSize;
    int minGridSize;
    tmpCollisionCount = 0;
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
    } else if (mode == 0) {

        // Pairwise detection O(n ^ 2)
        blockSize = 16;
        dim3 dimBlock(blockSize, blockSize);
        dim3 dimGrid((num_d_surface_voxels + dimBlock.x - 1) / dimBlock.x, (num_d_surface_voxels + dimBlock.y - 1) / dimBlock.y);
        gpu_update_attach<<<dimGrid, dimBlock>>>(d_surface_voxels, num_d_surface_voxels, watchDistance,
                this); // invoke two dimensional gpu threads 'CUDA C++ Programming
        // Guide', Nov 2019, P52.
        CUDA_CHECK_AFTER_CALL();
        VcudaDeviceSynchronize();
        
    } else if (mode == 1) {
        if (num_d_surface_voxels == 1) {
            return;
        }

        // sam:
        // bool needFullRebuild = false;  // is an input now

        if (d_collision_system->N != num_d_surface_voxels) {
            needFullRebuild = true;
            d_collision_system->set_num_objects_device(num_d_surface_voxels); // inform Collision System of new number of voxels. Allocate if needed.
        } else if (CurStepCount <= 1) {
            needFullRebuild = true;
        }
        // sam (end)

        // Tree based collision detection!
        int num_cols = 0;

        // copy position information for each voxel into the collision system.
        int num_SM, curDeviceId, gridSize, blockSize;
        cudaGetDevice(&curDeviceId);
        cudaDeviceGetAttribute(&num_SM, cudaDevAttrMultiProcessorCount, curDeviceId);
        blockSize = (num_d_surface_voxels-1)/num_SM;
        if (num_SM * blockSize < (num_d_surface_voxels)) {
            blockSize += 1;
        }
        if (blockSize > 256) {
            blockSize = 256;
            gridSize = (num_d_surface_voxels + 255)/256; 
        } else {
            gridSize = num_SM;
        }

        gpu_update_collision_system_pos_radius<<<gridSize, blockSize>>>(d_surface_voxels, num_d_surface_voxels, watchDistance, this);
        CUDA_CHECK_AFTER_CALL();
        VcudaDeviceSynchronize();

        // if number of surface voxels has changed, we need to re-init the collision detection tree.
        // Note that as voxels move, it makes sense to re-build the tree to improve the performance of tree traversal in the
        // find_collisions_device() method of the CollisionSystem, however rebuilding the tree is not necessary to have accurate collision detection.

        // if ( d_collision_system->N != num_d_surface_voxels || CurStepCount <= 1) {
        //     d_collision_system->N = num_d_surface_voxels;
        //     // d_collision_system->end = d_collision_system->start + num_d_surface_voxels; // sam
        //     //printf("Step %lu Updated Collision System surface voxel counts\n", CurStepCount);
        //     VcudaDeviceSynchronize();

        if ( needFullRebuild) { // sam
            d_collision_system->update_x_pos_ranks();
            d_collision_system->update_y_pos_ranks();
            d_collision_system->update_z_pos_ranks();
            d_collision_system->update_mortons();
            d_collision_system->build_tree();
        }

        if (CurStepCount%200 == 1 && CurStepCount > 100) {
            d_collision_system->update_x_pos_ranks();
        } else if (CurStepCount%200 == 2 && CurStepCount > 100) {
            d_collision_system->update_y_pos_ranks();
        } else if (CurStepCount%200 == 3 && CurStepCount > 100) {
            d_collision_system->update_z_pos_ranks();
        } else if (CurStepCount%200 == 4 && CurStepCount > 100) {
            d_collision_system->update_mortons();
        } else if (CurStepCount%200 == 5 && CurStepCount > 100) {
            d_collision_system->build_tree();
        }
        d_collision_system->update_bounding_boxes();
        num_cols = d_collision_system->find_collisions_device();

        if (num_cols == 0 ) { // no collisions were detected.
            return;
        } else if (num_cols < 0 ) {
            // ran out of pre-alloated memory for the collision system... :(
            printf(COLORCODE_BOLD_RED "CollisionSystem ran out of memory for collisions. Please set a higher value of MAX_COLLISIONS_PER_OBJECT\n" COLORCODE_RESET);
            assert(0);
        } else {
            int num_SM, curDeviceId, gridSize, blockSize;
            cudaGetDevice(&curDeviceId);
            cudaDeviceGetAttribute(&num_SM, cudaDevAttrMultiProcessorCount, curDeviceId);
            blockSize = (num_cols-1)/num_SM;
            if (num_SM * blockSize < (num_cols)) {
                blockSize += 1;
            }
            if (blockSize > 256) {
                blockSize = 256;
                gridSize = (num_cols + 255)/256; 
            } else {
                gridSize = num_SM;
            }
    
            gpu_update_sync_collisions<<<gridSize, blockSize>>>(d_surface_voxels, num_cols, watchDistance, this);
            CUDA_CHECK_AFTER_CALL();
            VcudaDeviceSynchronize();
            // printf("Step %lu with %d surface voxels and %d total voxels found %d real collisions and %d potential collisions\n", CurStepCount, num_d_surface_voxels, num_d_voxels, tmpCollisionCount,  num_cols);
            
            // for (int i = 0; i < 2* num_d_surface_voxels - 1; i++) {
            //     BoundingBox b = d_collision_system->bounding_boxes_d_ptr[i];
            //     printf("BoundingBox %d: [(%5.2f %5.2f), (%5.2f %5.2f), (%5.2f %5.2f)]\n", i, b.x_min, b.x_max, b.y_min, b.y_max, b.z_min, b.z_max);
            // }
            // printf("morton_numbers = [");
            // for (int i = 0; i < num_d_surface_voxels; i++) {
            //     printf("%lu, ", d_collision_system->mortons_d_ptr[i]);
            // }
            // printf("]\n\nleaf_parents = [");
            // for (int i = 0; i < num_d_surface_voxels; i++) {
            //     printf("%u, ", d_collision_system->leaf_parent_d_ptr[i]);
            // }
            // printf("]\n\ninternal_parent = [");
            // for (int i = 0; i < num_d_surface_voxels - 1; i++) {
            //     printf("%u, ", d_collision_system->internal_parent_d_ptr[i]);
            // }
            
            // printf("]\n\ninternal_childA = [");
            // for (int i = 0; i < num_d_surface_voxels - 1; i++) {
            //     printf("%u, ", d_collision_system->internal_childA_d_ptr[i]);
            // }
            
            // printf("]\n\ninternal_childB = [");
            // for (int i = 0; i < num_d_surface_voxels - 1; i++) {
            //     printf("%u, ", d_collision_system->internal_childB_d_ptr[i]);
            // }

            // printf("]\n\nBBoxFlags = [");
            // for (int i = 0; i < num_d_surface_voxels - 1; i++) {
            //     printf("%u, ", d_collision_system->internal_node_bbox_complete_flag_d_ptr[i]);
            // }
            
            // printf("]\n\n");

            // for (int obj_id = 0; obj_id < 10; obj_id++) {
            //     printf("Voxel %d: (%.4f %.4f %.4f %.4f)\n", obj_id, d_surface_voxels[obj_id]->pos.x, 
            //                                                         d_surface_voxels[obj_id]->pos.y,
            //                                                         d_surface_voxels[obj_id]->pos.z,
            //                                                         (float) d_surface_voxels[obj_id]->baseSizeAverage() * (float) COLLISION_ENVELOPE_RADIUS * 2.f);
            // }
            // printf("X   pos are sorted: %d\n", (int) thrust::is_sorted(thrust::device, d_collision_system->x_pos_d_ptr, d_collision_system->x_pos_d_ptr + d_collision_system->N));
            // printf("Y   pos are sorted: %d\n", (int) thrust::is_sorted(thrust::device, d_collision_system->y_pos_d_ptr, d_collision_system->y_pos_d_ptr + d_collision_system->N));
            // printf("Z   pos are sorted: %d\n", (int) thrust::is_sorted(thrust::device, d_collision_system->z_pos_d_ptr, d_collision_system->z_pos_d_ptr + d_collision_system->N));
            // printf("Mortons are sorted: %d\n", (int) thrust::is_sorted(thrust::device, d_collision_system->mortons_d_ptr, d_collision_system->mortons_d_ptr + d_collision_system->N));
            // assert(0);
            return;
        }
    }  else {
        printf("Please specify a mode in updateAttach(int mode)\n");
        assert(0); // Mode must be 0 or 1
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
    PRINT(this, "%f) regenerate surface voxels %d in %d. \n", currentTime, num_d_surface_voxels, num_d_voxels);
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
    d_targets.clear();
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

__device__ void VX3_VoxelyzeKernel::updateGroups() {
    // Sida: Dont update groups parallelly, instead, only update those pertain to collision (recorded in d_voxel_to_update_group and d_voxels_to_detach)
    if (d_voxels_to_detach.size()>0) {
        PRINT(this, "d_voxels_to_detach.size = %d.\n", d_voxels_to_detach.size());
        for (int i=0;i<d_voxels_to_detach.size();i++) {
            if (d_voxels_to_detach[i]->removed)
                continue;
            VX3_Voxel* voxel_to_detach = d_voxels_to_detach[i];
            // voxel_to_detach->removed = true;
            // Remove the original group, replacing by new groups.
            voxel_to_detach->d_group->removed = true;

            VX3_VoxelGroup *g = (VX3_VoxelGroup*)hamalloc(sizeof(VX3_VoxelGroup));
            g->deviceInit(this);
            // d_voxelgroups.push_back(g);
            g->d_voxels.push_back(voxel_to_detach);
            g->needUpdate = true;
            voxel_to_detach->d_group = g;
            d_voxel_to_update_group.push_back(voxel_to_detach);
            PRINT(this, "Detach voxel (%p) and remove group (%p (needUpdate %d)) with %d voxels in it.\n", voxel_to_detach, voxel_to_detach->d_group, voxel_to_detach->d_group->needUpdate, voxel_to_detach->d_group->d_voxels.size());
            // delete all the links as well
            for (int i = 0; i < 6; i++) {
                if (voxel_to_detach->links[i]) {
                    VX3_Voxel* neighbor = voxel_to_detach->adjacentVoxel((linkDirection)i);
                    VX3_VoxelGroup *g = (VX3_VoxelGroup*)hamalloc(sizeof(VX3_VoxelGroup));
                    g->deviceInit(this);
                    // d_voxelgroups.push_back(g);

                    g->d_voxels.push_back(neighbor);
                    g->needUpdate = true;
                    neighbor->d_group = g;
                    d_voxel_to_update_group.push_back(neighbor);
                    PRINT(this, "add to update waiting list: voxel (%p) group (%p) (from %p).\n", neighbor, g, this);

                    voxel_to_detach->links[i]->removed = true;
                    PRINT(this, "Remove link (%p) of detached voxel (%p).\n", voxel_to_detach->links[i], voxel_to_detach);
                    neighbor->links[oppositeDirection(i)] = NULL;
                    voxel_to_detach->links[i] = NULL;
                }
            }
        }
        d_voxels_to_detach.clear();
    }

    if (d_voxel_to_update_group.size()>0) {
        PRINT(this, "d_voxel_to_update_group.size = %d.\n", d_voxel_to_update_group.size());
        for (int i=0;i<d_voxel_to_update_group.size();i++) {
            if (d_voxel_to_update_group[i]->removed)
                continue; // removed voxel doesn't deserve a update.
            VX3_Voxel* v = d_voxel_to_update_group[i];

            PRINT(this, "i=%d, d_voxel_to_update_group.size()=%d.\n", i, d_voxel_to_update_group.size());
            PRINT(this, "(%p) to be update. %d voxels were in it's group (%p) (removed %d).\n", d_voxel_to_update_group[i], d_voxel_to_update_group[i]->d_group->d_voxels.size(), d_voxel_to_update_group[i]->d_group, d_voxel_to_update_group[i]->d_group->removed);
            
            d_voxel_to_update_group[i]->d_group->updateGroup(d_voxel_to_update_group[i]);
            PRINT(this, "after update, %d voxels are in the group (%p) (removed %d).\n", d_voxel_to_update_group[i]->d_group->d_voxels.size(), d_voxel_to_update_group[i]->d_group, d_voxel_to_update_group[i]->d_group->removed);
        }
        d_voxel_to_update_group.clear();
    }
}

__device__ void VX3_VoxelyzeKernel::surfaceGrow() {
    int minGridSize, blockSize;
    if (num_d_surface_voxels>0) {
        cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, gpu_surface_grow, 0,
                num_d_surface_voxels); // Dynamically calculate blockSize
        int gridSize_voxels = (num_d_surface_voxels + blockSize - 1) / blockSize;
        int blockSize_voxels = num_d_surface_voxels < blockSize ? num_d_surface_voxels : blockSize;
        gpu_surface_grow<<<gridSize_voxels, blockSize_voxels>>>(d_surface_voxels, num_d_surface_voxels);
        CUDA_CHECK_AFTER_CALL();
        VcudaDeviceSynchronize();
    }
}

__device__ bool VX3_VoxelyzeKernel::ThoroughValidationCheck() {
    VX3_OnlineTest tester;
    return tester.ThoroughTest(this);
}

/* Sub GPU Threads */
__global__ void gpu_surface_grow(VX3_Voxel** surface_voxels, int num) {
    int gindex = threadIdx.x + blockIdx.x * blockDim.x;
    if (gindex < num) {
        VX3_Voxel* v = surface_voxels[gindex];
        v->grow();
    }
}

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
__device__ bool is_neighbor(VX3_Voxel *voxel1, VX3_Voxel *voxel2, int depth) {
    if (voxel1 == voxel2) {
        return true;
    }
    if (voxel1==NULL || depth <= 0) { // cannot find in depth
        return false;
    }
    for (int i = 0; i < 6; i++) {
        if ( is_neighbor(voxel1->adjacentVoxel((linkDirection)i), voxel2, depth-1)) {
            return true;
        }
    }
    return false;
}

__device__ void handle_collision_attachment(VX3_Voxel *voxel1, VX3_Voxel *voxel2, double watchDistance, VX3_VoxelyzeKernel *k) {
    // if ((voxel1->ix==47 && voxel1->iy==101) || (voxel2->ix==47 && voxel2->iy==101)) {
    //     printf("break.\n");
    // }
    // if both of the voxels are fixed, no need to compute.
    if (voxel1->mat->fixed && voxel2->mat->fixed)
        return;

    // to exclude voxels already have link between them. check in depth of
    // 1, direct neighbor ignore the collision
    if (is_neighbor(voxel1, voxel2, 1)) {
        return;
    }

    // to check for distance: collision system provides a list of potentionly collided pairs, we need to check the distance here.
    VX3_Vec3D<double> diff = voxel1->pos - voxel2->pos;

    // Disable dynamical watch distance to save time
    // watchDistance = (voxel1->baseSizeAverage() + voxel2->baseSizeAverage()) * COLLISION_ENVELOPE_RADIUS * watchDistance;

    if (diff.x > k->staticWatchDistance || diff.x < -k->staticWatchDistance)
        return;
    if (diff.y > k->staticWatchDistance || diff.y < -k->staticWatchDistance)
        return;
    if (diff.z > k->staticWatchDistance || diff.z < -k->staticWatchDistance)
        return;

    if (diff.Length2() > k->staticWatchDistance_square)
        return;


    // calculate and store contact force, apply and clean in VX3_Voxel::force()
    VX3_Vec3D<> cache_contactForce1, cache_contactForce2;
    if (k->EnableCollision) {
        VX3_Collision collision(voxel1, voxel2);
        collision.updateContactForce();
        cache_contactForce1 = collision.contactForce(voxel1);
        cache_contactForce2 = collision.contactForce(voxel2);
        // voxel1->contactForce += cache_contactForce1;
        // voxel2->contactForce += cache_contactForce2;
        
        atomicAdd(&(voxel1->contactForce.x), cache_contactForce1.x);
        atomicAdd(&(voxel1->contactForce.y), cache_contactForce1.y);
        atomicAdd(&(voxel1->contactForce.z), cache_contactForce1.z);
        
        atomicAdd(&(voxel2->contactForce.x), cache_contactForce2.x);
        atomicAdd(&(voxel2->contactForce.y), cache_contactForce2.y);
        atomicAdd(&(voxel2->contactForce.z), cache_contactForce2.z);

        atomicAdd(&(k->tmpCollisionCount), 1);
        auto v1pos = voxel1->pos;
        auto v2pos = voxel2->pos;

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
    if (k->enableAttach) {
        if (k->d_attach_manager->attachWhileCollide(voxel1, voxel2)) {
            voxel1->contactForce -= cache_contactForce1;
            voxel2->contactForce -= cache_contactForce2;
        }
    }
    return;

}

/**
 *
 * Updates the collision system information about the surface voxels position and size.
 *
 * Casts position and radius information to floats for faster processing of collisions.
 *
 * @param surface_voxels list of pointers to surface voxels
 * @param num number of surface voxels
 * @param watchDistance  How close do two voxels need to be to track connections
 * @param k A device pointer to a VX3_VoxelyzeKernel
 */
__global__ void gpu_update_collision_system_pos_radius(VX3_Voxel **surface_voxels, int num, double watchDistance, VX3_VoxelyzeKernel *k) {
    unsigned int voxelId = threadIdx.x + blockIdx.x * blockDim.x;
    if (voxelId < num) {
        // printf("Updating Voxel %u (%u %u %u)\n",  voxelId, threadIdx.x, blockIdx.x, blockDim.x);
        float x, y, z, r;
        auto currentVoxel = surface_voxels[voxelId];
        auto pos = currentVoxel->pos;
        x = (float) pos.x;
        y = (float) pos.y;
        z = (float) pos.z;
        r = (float) currentVoxel->baseSizeAverage() * (float) COLLISION_ENVELOPE_RADIUS * 1.001f;
        assert (r > 0);
        k->d_collision_system->x_pos_d_ptr[voxelId] = x;
        k->d_collision_system->y_pos_d_ptr[voxelId] = y;
        k->d_collision_system->z_pos_d_ptr[voxelId] = z;
        k->d_collision_system->radius_d_ptr[voxelId] = r;
    }
}

__global__ void gpu_update_sync_collisions(VX3_Voxel **surface_voxels, int num, double watchDistance, VX3_VoxelyzeKernel *k) {
    unsigned int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid < num) {
        Collision col = k->d_collision_system->collisions_d_ptr[tid];
        VX3_Voxel *voxel1 = surface_voxels[col.a];
        VX3_Voxel *voxel2 = surface_voxels[col.b];
        if (voxel1->removed || voxel2->removed) {
            return;
        }
        handle_collision_attachment(voxel2, voxel1, watchDistance, k);
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
        // if (k->RandomizeCiliaEvery>0) { // sam
        //     surface_voxels[index]->CiliaForce = surface_voxels[index]->orient.RotateVec3D(
        //         surface_voxels[index]->baseCiliaForce * surface_voxels[index]->randCiliaCoef * -1);
        // }
        // else {
        surface_voxels[index]->CiliaForce = surface_voxels[index]->orient.RotateVec3D(
            surface_voxels[index]->baseCiliaForce + surface_voxels[index]->localSignal * surface_voxels[index]->shiftCiliaForce);
        // }

    }
}

// sam:
__global__ void gpu_update_brownian_motion(VX3_Voxel **surface_voxels, int num, int WorldSize, double seed, double currentTime, VX3_VoxelyzeKernel *k) {
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    if (index < num) {
        if (surface_voxels[index]->removed)
            return;
        if (surface_voxels[index]->mat->Cilia == 0)
            return;
        if (surface_voxels[index]->mat->TurnOnCiliaAfterThisManySeconds > k->currentTime)
            return;
        
        // randomize according to seed, timestep and original position in the grid (all of this is repeatable)
        int ix = surface_voxels[index]->indexX();
        int iy = surface_voxels[index]->indexY();
        int iz = surface_voxels[index]->indexZ();
        int randIndex = ix + WorldSize*iy + WorldSize*WorldSize*iz;

        curandState state;
        // curand_init(seed + currentTime, index, 0, &state);
        curand_init(seed + currentTime, randIndex, 0, &state);
        // surface_voxels[index]->randCiliaCoef = curand_uniform(&state);
        surface_voxels[index]->baseCiliaForce.x = 2*curand_uniform(&state)-1;
        surface_voxels[index]->baseCiliaForce.y = 2*curand_uniform(&state)-1;
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

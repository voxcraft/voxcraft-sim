//
// Created by Sida Liu
//  This class is the backbone of one simulation. At each time step, it starts multiple threads to handle calculation of all aspects.
//
#if !defined(VX3_VOXELYZE_KERNEL_H)
#define VX3_VOXELYZE_KERNEL_H
#include "VX3.cuh"
#include "VX_Sim.h" //readVXA
#include <map>

#include "VX3_vector.cuh"

#include "VX3_ForceField.h"

#include "VX3_Collision.h"
#include "VX3_Link.h"
#include "VX3_MaterialLink.h"
#include "VX3_Voxel.h"
#include "VX_Enums.h"
#include "VX3_AttachManager.h"
#include "VX3_GrowthManager.h"
#include "VX3_OnlineTest.h"

#include "../Cu-Collision-Detection/include/CollisionSystem.cuh"

/*
 * VX3_VoxelyzeKernel is a GPU mock class of CVoxelyze
 * Usage: setup a CVoxelyze, and use the constructor function to initialize a VX3_VoxelyzeKernel.
 */
class VX3_VoxelyzeKernel {
  public:
    /* Host methods */
    VX3_VoxelyzeKernel(CVX_Sim *In);
    VX3_VoxelyzeKernel() = default; // start from scratch, read VXA ourselves.

    void cleanup();

    __device__ void deviceInit();

    /* Cuda methods */
    __device__ bool doTimeStep(float dt = -1.0f);
    __device__ double recommendedTimeStep();
    __device__ void updateCurrentCenterOfMass();
    __device__ bool StopConditionMet();
    __device__ void updateTemperature();
    __device__ void updateAttach(int mode, bool needFullRebuild = false);  // sam
    __device__ void updateDetach();
    __device__ void regenerateSurfaceVoxels();
    __device__ VX3_MaterialLink *combinedMaterial(VX3_MaterialVoxel *mat1, VX3_MaterialVoxel *mat2);
    __device__ void computeFitness();
    __device__ void registerTargets();
    __device__ void computeTargetCloseness();
    __device__ void saveInitialPosition();

    // for Secondary Experiment
    __device__ void removeVoxels();
    __device__ void InitializeCenterOfMass();
    __device__ bool EarlyStopIfNoBotsRemain(); // sam
    __device__ void replenishMaterial(int start, int end, int step, int mat); // sam
    __device__ void convertMatIfSmallBody(int mat1, int mat2, bool convertSingletons); // sam
    __device__ void convertMatIfLargeBody(int mat1, int mat2); // sam
    
    __device__ void reInitAllGroups(); // sam

    __device__ bool addVoxel(int x, int y, int z, int mat); // sam

    // for Testing
    // check all the voxels, links and voxelgroups for validation.
    __device__ bool ThoroughValidationCheck();
    
    /* data */
    bool forceExit = false;
    char vxa_filename[256];
    double voxSize;            // lattice size
    double currentTime = 0.0f; // current time of the simulation in seconds
    double OptimalDt = 0.0f;
    double DtFrac;
    StopCondition StopConditionType;
    double StopConditionValue;
    unsigned long CurStepCount = 0;

    bool enableFloor = true;

    // Temperature:
    bool TempEnabled;                           // overall flag for temperature calculations
    bool VaryTempEnabled;                       // is periodic variation of temperature on?
    double TempBase, TempAmplitude, TempPeriod; // degress celcius
    // double currentTemperature; //updated based on time... (for phase 0... individual materials now have their own current temp

    VX3_Vec3D<double> currentCenterOfMass;
    VX3_Vec3D<double> initialCenterOfMass;

    std::vector<CVX_Link *> h_links;
    VX3_Link *d_links;
    int num_d_links;
    VX3_dVector<VX3_Link *> d_v_links;
    VX3_hdVector<VX3_Link *> hd_v_links;
    std::map<CVX_Link *, VX3_Link *> h_lookup_links;

    std::vector<CVX_Voxel *> h_voxels;
    VX3_Voxel *d_voxels;
    int num_d_voxels;
    int num_d_init_voxels;
    int MaxNewVoxelsAddedMidSim = 10000; // sam: pre-allocate memory for this many new voxels added mid sim
    VX3_Voxel **d_surface_voxels; // an array of pointer d_surface_voxels[i] -> d_voxels[j]
    int num_d_surface_voxels;
    std::map<CVX_Voxel *, VX3_Voxel *> h_lookup_voxels;

    std::vector<CVX_MaterialVoxel *> h_voxelMats;
    VX3_MaterialVoxel *d_voxelMats;
    int num_d_voxelMats;

    std::vector<CVX_MaterialLink *> h_linkMats;
    VX3_MaterialLink *d_linkMats;
    int num_d_linkMats;
    VX3_dVector<VX3_MaterialLink *>
        d_v_linkMats; // d_v_linkMats initially stores pointer -> d_linkMats[j], then during attachment, new linkMats will be push_back()ed.
    VX3_hdVector<VX3_MaterialLink *> hd_v_linkMats; // used to pass vector to kernel and initiate d_v_linkMats.

    // Collision Constants
    float boundingRadius; //(in voxel units) radius to collide a voxel at
    float watchDistance;  //(in voxel units) Distance between voxels (not including 2*boundingRadius for each voxel) to watch for collisions from.

    // bool* d_collisionsStale;
    CollisionSystem *d_collision_system;
    CollisionSystem *h_collision_system;
    VX3_dVector<VX3_Collision *> d_v_collisions;

    bool enableAttach;
    bool enableDetach;
    bool keepJustOneIfManyHaveSameGroupPosition = false; // sam
    bool EnableCollision = true;
    int CollisionMode = 1;
    int RecordStepSize = 0;
    int RecordLink = 0;
    int RecordVoxel = 0;
    int RecordFixedVoxels = 1; // sam
    int SurfaceVoxelsOnly = 1;

    // Safety Guard during the creation of new link
    int SafetyGuard = 500;

    // Force Field
    VX3_ForceField force_field;

    VX3_MathTreeToken fitness_function[1024];
    double fitness_score = 0;

    VX3_MathTreeToken AttachCondition[5][1024];
    VX3_MathTreeToken StopConditionFormula[1024];

    int collisionCount = 0;
    int tmpCollisionCount = 0;

    //Calculate Angle
    //A---B----C
    //A: currentCenterOfMass_history[0]
    //B: currentCenterOfMass_history[1]
    //C: currentCenterOfMass
    VX3_Vec3D<double> currentCenterOfMass_history[2];
    int angleSampleTimes = 0;
    double recentAngle = 0;

    int EnableTargetCloseness = 0;
    int SavePositionOfAllVoxels = 0;
    int EnableCilia = 0;
    double RandomizeCiliaEvery = 0;  // sam
    int EnableSignals = 0;
    double  targetCloseness = 0;
    VX3_dVector<VX3_Voxel*> d_targets;
    int numClosePairs = 0;
    bool isSurfaceChanged=false;
    double MaxDistInVoxelLengthsToCountAsPair=0.0;

    //Spatial Hash
    //index all surface voxels into grid, so we only need to compare neighbor grids for detecting collision
    //dx,dy,dz width height of one single grid. must larger than collision watch distance.
    VX3_Vec3D<> gridLowerBound;
    VX3_Vec3D<> gridUpperBound;
    VX3_Vec3D<> gridDelta;
    //number of grid on each side
    int lookupGrid_n = 10;
    //total number of grid
    int num_lookupGrids;
    VX3_dVector<VX3_Voxel*>* d_collisionLookupGrid;

    VX3_Vec3D<>* d_initialPosition = NULL;

    //for Secondary Experiment
    int SecondaryExperiment = 0;
    int SelfReplication = 0;
    double ReinitializeInitialPositionAfterThisManySeconds = 0.0;
    double SettleTimeBeforeNextRoundOfReplication = 0.0;  // sam
    bool InitialPositionReinitialized = true;  // sam
    int MinimumBotSize = 2; // sam

    double lastReplicationTime = 0.0; // sam

    int EnableExpansion=0;

    VX3_AttachManager* d_attach_manager;

    int mutexRotateSingleton=0;

    double lastRegenerateSurfaceTime = 0;
    
    // Using static watch distance and caching it improves performance.
    double staticWatchDistance = 0;
    double staticWatchDistance_square = 0;

    VX3_GrowthManager* d_growth_manager;

    int EnableSurfaceGrowth = 0;
    double SurfaceGrowth_Interval = 1;
    double SurfaceGrowth_activeTime = 0;
    double SurfaceGrowth_Rate = 1;
    int SurfaceGrowth_Growed = 0;
    RandomGenerator* randomGenerator;

    VX3_dVector<VX3_VoxelGroup *> d_voxelgroups;
    __device__ void updateGroups();

    __device__ void surfaceGrow();

    // To remember what voxels their groups should be update. and at the end of the timestep, update them one-by-one, sequentially, not in parallel.
    VX3_dVector<VX3_Voxel *> d_voxel_to_update_group;
    VX3_dVector<VX3_Voxel *> d_voxel_to_absorb;

    bool VerboseMode;
    bool SkipThoroughTest;
    unsigned int ThoroughTestStepSize;
    unsigned int ThoroughTestStartAt;

    int GPU_id;
};

#endif // VX3_VOXELYZE_KERNEL_H

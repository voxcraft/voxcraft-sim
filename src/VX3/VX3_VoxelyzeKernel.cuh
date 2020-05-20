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

    /* Cuda methods */
    __device__ bool doTimeStep(float dt = -1.0f);
    __device__ double recommendedTimeStep();
    __device__ void updateCurrentCenterOfMass();
    __device__ bool StopConditionMet();
    __device__ void updateTemperature();
    __device__ void syncVectors();
    __device__ void updateAttach();
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

    /* data */
    bool forceExit = false;
    char vxa_filename[256];
    double voxSize;            // lattice size
    double currentTime = 0.0f; // current time of the simulation in seconds
    double OptimalDt = 0.0f;
    double DtFrac;
    StopCondition StopConditionType;
    double StopConditionValue;
    unsigned long CurStepCount = 0.0f;

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
    VX3_dVector<VX3_Collision *> d_v_collisions;

    bool enableAttach;
    bool enableDetach;
    bool EnableCollision = true;
    int RecordStepSize = 0;
    int RecordLink = 0;
    int RecordVoxel = 0;

    // Safety Guard during the creation of new link
    int SafetyGuard = 500;

    // Force Field
    VX3_ForceField force_field;

    VX3_MathTreeToken fitness_function[1024];
    double fitness_score = 0;

    VX3_MathTreeToken AttachCondition[5][1024];
    VX3_MathTreeToken StopConditionFormula[1024];

    int collisionCount = 0;

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
    double ReinitializeInitialPositionAfterThisManySeconds = 0.0;
    bool InitialPositionReinitialized = false;

    int EnableExpansion=0;
};

#endif // VX3_VOXELYZE_KERNEL_H

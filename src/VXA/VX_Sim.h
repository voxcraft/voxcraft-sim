/*******************************************************************************
Copyright (c) 2010, Jonathan Hiller (Cornell University)
If used in publication cite "J. Hiller and H. Lipson "Dynamic Simulation of Soft Heterogeneous Objects" In press. (2011)"

This file is part of Voxelyze.
Voxelyze is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Voxelyze is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
See <http://www.opensource.org/licenses/lgpl-3.0.html> for license details.
*******************************************************************************/

//Branched copy placeholder

#ifndef VX_SIM_H
#define VX_SIM_H

#include "VX_Environment.h"
#include "VX_MeshUtil.h"
#include "VX_Enums.h"
#include <deque>

#define VX2 //use VX2 instead of VX1

#include "Voxelyze.h"
#include "VX_Material.h"

#define MIN_TEMP_FACTOR 0.1 /*The smallest a voxel can go. (for stability reasons)*/
#define HISTORY_SIZE 500

struct SimState { //Information about current simulation state:
	void Clear() {CurCM = TotalObjDisp = Vec3D<>(0,0,0); NormObjDisp = MaxVoxDisp = MaxVoxVel = MaxVoxKinE = MaxBondStrain = MaxBondStress = MaxBondStrainE = TotalObjKineticE = TotalObjStrainE = MaxPressure = MinPressure = 0.0;}
	Vec3D<> CurCM;
	Vec3D<> TotalObjDisp; //a vector total of the magnitude of displacements
	vfloat NormObjDisp; //reduced to a scalar (magnitude) 
	vfloat MaxVoxDisp, MaxVoxVel, MaxVoxKinE, MaxBondStrain, MaxBondStress, MaxBondStrainE, MaxPressure, MinPressure;
	vfloat TotalObjKineticE, TotalObjStrainE;
};

//!Dynamic simulation class for time simulation of voxel objects.
class CVX_Sim
{
public:
	CVX_Sim(void); //!< Constructor
	~CVX_Sim(void); //!< Destructor

	CVoxelyze Vx;
	void CopyMat(CVXC_Material* pOld, CVX_Material* pNew); //copies parameters from pOld to pNew

	bool LoadVXAFile(std::string filename, std::string* pRetMsg = NULL);
	bool ReadVXA(CXML_Rip* pXML, std::string* RetMessage = NULL);
	bool ReadXML(CXML_Rip* pXML, std::string* RetMessage = NULL);

	//input information
	CVX_Environment* pEnv; //!< Pointer to the physical environment information. This variable is set on import() and should not be manually changed.
	CVX_Object LocalVXC; //!< Local copy of the voxel object. This copy is stored to ensure it never changes throughout the simulation. This variable is set on import() and should not be manually changed.

	//Simulation Management
	bool Import(CVX_Environment* pEnvIn = NULL, CMesh* pSurfMeshIn = NULL, std::string* RetMessage = NULL); //!< Imports a physical environment into the simulator.

	bool IsInitalized(void) const {return Initalized;} //!< Returns true if a valid environment has been imported.

	void ClearAll(void); //!< Clears all environment-specific information form the simulation.

	//Integration/simulation running
	vfloat DtFrac = 0.9; //percent of maximum dt to use
	vfloat OptimalDt; //calculated optimal dt
	vfloat dt; //actual seconds per timestep
	vfloat CurTime; //The current absolute time of the simulation in seconds.
	int CurStepCount;

	//Simulator features:

	void EnableFeature(const int VXSFEAT, bool Enabled=true);
	bool IsFeatureEnabled(const int VXSFEAT);

	//Damping:
	void SetBondDampZ(vfloat BondDampZIn) {BondDampingZ = BondDampZIn; for (int i=0; i<Vx.materialCount(); i++) Vx.material(i)->setInternalDamping(BondDampZIn);} //!< Sets the damping ratio for connected voxels. When this is non-zero, each voxel is damped (based on its mass and stiffness) according to its relative velocity to the other voxel in each bond. Range is [0.0 to 1.0]. Values greater than 1.0 may cause numerical instability.
	void SetCollisionDampZ(vfloat ColDampZIn) {ColDampingZ = ColDampZIn; 	for (int i=0; i<Vx.materialCount(); i++) Vx.material(i)->setCollisionDamping(ColDampZIn);} //!< Sets the damping ratio for voxels in colliding state. When this is non-zero, each voxel is damped (based on its mass and stiffness) according to the penetration velocity. Range is [0.0 to 1.0]. Values greater than 1.0 may cause numerical instability.
	void SetSlowDampZ(vfloat SlowDampIn) {SlowDampingZ = SlowDampIn;  for (int i=0; i<Vx.materialCount(); i++) Vx.material(i)->setGlobalDamping(SlowDampIn);} //!< Sets the damping ratio that slows downs voxels. When this is non-zero, each voxel is damped (based on its mass and stiffness) to ground. Range is [0.0 to 1.0]. Values greater than 1.0 may cause numerical instability.
	
	vfloat GetBondDampZ(void){if (Vx.materialCount()>0) return Vx.material(0)->internalDamping(); else return BondDampingZ;} //!< Returns the current bond damping.
	vfloat GetCollisionDampZ(void) {if (Vx.materialCount()>0) return Vx.material(0)->collisionDamping(); else return ColDampingZ;} //!< Returns the current collision damping.
	vfloat GetSlowDampZ(void) {if (Vx.materialCount()>0) return Vx.material(0)->globalDamping(); else return SlowDampingZ;} //!< Returns the current voxel slowing damping.

	void SetGravityAccel(float grav);
	float GetGravityAccel(void);

	//Material Blending
	Vec3D<> MixRadius;
	MatBlendModel BlendModel;
	vfloat PolyExp; //polynomial exponent if using polynomial model

	//Equlibrium mode
	void ZeroAllMotion(void);
	bool MotionZeroed;

	//Stop conditions
	void SetStopConditionType(StopCondition StopConditionTypeIn = SC_NONE) {StopConditionType = StopConditionTypeIn;}
	void SetStopConditionValue(vfloat StopConditionValueIn = 0.0) {StopConditionValue = StopConditionValueIn;} 
	StopCondition GetStopConditionType(void){return StopConditionType;}
	vfloat GetStopConditionValue(void){return StopConditionValue;}
	bool StopConditionMet(void); //have we met the stop condition yet?

	//Information about current state:
	SimState SS;
	int StatToCalc;

	Vec3D<> GetCM(void);
	Vec3D<> IniCM; //initial center of mass

	int GetNumTouchingFloor();

	Vec3D<> GetSumForce(CVX_FRegion* pRegion); //returns total force on a region
	vfloat GetSumForceDir(CVX_FRegion* pRegion); //returns total force on a region in the direction of its displacement
	Vec3D<> GetAvgDisplace(CVX_FRegion* pRegion); //returns the average displacement in x,y,z of a region

		//temperature
	void UpdateMatTemps(void); //updates expansions for each material
	void UpdateMuMemory(void); //updates array for turning poissons ratio on and off.

protected:
	bool Initalized; //!< Flag to denote if simulation is runnable. True if there is an environement successfully loaded, false otherwise.

	//Integration
	//bool Integrate();
	bool UpdateStats(std::string* pRetMessage = NULL); //returns false if simulation diverged...

	StopCondition StopConditionType;
	vfloat StopConditionValue;

	void EnableEquilibriumMode(bool Enabled);
	void EnableVolumeEffects(bool Enabled);
	bool IsVolumeEffectsEnabled();

	int CurSimFeatures;

	//Damping:
	vfloat BondDampingZ; //Damping factor zeta (1 = critical damping, 0 = no damping...
	vfloat ColDampingZ;	//Damping factor zeta (1 = critical damping, 0 = no damping...
	vfloat SlowDampingZ; //Damping factor zeta (1 = critical damping, 0 = no damping...

	//Equilibrium mode
	bool KineticEDecreasing(void); //returns true if kinetic energy of the structure is decreasing
	vfloat MemBondDampZ, MemSlowDampingZ;
	//bool MemMaxVelEnabled;

	void ClearHistories() {std::fill(KinEHistory.begin(), KinEHistory.end(), -1); std::fill(TotEHistory.begin(), TotEHistory.end(), -1); std::fill(MaxMoveHistory.begin(), MaxMoveHistory.end(), -1);}
	std::deque<vfloat> KinEHistory; //value [0] is the newest...
	std::deque<vfloat> TotEHistory;
	std::deque<vfloat> MaxMoveHistory;

	std::vector<float> muMemory;

	public:
	CMesh* ImportSurfMesh; //local copy of any imported geometry surface mesh (optional)

	friend class VX3_SimulationManager;
	friend class VX3_VoxelyzeKernel;

};

#endif //VX_SIM_H

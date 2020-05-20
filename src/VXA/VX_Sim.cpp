/*******************************************************************************
Copyright (c) 2010, Jonathan Hiller (Cornell University)
If used in publication cite "J. Hiller and H. Lipson "Dynamic Simulation of Soft Heterogeneous Objects" In press. (2011)"

This file is part of Voxelyze.
Voxelyze is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your option) any later version. Voxelyze is distributed in the hope
that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU Lesser General Public License for more details. See <http://www.opensource.org/licenses/lgpl-3.0.html> for license details.
*******************************************************************************/

#include "VX_Sim.h"
#include "VXS_Bond.h"
#include "VXS_Voxel.h"

#include <sstream>

#include "VX_Material.h"
#include "VX_Voxel.h"

CVX_Sim::CVX_Sim(void) // : VoxelInput(this), BondInput(this) // : out("Logfile.txt", std::ios::ate)
{
    ImportSurfMesh = NULL;
    MotionZeroed = false;

    CurSimFeatures = VXSFEAT_PLASTICITY | VXSFEAT_FAILURE;

    MixRadius = Vec3D<>(0.0, 0.0, 0.0);
    BlendModel = MB_LINEAR;
    PolyExp = 1.0;

    KinEHistory.resize(HISTORY_SIZE, -1.0);
    TotEHistory.resize(HISTORY_SIZE, -1.0);
    MaxMoveHistory.resize(HISTORY_SIZE, -1.0);

    SetStopConditionType();
    SetStopConditionValue();

    StatToCalc = CALCSTAT_ALL;

    ClearAll();
}

CVX_Sim::~CVX_Sim(void) { ClearAll(); }

bool CVX_Sim::LoadVXAFile(std::string filename, std::string *pRetMsg) {
    CXML_Rip XML;
    if (!XML.LoadFile(filename, pRetMsg))
        return false;
    ReadVXA(&XML, pRetMsg);
    return true;
}

/*! The environment should have been previously initialized and linked with a single voxel object.
This function sets or resets the entire simulation with the new environment.
@param[in] pEnvIn A pointer to initialized CVX_Environment to import into the simulator.
@param[out] RetMessage A pointer to initialized string. Output information from the Import function is appended to this string.
*/
bool CVX_Sim::Import(CVX_Environment *pEnvIn, CMesh *pSurfMeshIn, std::string *RetMessage) {
    if (pEnvIn != NULL)
        pEnv = pEnvIn;

    Initalized = false;
    Vx.clear();
    LocalVXC = *pEnv->pObj; // make a copy of the reference digital object!

    EnableFeature(VXSFEAT_GRAVITY, pEnv->IsGravityEnabled());
    SetGravityAccel(pEnv->GetGravityAccel());
    EnableFeature(VXSFEAT_FLOOR, pEnv->IsFloorEnabled());
    EnableFeature(VXSFEAT_TEMPERATURE, pEnv->IsTempEnabled());
    EnableFeature(VXSFEAT_TEMPERATURE_VARY, pEnv->IsTempVaryEnabled());

    // transfer the materials in (and map indices)...
    std::vector<CVX_Material *> VxcToVx2MatIndexMap;                  // size of VXC materials, index is temporary to VX2 mat
    VxcToVx2MatIndexMap.resize(LocalVXC.GetNumMaterials() - 1, NULL); // skip erase
    muMemory.clear();
    for (int i = 1; i < LocalVXC.GetNumMaterials(); i++) { // for each material
        if (LocalVXC.GetBaseMat(i)->GetMatType() == SINGLE)
            VxcToVx2MatIndexMap[i - 1] = Vx.addMaterial();
        CopyMat(LocalVXC.GetBaseMat(i), VxcToVx2MatIndexMap[i - 1]);
        VxcToVx2MatIndexMap[i - 1]->setInternalDamping(BondDampingZ);
        VxcToVx2MatIndexMap[i - 1]->setGlobalDamping(SlowDampingZ);
        VxcToVx2MatIndexMap[i - 1]->setCollisionDamping(ColDampingZ);

        muMemory.push_back(LocalVXC.GetBaseMat(i)->GetPoissonsRatio()); // remember for toggleing volume effects on and off.
    }

    // add the voxels
    Vx.setVoxelSize(LocalVXC.GetLatDimEnv().x);
    int x, y, z;
    std::vector<CVX_Voxel *> VoxList;                     // local list of all voxel pointers
    for (int i = 0; i < LocalVXC.GetStArraySize(); i++) { // for each voxel in the array
        int VxcMatIndex = LocalVXC.GetMat(i) - 1;
        if (VxcMatIndex >= 0) {
            LocalVXC.GetXYZNom(&x, &y, &z, i);
            if (VxcMatIndex >= VxcToVx2MatIndexMap.size()) {
                printf("ERROR: No such material.\n");
            }
            VoxList.push_back(Vx.setVoxel(VxcToVx2MatIndexMap[VxcMatIndex], x, y, z));
            // PhaseOffset
            VoxList.back()->phaseOffset = pEnv->pObj->Structure.GetPhaseOffset(VoxList.size() - 1);
            // BaseCiliaForce
            VoxList.back()->baseCiliaForce = pEnv->pObj->Structure.GetBaseCiliaForce(VoxList.size() - 1);
            // ShiftCiliaForce
            VoxList.back()->shiftCiliaForce = pEnv->pObj->Structure.GetShiftCiliaForce(VoxList.size() - 1);
        }
    }

    // set any boundary conditions (this can be optimized much better
    int NumBCs = pEnv->GetNumBCs();
    std::vector<int> Sizes(NumBCs, 0);
    for (int i = 0; i < NumBCs; i++)
        Sizes[i] = pEnv->GetNumTouching(i); // count the number of voxels touching each bc
    Vec3D<> BCsize = pEnv->pObj->GetLatDimEnv() / 2.0;
    Vec3D<> WSSize = pEnv->pObj->GetWorkSpace();
    int numVox = VoxList.size();
    for (int i = 0; i < numVox; i++) {
        CVX_Voxel *pThisVox = VoxList[i];
        Vec3D<> ThisPos = pThisVox->position() + LocalVXC.GetLatDimEnv() / 2;
        for (int j = 0; j < NumBCs; j++) { // go through each primitive defined as a constraint!
            CVX_FRegion *pCurBc = pEnv->GetBC(j);
            char ThisDofFixed = pCurBc->DofFixed;
            if (pCurBc->GetRegion()->IsTouching(&ThisPos, &BCsize, &WSSize)) { // if this point is within
                if (IS_FIXED(DOF_X, ThisDofFixed))
                    pThisVox->external()->setDisplacement(X_TRANSLATE, pCurBc->Displace.x);
                if (IS_FIXED(DOF_Y, ThisDofFixed))
                    pThisVox->external()->setDisplacement(Y_TRANSLATE, pCurBc->Displace.y);
                if (IS_FIXED(DOF_Z, ThisDofFixed))
                    pThisVox->external()->setDisplacement(Z_TRANSLATE, pCurBc->Displace.z);
                if (IS_FIXED(DOF_TX, ThisDofFixed))
                    pThisVox->external()->setDisplacement(X_ROTATE, pCurBc->AngDisplace.x);
                if (IS_FIXED(DOF_TY, ThisDofFixed))
                    pThisVox->external()->setDisplacement(Y_ROTATE, pCurBc->AngDisplace.y);
                if (IS_FIXED(DOF_TZ, ThisDofFixed))
                    pThisVox->external()->setDisplacement(Z_ROTATE, pCurBc->AngDisplace.z);

                pThisVox->external()->addForce((Vec3D<float>)(pThisVox->external()->force() + pCurBc->Force / Sizes[j]));
                pThisVox->external()->addMoment((Vec3D<float>)(pThisVox->external()->moment() + pCurBc->Torque / Sizes[j]));
                // pThisVox->setExternalForce((Vec3D<float>)(pThisVox->externalForce() + pCurBc->Force/Sizes[j]));
                // pThisVox->setExternalMoment((Vec3D<float>)(pThisVox->externalMoment() + pCurBc->Torque/Sizes[j]));
            }
        }
    }

    OptimalDt = Vx.recommendedTimeStep(); // to set up dialogs parameter ranges, we need this before the first iteration.
    if (IsFeatureEnabled(VXSFEAT_TEMPERATURE))
        UpdateMatTemps();
    EnableVolumeEffects(IsFeatureEnabled(VXSFEAT_VOLUME_EFFECTS));

    Initalized = true;

    return true;
}

bool CVX_Sim::ReadVXA(CXML_Rip *pXML, std::string *RetMessage) // pointer to VXA element
{
    // pObj->ClearMatter();
    std::string ThisVersion = "1.1";
    std::string Version;
    pXML->GetElAttribute("Version", &Version);
    if (atof(Version.c_str()) > atof(ThisVersion.c_str()))
        if (RetMessage)
            *RetMessage +=
                "Attempting to open newer version of VXA file. Results may be unpredictable.\nUpgrade to newest version of VoxCAD.\n";

    if (pXML->FindElement("Simulator")) {
        ReadXML(pXML);
        pXML->UpLevel();
    }

    // load environment
    if (pEnv && pXML->FindElement("Environment")) {
        pEnv->ReadXML(pXML);
        pXML->UpLevel();
    }

    // Load VXC if pObj is valid...
    if (pEnv->pObj && (pXML->FindElement("VXC") || pXML->FindElement("DMF"))) {
        pEnv->pObj->ReadXML(pXML, false, RetMessage);
        pXML->UpLevel();
    }
    return true;
}

bool CVX_Sim::ReadXML(CXML_Rip *pXML, std::string *RetMessage) {
    int tmpInt;
    vfloat tmpVFloat;
    bool tmpBool;

    if (pXML->FindElement("Integration")) {
        // if (pXML->FindLoadElement("Integrator", &tmpInt)) CurIntegrator = (IntegrationType)tmpInt; else CurIntegrator = I_EULER;
        if (!pXML->FindLoadElement("DtFrac", &DtFrac))
            DtFrac = (vfloat)0.9;
        pXML->UpLevel();
    }

    if (pXML->FindElement("Damping")) {
        float tmp;
        if (!pXML->FindLoadElement("BondDampingZ", &tmp))
            SetBondDampZ(0.1);
        else
            SetBondDampZ(tmp); // BondDampingZ = 0.1;
        if (!pXML->FindLoadElement("ColDampingZ", &tmp))
            SetCollisionDampZ(1.0);
        else
            SetCollisionDampZ(tmp); // ColDampingZ = 1.0;
        if (!pXML->FindLoadElement("SlowDampingZ", &tmp))
            SetSlowDampZ(1.0);
        else
            SetSlowDampZ(tmp); // SlowDampingZ = 1.0;
        pXML->UpLevel();
    }

    if (pXML->FindElement("Collisions")) {
        if (!pXML->FindLoadElement("SelfColEnabled", &tmpBool))
            tmpBool = false;
        EnableFeature(VXSFEAT_COLLISIONS, tmpBool);
        // if (pXML->FindLoadElement("ColSystem", &tmpInt)) CurColSystem = (ColSystem)tmpInt; else CurColSystem = COL_SURFACE_HORIZON;
        // if (!pXML->FindLoadElement("CollisionHorizon", &CollisionHorizon)) CollisionHorizon = (vfloat)2.0;
        pXML->UpLevel();
    }

    if (pXML->FindElement("Features")) {
        if (!pXML->FindLoadElement("MaxVelLimitEnabled", &tmpBool))
            tmpBool = false;
        EnableFeature(VXSFEAT_MAX_VELOCITY, false); // EnableFeature(VXSFEAT_MAX_VELOCITY, tmpBool);
        // if (!pXML->FindLoadElement("MaxVoxVelLimit", &MaxVoxVelLimit)) MaxVoxVelLimit = (vfloat)0.1;
        if (!pXML->FindLoadElement("BlendingEnabled", &tmpBool))
            tmpBool = false;
        EnableFeature(VXSFEAT_BLENDING, tmpBool);

        if (pXML->FindLoadElement("MixRadius", &tmpVFloat))
            MixRadius = Vec3D<>(tmpVFloat, tmpVFloat, tmpVFloat); // look for legacy first
        else {
            if (!pXML->FindLoadElement("XMixRadius", &MixRadius.x))
                MixRadius.x = 0;
            if (!pXML->FindLoadElement("YMixRadius", &MixRadius.y))
                MixRadius.y = 0;
            if (!pXML->FindLoadElement("ZMixRadius", &MixRadius.z))
                MixRadius.z = 0;
        }

        if (pXML->FindLoadElement("BlendModel", &tmpInt))
            BlendModel = (MatBlendModel)tmpInt;
        else
            BlendModel = MB_LINEAR;
        if (!pXML->FindLoadElement("PolyExp", &PolyExp))
            PolyExp = 1.0;

        if (!pXML->FindLoadElement("FluidDampEnabled", &tmpBool))
            tmpBool = false; // do nothing for now...
        if (!pXML->FindLoadElement("VolumeEffectsEnabled", &tmpBool))
            tmpBool = false;
        EnableFeature(VXSFEAT_VOLUME_EFFECTS, tmpBool);
        if (!pXML->FindLoadElement("EnforceLatticeEnabled", &tmpBool))
            tmpBool = false; // do nothing for now...
        pXML->UpLevel();
    }

    if (pXML->FindElement("StopCondition")) {
        if (pXML->FindLoadElement("StopConditionType", &tmpInt))
            SetStopConditionType((StopCondition)tmpInt);
        else
            SetStopConditionType();
        if (pXML->FindLoadElement("StopConditionValue", &tmpVFloat))
            SetStopConditionValue(tmpVFloat);
        else
            SetStopConditionValue();
        pXML->UpLevel();
    }

    if (pXML->FindElement("EquilibriumMode")) {
        if (!pXML->FindLoadElement("EquilibriumModeEnabled", &tmpBool))
            tmpBool = false;
        if (tmpBool && !IsFeatureEnabled(VXSFEAT_EQUILIBRIUM_MODE))
            EnableFeature(VXSFEAT_EQUILIBRIUM_MODE, true);
        // if (EquilibriumModeEnabled) EnableEquilibriumMode(true); //so it can set up energy history if necessary
        pXML->UpLevel();
    }

    // MeshAutoGenerated=true;
    if (pXML->FindElement("SurfMesh")) {
        if (pXML->FindElement("CMesh")) {
            if (!ImportSurfMesh)
                ImportSurfMesh = new CMesh;
            // MeshAutoGenerated=false;
            ImportSurfMesh->ReadXML(pXML);
            pXML->UpLevel();
        }
        pXML->UpLevel();
    }

    return true; // ReadAdditionalSimXML(pXML, RetMessage);
}

void CVX_Sim::ClearAll(void) // Reset all initialized variables
{
    Initalized = false;
    LocalVXC.ClearMatter();

    // This should be all the stuff set by "Import()"

    ClearHistories();

    dt = (vfloat)0.0; // calculated per-step
    CurTime = (vfloat)0.0;
    CurStepCount = 0;

    SS.Clear();
    IniCM = Vec3D<>(0, 0, 0);

    delete ImportSurfMesh;
    ImportSurfMesh = NULL;
}

void CVX_Sim::EnableFeature(const int VXSFEAT, bool Enabled) {
    if (Enabled)
        CurSimFeatures |= VXSFEAT;
    else
        CurSimFeatures &= ~VXSFEAT;

    switch (VXSFEAT) {
    case VXSFEAT_GRAVITY:
        if (Enabled)
            Vx.setGravity(1.0);
        else
            Vx.setGravity(0.0);
        break;
    case VXSFEAT_FLOOR:
        Vx.enableFloor(Enabled);
        break;
    case VXSFEAT_COLLISIONS:
        Vx.enableCollisions(Enabled);
        break;
    case VXSFEAT_EQUILIBRIUM_MODE:
        EnableEquilibriumMode(Enabled);
        break;
    case VXSFEAT_TEMPERATURE:
        if (pEnv)
            pEnv->EnableTemp(Enabled);
        UpdateMatTemps();
        break;
    case VXSFEAT_TEMPERATURE_VARY:
        if (pEnv)
            pEnv->EnableTempVary(Enabled);
        break;
    case VXSFEAT_VOLUME_EFFECTS:
        EnableVolumeEffects(Enabled);
        break;
    }
}

bool CVX_Sim::IsFeatureEnabled(const int VXSFEAT) {
    switch (VXSFEAT) {
    case VXSFEAT_GRAVITY:
        return Vx.gravity() != 0.0f;
        break;
    case VXSFEAT_FLOOR:
        return Vx.isFloorEnabled();
        break;
    case VXSFEAT_COLLISIONS:
        return Vx.isCollisionsEnabled();
        break; //|| stuffity stuff
    }
    return (CurSimFeatures & VXSFEAT);
}

void CVX_Sim::CopyMat(CVXC_Material *pOld, CVX_Material *pNew) // copies parameters from pOld to pNew
{
    pNew->matid = pOld->matid;
    pNew->isPaceMaker = (bool)pOld->isPaceMaker;
    pNew->PaceMakerPeriod = pOld->PaceMakerPeriod;
    pNew->isElectricalActive = (bool)pOld->isElectricalActive;
    pNew->signalValueDecay = pOld->signalValueDecay;
    pNew->signalTimeDelay = pOld->signalTimeDelay;
    pNew->inactivePeriod = pOld->inactivePeriod;
    pNew->isMeasured = pOld->isMeasured;

    pNew->RemoveFromSimulationAfterThisManySeconds = pOld->RemoveFromSimulationAfterThisManySeconds;
    pNew->TurnOnThermalExpansionAfterThisManySeconds = pOld->TurnOnThermalExpansionAfterThisManySeconds;
    pNew->TurnOnCiliaAfterThisManySeconds = pOld->TurnOnCiliaAfterThisManySeconds;
	
    pNew->isTarget = (bool)pOld->isTarget;
    pNew->fixed = (bool)pOld->Fixed;
    pNew->sticky = (bool)pOld->sticky;
    pNew->Cilia = pOld->Cilia;

    pNew->setName(pOld->GetName().c_str());
    pNew->setColor(pOld->GetRedi(), pOld->GetGreeni(), pOld->GetBluei(), pOld->GetAlphai());
    switch (pOld->GetMatModel()) {
    case MDL_LINEAR:
        pNew->setModelLinear(pOld->GetElasticMod());
        break;
    case MDL_LINEAR_FAIL:
        pNew->setModelLinear(pOld->GetElasticMod(), pOld->GetFailStress());
        break;
    case MDL_BILINEAR:
        pNew->setModelBilinear(pOld->GetElasticMod(), pOld->GetPlasticMod(), pOld->GetYieldStress(), pOld->GetFailStress());
        break;
    case MDL_DATA: {
        std::vector<float> tmpStress, tmpStrain;
        int numPts = pOld->GetDataPointCount();
        for (int i = 0; i < numPts; i++) {
            tmpStress.push_back(pOld->GetStressData(i));
            tmpStrain.push_back(pOld->GetStrainData(i));
        }
        pNew->setModel(numPts, &(tmpStrain[0]), &(tmpStress[0]));
        break;
    }
    }

    pNew->setPoissonsRatio(pOld->GetPoissonsRatio());
    pNew->setDensity(pOld->GetDensity());
    pNew->setCte(pOld->GetCTE());
    pNew->setStaticFriction(pOld->GetuStatic());
    pNew->setKineticFriction(pOld->GetuDynamic());
    pNew->setGlobalDamping(GetSlowDampZ());
    pNew->setInternalDamping(GetBondDampZ());
    pNew->setCollisionDamping(GetCollisionDampZ());
}

void CVX_Sim::EnableEquilibriumMode(bool Enabled) {
    if (Enabled) {
        MemBondDampZ = BondDampingZ;
        MemSlowDampingZ = SlowDampingZ;

        for (int i = 0; i < Vx.materialCount(); i++) { // set for each material
            CVX_Material *pMat = Vx.material(i);
            if (i == 0) {
                MemBondDampZ = pMat->internalDamping();
                MemSlowDampingZ = pMat->globalDamping();
            }
            pMat->setInternalDamping(1.0);
            pMat->setGlobalDamping(0.0);
        }
    } else {
        for (int i = 0; i < Vx.materialCount(); i++) { // set for each material
            CVX_Material *pMat = Vx.material(i);
            pMat->setInternalDamping(MemBondDampZ);
            pMat->setGlobalDamping(MemSlowDampingZ);
        }
    }
}

void CVX_Sim::EnableVolumeEffects(bool Enabled) {
    if (Vx.materialCount() != muMemory.size())
        return;

    if (Enabled) {
        for (int i = 0; i < Vx.materialCount(); i++)
            Vx.material(i)->setPoissonsRatio(muMemory[i]);
    } else {
        for (int i = 0; i < Vx.materialCount(); i++)
            Vx.material(i)->setPoissonsRatio(0);
    }

    OptimalDt = Vx.recommendedTimeStep();
}

bool CVX_Sim::IsVolumeEffectsEnabled() {
    for (int i = 0; i < Vx.materialCount(); i++)
        Vx.material(i)->setPoissonsRatio(muMemory[i]);

    return false;
}

void CVX_Sim::UpdateMatTemps(void) // updates expansions for each material
{
    for (int iz = Vx.indexMinZ(); iz <= Vx.indexMaxZ(); iz++) {
        for (int iy = Vx.indexMinY(); iy <= Vx.indexMaxY(); iy++) {
            for (int ix = Vx.indexMinX(); ix <= Vx.indexMaxX(); ix++) {
                CVX_Voxel *pV = Vx.voxel(ix, iy, iz);
                float thisTemp = 0;
                if (pV != NULL) {
                    if (IsFeatureEnabled(VXSFEAT_TEMPERATURE))
                        pV->setTemperature(pEnv->UpdateCurTemp(CurTime, &LocalVXC) - pEnv->GetTempBase()); // pEnv->GetTempAmplitude());
                    else
                        pV->setTemperature(0);
                }
            }
        }
    }
}

void CVX_Sim::SetGravityAccel(float grav) { Vx.setGravity(-grav / 9.80665); }

float CVX_Sim::GetGravityAccel(void) { return -9.80665 * Vx.gravity(); }

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// below are not used in new system.
// void CVX_Sim::SaveVXAFile(std::string filename)
// {
// 	CXML_Rip XML;
// 	WriteVXA(&XML);
// 	XML.SaveFile(filename);
// }

// void CVX_Sim::WriteVXA(CXML_Rip* pXML)
// {
// 	pXML->DownLevel("VXA");
// 	pXML->SetElAttribute("Version", "1.1");
// 	WriteXML(pXML);
// 	pEnv->WriteXML(pXML);
// 	pEnv->pObj->WriteXML(pXML);
// 	pXML->UpLevel();
// }
// void CVX_Sim::WriteXML(CXML_Rip* pXML)
// {
// 	pXML->DownLevel("Simulator");
// 		pXML->DownLevel("Integration");
// 		pXML->Element("Integrator", 0); //0 = euler in older versions
// 		pXML->Element("DtFrac", DtFrac);
// 		pXML->UpLevel();

// 		pXML->DownLevel("Damping");
// 		pXML->Element("BondDampingZ", BondDampingZ); //BondDampingZ);
// 		pXML->Element("ColDampingZ", ColDampingZ); //ColDampingZ);
// 		pXML->Element("SlowDampingZ", SlowDampingZ); //SlowDampingZ);
// //		pXML->Element("BondDampingZ", Vx.material(0)->internalDamping()); //BondDampingZ);
// //		pXML->Element("ColDampingZ", Vx.material(0)->collisionDamping()); //ColDampingZ);
// //		pXML->Element("SlowDampingZ", Vx.material(0)->globalDamping()); //SlowDampingZ);
// 		pXML->UpLevel();

// 		pXML->DownLevel("Collisions");
// 		pXML->Element("SelfColEnabled", IsFeatureEnabled(VXSFEAT_COLLISIONS));
// 		pXML->Element("ColSystem", COL_SURFACE_HORIZON);
// 		pXML->Element("CollisionHorizon", 3.0);
// 		pXML->UpLevel();

// 		pXML->DownLevel("Features");
// //		pXML->Element("MaxVelLimitEnabled", IsFeatureEnabled(VXSFEAT_MAX_VELOCITY));
// //		pXML->Element("MaxVoxVelLimit", MaxVoxVelLimit);
// 		pXML->Element("BlendingEnabled", IsFeatureEnabled(VXSFEAT_BLENDING));
// 		pXML->Element("XMixRadius", MixRadius.x);
// 		pXML->Element("YMixRadius", MixRadius.y);
// 		pXML->Element("ZMixRadius", MixRadius.z);
// 		pXML->Element("BlendModel", BlendModel);
// 		pXML->Element("PolyExp", PolyExp);
// //		pXML->Element("FluidDampEnabled", IsFeatureEnabled(VXSFEAT_));
// 		pXML->Element("VolumeEffectsEnabled", IsFeatureEnabled(VXSFEAT_VOLUME_EFFECTS));
// //		pXML->Element("EnforceLatticeEnabled", IsFeatureEnabled(VXSFEAT_));
// 		pXML->UpLevel();

// 		pXML->DownLevel("StopCondition");
// 		pXML->Element("StopConditionType", (int)StopConditionType);
// 		pXML->Element("StopConditionValue", StopConditionValue);
// 		pXML->UpLevel();

// 		pXML->DownLevel("EquilibriumMode");
// 		pXML->Element("EquilibriumModeEnabled", IsFeatureEnabled(VXSFEAT_EQUILIBRIUM_MODE));
// 	//	pXML->Element("StopConditionValue", StopConditionValue);
// 		pXML->UpLevel();

// 		if (ImportSurfMesh){
// 			pXML->DownLevel("SurfMesh");
// 			ImportSurfMesh->WriteXML(pXML, true);
// 			pXML->UpLevel();
// 		}

// //		WriteAdditionalSimXML(pXML);
// 	pXML->UpLevel();

// }
// CVX_Sim& CVX_Sim::operator=(const CVX_Sim& rSim) //overload "="
// {
// 	//TODO: set everything sensible equal.

// 	return *this;
// }

// void CVX_Sim::ResetSimulation(void)
// {
// 	Vx.resetTime();

// 	dt = (vfloat)0.0; //calculated per-step
// 	CurTime = (vfloat)0.0;
// 	CurStepCount = 0;

// }

// /*! Given the current state of the simulation (Voxel positions and velocities) and information about the current environment, advances
// the simulation by the maximum stable timestep. The integration scheme denoted by the CurIntegrator member variable is used. Calculates
// some relevant system statistics such as maximum displacements and velocities and total force. Returns true if the time step was
// successful, false otherwise.
// @param[out] pRetMessage Pointer to an initialized string. Messages generated in this function will be appended to the string.
// */
// bool CVX_Sim::TimeStep(std::string* pRetMessage)
// {
// 	if(IsFeatureEnabled(VXSFEAT_TEMPERATURE_VARY)) UpdateMatTemps(); //updates the temperatures

// 	if (IsFeatureEnabled(VXSFEAT_VOLUME_EFFECTS)) dt = DtFrac*Vx.recommendedTimeStep();
// 	else dt = DtFrac*OptimalDt;

// 	CurTime += dt; //keep track of time!
// 	CurStepCount++; //keep track of current step...

// 		//update information to calculate
// 	switch (GetStopConditionType()){ //may need to calculate certain items depending on stop condition
// 		case SC_CONST_MAXENERGY: StatToCalc |= CALCSTAT_KINE; StatToCalc |= CALCSTAT_STRAINE; break;
// 		case SC_MIN_KE: StatToCalc |= CALCSTAT_KINE; break;
// 		case SC_MIN_MAXMOVE: StatToCalc |= CALCSTAT_VEL; break;
// 	}
// 	if (IsFeatureEnabled(VXSFEAT_EQUILIBRIUM_MODE)) StatToCalc |= CALCSTAT_KINE;

// 	if (!Vx.doTimeStep(dt)){
// 		if (pRetMessage) *pRetMessage = "Simulation Diverged. Please reduce forces or accelerations.\n";
// 		return false;
// 	}

// 	if (IsFeatureEnabled(VXSFEAT_EQUILIBRIUM_MODE) && KineticEDecreasing()){
// 		ZeroAllMotion(); MotionZeroed = true;}
// 	else MotionZeroed = false;
// 	UpdateStats(pRetMessage);

// 	return true;
// }

// void CVX_Sim::ZeroAllMotion(void)
// {
// 	const std::vector<CVX_Voxel*>* pVoxList = Vx.voxelList();
// 	for (std::vector<CVX_Voxel*>::const_iterator it = pVoxList->begin(); it != pVoxList->end(); it++){
// 		(*it)->haltMotion();
// 	}

// }

// bool CVX_Sim::StopConditionMet(void) //have we met the stop condition yet?
// {
// 	if (CurStepCount<2*HISTORY_SIZE) return false;

// 	int numJump; //how many timesteps to look back in order to have 10 data points within the history length
// 	vfloat fNumVoxInv;
// 	if (StopConditionType==SC_CONST_MAXENERGY || StopConditionType==SC_MIN_KE || StopConditionType==SC_MIN_MAXMOVE){
// 		fNumVoxInv = 1.0/(float)Vx.voxelCount();
// 		numJump = HISTORY_SIZE/10;
// 	}

// 	switch(StopConditionType){
// 		case SC_NONE: return false;
// 		case SC_MAX_TIME_STEPS: return (CurStepCount>(int)(StopConditionValue+0.5))?true:false;
// 		case SC_MAX_SIM_TIME: return CurTime>StopConditionValue?true:false;
// 		case SC_TEMP_CYCLES:  return CurTime>pEnv->GetTempPeriod()*StopConditionValue?true:false;
// 		case SC_CONST_MAXENERGY:{
// 			vfloat IniTotVal = TotEHistory[0];
// 			for (int i=numJump; i<HISTORY_SIZE; i+=numJump){
// 				if (TotEHistory[i] == -1) return false;
// 				if (abs(TotEHistory[i]-IniTotVal)*fNumVoxInv > 0.001*StopConditionValue) return false;
// 			}
// 			return true;
// 		  }
// 		case SC_MIN_KE:{
// 			for (int i=0; i<HISTORY_SIZE; i+=numJump){
// 				if (KinEHistory[i] == -1) return false;
// 				if (KinEHistory[i]*fNumVoxInv > 0.001*StopConditionValue) return false;
// 			}
// 			return true;
// 		  }
// 		case SC_MIN_MAXMOVE:{
// 			for (int i=0; i<HISTORY_SIZE; i+=numJump){
// 				if (MaxMoveHistory[i] == -1) return false;
// 				if (MaxMoveHistory[i] > 0.001*StopConditionValue) return false;
// 			}
// 			return true;
// 		}

// 		default: return false;
// 	}
// }

// bool CVX_Sim::UpdateStats(std::string* pRetMessage) //updates simulation state (SS)
// {
// 	//if (SelfColEnabled) StatToCalc |= CALCSTAT_VEL; //always need velocities if self collisition is enabled
// 	if (IsFeatureEnabled(VXSFEAT_COLLISIONS)) StatToCalc |= CALCSTAT_VEL; //always need velocities if self collisition is enabled
// 	if (StatToCalc == CALCSTAT_NONE) return true;
// 	bool CCom=StatToCalc&CALCSTAT_COM, CDisp=StatToCalc&CALCSTAT_DISP, CVel=StatToCalc & CALCSTAT_VEL, CKinE=StatToCalc&CALCSTAT_KINE,
// CStrE=StatToCalc&CALCSTAT_STRAINE, CEStrn=StatToCalc&CALCSTAT_ENGSTRAIN, CEStrs=StatToCalc&CALCSTAT_ENGSTRESS,
// CPressure=StatToCalc&CALCSTAT_PRESSURE;

// 	if (CCom) SS.CurCM = GetCM(); //calculate center of mass

// 	//update the overall statisics (can't do this within threaded loops and be safe without mutexes...
// 	vfloat tmpMaxVoxDisp2 = 0, tmpMaxVoxVel2 = 0, tmpMaxVoxKineticE = 0, tmpMaxVoxStrainE = 0, tmpMaxPressure = -FLT_MAX, tmpMinPressure
// = FLT_MAX; 	vfloat tmpMaxBondStrain=0, tmpMaxBondStress=0, tmpTotalObjKineticE = 0, tmpTotalObjStrainE=0; 	Vec3D<> tmpTotalObjDisp(0,0,0);

// 	if (CDisp || CVel || CKinE || CPressure){
// 		int nVox = Vx.voxelCount();// NumVox();
// 		for (int i=0; i<nVox; i++){ //for each voxel
// 			const CVX_Voxel* it = Vx.voxel(i); //pointer to this voxel

// 			if (CDisp) { //Displacements
// 				tmpTotalObjDisp += it->velocity().Abs()*dt; //keep track of displacements on global object
// 				const float ThisMaxVoxDisp2 = it->displacement().Length2();
// 				if (ThisMaxVoxDisp2 > tmpMaxVoxDisp2) tmpMaxVoxDisp2 = ThisMaxVoxDisp2;
// 			}

// 			if (CVel) { //Velocities
// 				const vfloat ThisMaxVoxVel2 = it->velocity().Length2(); //it->GetCurVel().Length2();
// 				if (ThisMaxVoxVel2 > tmpMaxVoxVel2) tmpMaxVoxVel2 = ThisMaxVoxVel2;
// 			}
// 			if (CKinE) { // kinetic energy
// 				const vfloat ThisMaxKineticE = it->kineticEnergy(); // kineticEnergy(); // it->GetCurKineticE();
// 				if (ThisMaxKineticE > tmpMaxVoxKineticE) tmpMaxVoxKineticE = ThisMaxKineticE;
// 				tmpTotalObjKineticE += ThisMaxKineticE; //keep track of total kinetic energy
// 			}
// 			if (CPressure){
// 				const vfloat ThisPressure = it->pressure();
// 				if (ThisPressure > tmpMaxPressure) tmpMaxPressure = ThisPressure;
// 				if (ThisPressure < tmpMinPressure) tmpMinPressure = ThisPressure;

// 			}
// 		}

// 		if (CDisp){ //Update SimState (SS)
// 			tmpTotalObjDisp /= nVox;
// 			SS.TotalObjDisp = tmpTotalObjDisp;
// 			SS.NormObjDisp = tmpTotalObjDisp.Length();
// 			SS.MaxVoxDisp = sqrt(tmpMaxVoxDisp2);
// 		}

// 		if (CVel) SS.MaxVoxVel = sqrt(tmpMaxVoxVel2);
// 		if (CKinE) {
// 			SS.MaxVoxKinE = tmpMaxVoxKineticE;
// 			SS.TotalObjKineticE = tmpTotalObjKineticE;
// 		}
// 		if (CPressure){
// 			SS.MaxPressure = tmpMaxPressure;
// 			SS.MinPressure = tmpMinPressure;
// 		}
// 	}

// 	if (CStrE || CEStrn || CEStrs){
// 		int nLink = Vx.linkCount();
// 		for (int i=0; i<nLink; i++){ //for each voxel
// 			CVX_Link* it = Vx.link(i);
// 	//	for (std::vector<CVXS_BondInternal>::iterator it = BondArrayInternal.begin(); it != BondArrayInternal.end(); it++){
// 			if (CStrE){
// 				const vfloat ThisMaxStrainE =  it->strainEnergy(); // GetStrainEnergy();
// 				if (ThisMaxStrainE > tmpMaxVoxStrainE) tmpMaxVoxStrainE = ThisMaxStrainE;
// 				tmpTotalObjStrainE += ThisMaxStrainE;
// 			}

// 			if (CEStrn && it->axialStrain() > tmpMaxBondStrain) tmpMaxBondStrain = it->axialStrain(); //shouldn't these pull from
// bonds? would make more sense... 			if (CEStrs && it->axialStress() > tmpMaxBondStress) tmpMaxBondStress = it->axialStress();

// //			if (CEStrn && it->GetEngStrain() > tmpMaxBondStrain) tmpMaxBondStrain = it->GetEngStrain(); //shouldn't these pull from
// bonds? would make more sense...
// //			if (CEStrs && it->GetEngStress() > tmpMaxBondStress) tmpMaxBondStress = it->GetEngStress();
// 		}

// 		//Updata SimState (SS)
// 		if (CStrE){
// 			SS.MaxBondStrainE = tmpMaxVoxStrainE;
// 			SS.TotalObjStrainE = tmpTotalObjStrainE;
// 		}

// 		if (CEStrn) SS.MaxBondStrain = tmpMaxBondStrain;
// 		if (CEStrs) SS.MaxBondStress = tmpMaxBondStress;

// 	}

// 	//update histories
// 	MaxMoveHistory.push_front(CVel ? SS.MaxVoxVel*dt : -1.0); MaxMoveHistory.pop_back();
// 	KinEHistory.push_front(CKinE ? SS.TotalObjKineticE : -1.0); KinEHistory.pop_back();
// 	TotEHistory.push_front((CStrE && CKinE) ? SS.TotalObjKineticE + SS.TotalObjStrainE : -1.0); TotEHistory.pop_back();

// 	return true;
// }
// void CVX_Sim::UpdateMuMemory(void) //updates array for turning poissons ratio on and off.
// {
// 	bool wasVolEnabled = IsFeatureEnabled(VXSFEAT_VOLUME_EFFECTS);
// 	if (!wasVolEnabled) EnableFeature(VXSFEAT_VOLUME_EFFECTS);

// 	muMemory.clear();
// 	for (int i=0; i<Vx.materialCount(); i++){ //for each material
// 		muMemory.push_back(Vx.material(i)->poissonsRatio()); //remember for toggleing volume effects on and off.
// 	}

// 	if (!wasVolEnabled) EnableFeature(VXSFEAT_VOLUME_EFFECTS, false);
// }

// Vec3D<> CVX_Sim::GetCM(void)
// {
// 	vfloat TotalMass = 0;
// 	Vec3D<> Sum(0,0,0);
// 	int nVox = Vx.voxelCount();
// 	for (int i=0; i<nVox; i++){
// 		CVX_Voxel* it = Vx.voxel(i);
// 		vfloat ThisMass = it->material()->mass();
// 		Sum += it->position()*ThisMass;
// 		TotalMass += ThisMass;
// 	}

// 	return Sum/TotalMass;
// }

// int CVX_Sim::GetNumTouchingFloor()
// {
// 	int NumTouching = 0;

// 	int LocNumVox = Vx.voxelCount();
// 	for (int i=0; i<LocNumVox; i++){
// 		if (Vx.voxel(i)->floorPenetration() > 0) NumTouching++;
// 	}
// 	return NumTouching;
// }

// bool CVX_Sim::KineticEDecreasing(void)
// {
// 	 if (KinEHistory[0]+KinEHistory[1]+KinEHistory[2] < KinEHistory[3]+KinEHistory[4]+KinEHistory[5] && !(KinEHistory[0] == 0 ||
// KinEHistory[1] == 0 || KinEHistory[2] == 0 || KinEHistory[3] == 0 || KinEHistory[4] == 0 || KinEHistory[5] == 0)) return true; 	 else
// return false;
// }

// Vec3D<> CVX_Sim::GetSumForce(CVX_FRegion* pRegion)
// {
// 	return Vec3D<>(0,0,0);
// }

// vfloat CVX_Sim::GetSumForceDir(CVX_FRegion* pRegion)
// {
// 	//right now only fixed regions... (forced regions should be zero!)
// 	//get force only in dircetion of pull!
// 	Vec3D<> Res = GetSumForce(pRegion);

// 	Vec3D<> Dir = pRegion->Displace;
// 	if (Dir.Length2() == 0) return Res.Length(); //return magnitude of no direction...
// 	else {
// 		Dir.Normalize();
// 		return Res.Dot(Dir);
// 	}

// }

// Vec3D<> CVX_Sim::GetAvgDisplace(CVX_FRegion* pRegion) //returns the average displacement in x,y,z of a region
// {
// 	return Vec3D<>(0,0,0);
// }
// int CVX_Sim::NumYielded(void)
// {
// 	int NumYieldRet = 0;
// 	int NumBondNow = Vx.linkCount();
// 	for (int i=0; i<NumBondNow; i++)
// 		if (Vx.link(i)->isYielded()) NumYieldRet++;

// 	return NumYieldRet;
// }

// int CVX_Sim::NumBroken(void)
// {
// 	int NumBrokenRet = 0;
// 	int NumBondNow = Vx.linkCount();
// 	for (int i=0; i<NumBondNow; i++)
// 		if (Vx.link(i)->isFailed()) NumBrokenRet++;

// 	return NumBrokenRet;

// }

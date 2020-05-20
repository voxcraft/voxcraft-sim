/*******************************************************************************
Copyright (c) 2010, Jonathan Hiller (Cornell University)
If used in publication cite "J. Hiller and H. Lipson "Dynamic Simulation of Soft Heterogeneous Objects" In press. (2011)"

This file is part of Voxelyze.
Voxelyze is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Voxelyze is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
See <http://www.opensource.org/licenses/lgpl-3.0.html> for license details.
*******************************************************************************/

#include "VXS_Voxel.h"
#include "VXS_Bond.h"
#include "VX_Sim.h"

CVXS_Voxel::CVXS_Voxel(CVX_Sim* pSimIn, int SIndexIn, int XIndexIn, int MatIndexIn, Vec3D<>& NominalPositionIn, vfloat OriginalScaleIn) : CVXt_Voxel(pSimIn, SIndexIn, XIndexIn, MatIndexIn, NominalPositionIn, OriginalScaleIn)
{
	ResetVoxel(); //sets state variables to zero


	ExternalInputScale = 1.0;

//	CornerPosCur = Vec3D<>(OriginalScaleIn/2, OriginalScaleIn/2, OriginalScaleIn/2);
//	CornerNegCur = Vec3D<>(-OriginalScaleIn/2, -OriginalScaleIn/2, -OriginalScaleIn/2);

	SetColor(0,0,0,1);


}

CVXS_Voxel::~CVXS_Voxel(void)
{
}

CVXS_Voxel& CVXS_Voxel::operator=(const CVXS_Voxel& VIn)
{
	CVXt_Voxel::operator=(VIn);
	
	ExternalInputScale=VIn.ExternalInputScale;
	
	InputForce = Vec3D<>(0,0,0);

	Pos = VIn.Pos;
	LinMom = VIn.LinMom;
	Angle = VIn.Angle;
	AngMom = VIn.AngMom;
	Vel = VIn.Vel;
	KineticEnergy = VIn.KineticEnergy;
	AngVel = VIn.AngVel;
	Pressure = VIn.Pressure;
	Scale=VIn.Scale;

	StaticFricFlag = VIn.StaticFricFlag;
	VYielded = VIn.VYielded;
	VBroken = VIn.VBroken;

	ColBondInds = VIn.ColBondInds;
	UpdateColBondPointers();

	m_Red = VIn.m_Red;
	m_Green = VIn.m_Green;
	m_Blue = VIn.m_Blue;
	m_Trans = VIn.m_Trans;

//	SizeCurrent = VIn.SizeCurrent;
	CornerPosCur = VIn.CornerPosCur;
	CornerNegCur = VIn.CornerNegCur;

	ForceCurrent = VIn.ForceCurrent;

	poissonsstrain = VIn.poissonsstrain;

	return *this;
}

void CVXS_Voxel::ResetVoxel(void) //resets this voxel to its default (imported) state.
{
	LinMom = Vec3D<double>(0,0,0);
	Angle = Quat3D<double>(1.0, 0, 0, 0);
	AngMom = Vec3D<double>(0,0,0);
//	Scale = 0;
	Vel = Vec3D<>(0,0,0);
	KineticEnergy = 0;
	AngVel = Vec3D<>(0,0,0);
	Pressure=0;


	Pos = GetNominalPosition(); //only position and size need to be set
	Scale = GetNominalSize();


	InputForce = Vec3D<>(0,0,0); //?


	StaticFricFlag = false;
	VYielded = false;
	VBroken = false;

//	SizeCurrent = Vec3D<>(Scale, Scale, Scale);
	CornerPosCur = Vec3D<>(Scale/2, Scale/2, Scale/2);
	CornerNegCur = Vec3D<>(-Scale/2, -Scale/2, -Scale/2);

	ForceCurrent = Vec3D<>(0,0,0);
//	StrainPosDirsCur = Vec3D<>(0,0,0);
//	StrainNegDirsCur = Vec3D<>(0,0,0);

	poissonsstrain=Vec3D<>(0,0,0);
}


bool CVXS_Voxel::LinkColBond(int CBondIndex) //simulation bond index...
{
	//if (!pSim || CBondIndex >= pSim->BondArrayCollision.size()) return false;

	//ColBondInds.push_back(CBondIndex);
	//ColBondPointers.push_back(&(pSim->BondArrayCollision[CBondIndex]));

	return true;
}

void CVXS_Voxel::UnlinkColBonds(void)
{
	ColBondInds.clear();
	ColBondPointers.clear();
}


void CVXS_Voxel::UpdateColBondPointers() //updates all links (pointers) to bonds according top current p_Sim
{
	//int NumColBonds = ColBondInds.size();
	//if (NumColBonds == 0) return;
	//ColBondPointers.resize(NumColBonds);
	//for (int i=0; i<NumColBonds; i++){
	//	ColBondPointers[i] = &(pSim->BondArrayCollision[ColBondInds[i]]);
	//}
}





//http://klas-physics.googlecode.com/svn/trunk/src/general/Integrator.cpp (reference)
void CVXS_Voxel::EulerStep()
{
	double dt = pSim->dt;
	//bool EqMode = p_Sim->IsEquilibriumEnabled();
	if (IS_ALL_FIXED(DofFixed) & !pSim->IsFeatureEnabled(VXSFEAT_VOLUME_EFFECTS)){ //if fixed, just update the position and forces acting on it (for correct simulation-wide summing
		LinMom = Vec3D<double>(0,0,0);
		Pos = NominalPosition + ExternalInputScale*ExternalDisp;
		AngMom = Vec3D<double>(0,0,0);
		//Angle.FromRotationVector(Vec3D<double>(ExternalInputScale*ExternalTDisp));
		Angle = Quat3D<double>(Vec3D<double>(ExternalInputScale*ExternalTDisp));

	}
	else {
		Vec3D<> ForceTot = CalcTotalForce(); //TotVoxForce;

		//DISPLACEMENT
		LinMom = LinMom + ForceTot*dt;
		Vec3D<double> Disp(LinMom*(dt*_massInv)); //vector of what the voxel moves

		if (pSim->IsFeatureEnabled(VXSFEAT_FLOOR) && GetCurGroundPenetration() > 0){
			vfloat work = ForceTot.x*Disp.x + ForceTot.y*Disp.y; //F dot disp
			if(KineticEnergy + work < 0){ //change of direction
				StaticFricFlag = true;
			}
		}
		if (StaticFricFlag){
			LinMom.x = LinMom.y = 0;
			Disp.x = Disp.y = 0;
		}

//		vfloat deltaKE = force dot displacement

		//if(pSim->IsFeatureEnabled(VXSFEAT_MAX_VELOCITY)){ //check to make sure we're not going over the speed limit!
		//	vfloat DispMag = Disp.Length();
		//	vfloat MaxDisp = pSim->GetMaxVoxVelLimit()*NominalSize; // p_Sim->pEnv->pObj->GetLatticeDim();
		//	if (DispMag>MaxDisp) Disp *= (MaxDisp/DispMag);
		//}
		Pos += Disp; //update position (source of noise in float mode???

		if (IS_FIXED(DOF_X, DofFixed)){Pos.x = NominalPosition.x + ExternalInputScale*ExternalDisp.x; LinMom.x = 0;}
		if (IS_FIXED(DOF_Y, DofFixed)){Pos.y = NominalPosition.y + ExternalInputScale*ExternalDisp.y; LinMom.y = 0;}
		if (IS_FIXED(DOF_Z, DofFixed)){Pos.z = NominalPosition.z + ExternalInputScale*ExternalDisp.z; LinMom.z = 0;}

		//ANGLE
		Vec3D<> TotVoxMoment = CalcTotalMoment(); //debug
		
		AngMom = AngMom + TotVoxMoment*dt;

//		if (pSim->IsFeatureEnabled(VXSFEAT_IGNORE_ANG_DAMP)){
//			//AngMom *= (1-pSim->GetSlowDampZ());
//			AngMom /= 1.01;
//		}
//		else {
			if (pSim->IsFeatureEnabled(VXSFEAT_VOLUME_EFFECTS)) AngMom /= 1.01; //TODO: remove angmom altogehter???
			else {
				vfloat AngMomFact = (1 - 10*pSim->GetSlowDampZ() * _inertiaInv *_2xSqIxExSxSxS*dt);
				AngMom *= AngMomFact; 
			}
//		}

		//convert Angular velocity to quaternion form ("Spin")
		Vec3D<double> dSAngVel(AngMom * _inertiaInv);

//old version
//		Quat3D<double> Spin = 0.5 * Quat3D<double>(0, dSAngVel.x, dSAngVel.y, dSAngVel.z) * Angle; //current "angular velocity"
//		Angle += Quat3D<double>(Spin*dt); //see above
//		Angle.NormalizeFast(); //Through profiling, quicker to normalize every time than check to see if needed then do it...

//from http://physicsforgames.blogspot.com/2010/02/quaternions.html (superior!)
		Quat3D<double> deltaQ;
		Vec3D<double> theta = 0.5*dt*dSAngVel;
		double thetaMag2 = theta.Length2();
		double s;
		if (thetaMag2*thetaMag2 / 24.0 < DBL_EPSILON){
			deltaQ.w = 1.0 - thetaMag2 / 2.0;
			s = 1.0 - thetaMag2 / 6.0;
		 }
		 else {
			  float thetaMag = sqrt(thetaMag2);
			  deltaQ.w = cos(thetaMag);
			  s = sin(thetaMag) / thetaMag;
		 }
		 deltaQ.x = theta.x * s;
		 deltaQ.y = theta.y * s;
		 deltaQ.z = theta.z * s;
		 Angle = deltaQ*Angle;


	//	TODO: Only constrain fixed angles if one is non-zero! (support symmetry boundary conditions while still only doing this calculation) (only works if all angles are constrained for now...)
		if (IS_FIXED(DOF_TX, DofFixed) && IS_FIXED(DOF_TY, DofFixed) && IS_FIXED(DOF_TZ, DofFixed)){
//			Angle.FromRotationVector(Vec3D<double>(ExternalInputScale*ExternalTDisp));
			Angle = Quat3D<double>(Vec3D<double>(ExternalInputScale*ExternalTDisp));
			AngMom = Vec3D<>(0,0,0);
		}
	}

	//SCALE
	//	ScaleMom = ScaleMom + CalcTotalScaleForce()*p_Sim->dt;
	vfloat TempFact = 1.0;
	if(pSim->IsFeatureEnabled(VXSFEAT_TEMPERATURE)){ 
		//TempFact = (1+(p_Sim->pEnv->CurTemp-p_Sim->pEnv->TempBase)*GetCTE()); //LocalVXC.GetBaseMat(VoxArray[i].MatIndex)->GetCTE());
		
//		vfloat ThisTemp = p_Sim->pEnv->pObj->GetBaseMat(GetMaterial())->GetCurMatTemp();
//		vfloat ThisTemp = pSim->pEnv->pObj->GetBaseMat(GetMaterialIndex())->GetCurMatTemp();
		vfloat ThisTemp = _pMat->GetCurMatTemp();

		vfloat ThisCTE = GetCTE();
		vfloat TempBase =  pSim->pEnv->GetTempBase();

		TempFact = (1+(ThisTemp - TempBase)*ThisCTE);	//To allow selective temperature actuation for each different material
	}
	if (TempFact < MIN_TEMP_FACTOR) TempFact = MIN_TEMP_FACTOR;
	Scale = TempFact*NominalSize;

	//Recalculate secondary:
	AngVel = AngMom * _inertiaInv;
	Vel = LinMom * _massInv;
	if(pSim->StatToCalc & CALCSTAT_KINE) KineticEnergy = 0.5*Mass*Vel.Length2() + 0.5*Inertia*AngVel.Length2(); //1/2 m v^2
	if(pSim->StatToCalc & CALCSTAT_PRESSURE) {
//		vfloat VolumetricStrain = GetVoxelStrain(AXIS_X) + GetVoxelStrain(AXIS_Y) + GetVoxelStrain(AXIS_Z);
		vfloat VolumetricStrain = strain(false).x+strain(false).y+strain(false).z;
		Pressure = - _pMat->GetElasticMod()*VolumetricStrain/(3*(1-2*_pMat->GetPoissonsRatio())); //http://www.colorado.edu/engineering/CAS/courses.d/Structures.d/IAST.Lect05.d/IAST.Lect05.pdf
	}

	updatePoissonsstrain();

}


void CVXS_Voxel::SetColor(float r, float g, float b, float a)
{
	m_Red = r;
	m_Green = g;
	m_Blue = b;
	m_Trans = a;
}	

//void CVXS_Voxel::SetStrainDir(BondDir Bond, vfloat StrainIn)
//{
//	switch (Bond){
//	case BD_PX: StrainPosDirsCur.x = StrainIn; break;
//	case BD_PY: StrainPosDirsCur.y = StrainIn; break;
//	case BD_PZ: StrainPosDirsCur.z = StrainIn; break;
//	case BD_NX: StrainNegDirsCur.x = StrainIn; break;
//	case BD_NY: StrainNegDirsCur.y = StrainIn; break;
//	case BD_NZ: StrainNegDirsCur.z = StrainIn; break;
//	}
//}
//
//vfloat CVXS_Voxel::GetVoxelStrain(Axis DesiredAxis)
//{
//	bool pd, nd; //positive and negative directions
//	switch (DesiredAxis){
//		case AXIS_X:
//			pd = InternalBondPointers[BD_PX]!=NULL, nd = InternalBondPointers[BD_NX]!=NULL;
//			if (!pd && !nd) return 0;
//			else if (pd && !nd) return StrainPosDirsCur.x;
//			else if (!pd && nd) return StrainNegDirsCur.x;
//			else return 0.5*(StrainPosDirsCur.x + StrainNegDirsCur.x);
//			break;
//		case AXIS_Y:
//			pd = InternalBondPointers[BD_PY]!=NULL, nd = InternalBondPointers[BD_NY]!=NULL;
//			if (!pd && !nd) return 0;
//			else if (pd && !nd) return StrainPosDirsCur.y;
//			else if (!pd && nd) return StrainNegDirsCur.y;
//			else return 0.5*(StrainPosDirsCur.y + StrainNegDirsCur.y);
//			break;
//		case AXIS_Z:
//			pd = InternalBondPointers[BD_PZ]!=NULL, nd = InternalBondPointers[BD_NZ]!=NULL;
//			if (!pd && !nd) return 0;
//			else if (pd && !nd) return StrainPosDirsCur.z;
//			else if (!pd && nd) return StrainNegDirsCur.z;
//			else return 0.5*(StrainPosDirsCur.z + StrainNegDirsCur.z);
//			break;
//		default: return 0;
//
//	}
//
//
//}



Vec3D<> CVXS_Voxel::CalcTotalForce()
{
////	THE NEXT optimization target
//	//INTERNAL forces
//	Vec3D<> TotalForce = Vec3D<>(0,0,0);
//
////	Vec3D<> TotalForceB = Vec3D<>(0,0,0);
//	for (int i=0; i<6; i++){
//		if (InternalBondPointers[i]){
//			if (i%2==0){
//				TotalForce +=InternalBondPointers[i]->GetForce1();
//			}
//			else {
//				TotalForce +=InternalBondPointers[i]->GetForce2();
//			}
//		}
//	}
//	TotalForce = Angle.RotateVec3D(TotalForce);
//
//	TotalForce += -pSim->GetSlowDampZ() * Vel * _2xSqMxExS; //(2*sqrt(Mass*GetEMod()*Scale.x));
//
//	//Forces from collision bonds: To optimize!
//	if (pSim->IsFeatureEnabled(VXSFEAT_COLLISIONS)){
//		int NumColBond = ColBondPointers.size();
//		for (int i=0; i<NumColBond; i++){
//			if (IAmVox2Col(i)) TotalForce += ColBondPointers[i]->GetForce2();
//			else TotalForce += ColBondPointers[i]->GetForce1();
//		}
//	}
//
//	//Forced from input bond
//	TotalForce -= InputForce;
//
//
//
//	//From gravity
//	if (pSim->IsFeatureEnabled(VXSFEAT_GRAVITY))
//		TotalForce.z += Mass*pSim->pEnv->GetGravityAccel();
//
//	//EXTERNAL forces
//	TotalForce += ExternalInputScale*ExternalForce; //add in any external forces....
//
//
//	//if (pSim->IsFeatureEnabled(VXSFEAT_VOLUME_EFFECTS)){
//
//	//	updatePoissonsstrain();
//	//	Vec3D<> CurLocStrain = poissonsstrain;
////
//////		vfloat mu = GetPoisson();
////
//////		bool px = (InternalBondPointers[BD_PX] != NULL);
//////		bool nx = (InternalBondPointers[BD_NX] != NULL);
//////		bool py = (InternalBondPointers[BD_PY] != NULL);
//////		bool ny = (InternalBondPointers[BD_NY] != NULL);
//////		bool pz = (InternalBondPointers[BD_PZ] != NULL);
//////		bool nz = (InternalBondPointers[BD_NZ] != NULL);
////
//////		bool Tx = px && nx || ((px || nx) && (IS_FIXED(DOF_X, DofFixed) || ExternalForce.x != 0)); //if bond on both sides or pulling against a fixed or forced constraint
//////		bool Ty = py && ny || ((py || ny) && (IS_FIXED(DOF_Y, DofFixed) || ExternalForce.y != 0)); //if bond on both sides or pulling against a fixed or forced constraint
//////		bool Tz = pz && nz || ((pz || nz) && (IS_FIXED(DOF_Z, DofFixed) || ExternalForce.z != 0)); //if bond on both sides or pulling against a fixed or forced constraint
//////		CurLocStrain = strain(true);
////
//////		if (Tx){
//////			CurLocStrain.x = GetVoxelStrain(AXIS_X);
////////			if (px && !nx) CurLocStrain.x = StrainPosDirsCur.x;
////////			else if (!px && nx) CurLocStrain.x = StrainNegDirsCur.x;
////////			else CurLocStrain.x = 0.5*(StrainPosDirsCur.x + StrainNegDirsCur.x);
//////		}
//////		if (Ty){
//////			CurLocStrain.y = GetVoxelStrain(AXIS_Y);
//////
////////			if (py && !ny) CurLocStrain.y = StrainPosDirsCur.y;
////////			else if (!py && ny) CurLocStrain.y = StrainNegDirsCur.y;
////////			else CurLocStrain.y = 0.5*(StrainPosDirsCur.y + StrainNegDirsCur.y);
//////		}
//////		if (Tz){
//////			CurLocStrain.z = GetVoxelStrain(AXIS_Z);
//////
////////			if (pz && !nz) CurLocStrain.z = StrainPosDirsCur.z;
////////			else if (!pz && nz) CurLocStrain.z = StrainNegDirsCur.z;
////////			else CurLocStrain.z = 0.5*(StrainPosDirsCur.z + StrainNegDirsCur.z);
//////		}
////
//////		if (!Tx) CurLocStrain.x=0;
//////		if (!Ty) CurLocStrain.y=0;
//////		if (!Tz) CurLocStrain.z=0;
////
////
////		//account for poissons (simple?)
////		//Vec3D<> tmp(	CurLocStrain.x + pow(1+CurLocStrain.y + CurLocStrain.z, -mu)-1, 
////		//				CurLocStrain.y + pow(1+CurLocStrain.x + CurLocStrain.z, -mu)-1, 
////		//				CurLocStrain.z + pow(1+CurLocStrain.x + CurLocStrain.y, -mu)-1);
////
//////		Vec3D<> tmp2 = CurLocStrain;
//////		if (!Tx && !Ty && !Tz) tmp2 = Vec3D<>(0,0,0); //if nothing pushing or pulling, no strain on this bond!
//////		else if (!Tx && Ty && Tz) tmp2.x = pow(1+tmp2.y + tmp2.z, -mu)-1;
//////		else if (Tx && !Ty && Tz) tmp2.y = pow(1+tmp2.x + tmp2.z, -mu)-1; //??
//////		else if (Tx && Ty && !Tz) tmp2.z = pow(1+tmp2.x + tmp2.y, -mu)-1;
//////		else if (!Tx && !Ty && Tz) tmp2.x = tmp2.y = pow(1+tmp2.z, -mu)-1;
//////		else if (!Tx && Ty && !Tz) tmp2.x = tmp2.z = pow(1+tmp2.y, -mu)-1;
//////		else if (Tx && !Ty && !Tz) tmp2.y = tmp2.z = pow(1+tmp2.x, -mu)-1;
////		//else if (Tx && Ty && Tz) //we already have everything!
////
//////		poissonsstrain = tmp; //store to member
//////		CurLocStrain = tmp;
////
////		//TODO: get down to a force for this bond, then dump it (and current stiffness) to the bond
//////		for (int i=0; i<NumLocBond; i++){
////		for (int i=0; i<6; i++){
////////			CVXS_Bond* pThisBond = GetBond(i);
////			CVXS_Bond* pThisBond = InternalBondPointers[i];
////			if (!pThisBond) continue;
////
////			bool IAmVox1 = !IAmInternalVox2(i); //IsMe(pThisBond->GetpV1()); //otherwise vox 2 of the bond
////			switch (pThisBond->GetBondAxis()){
////			case AXIS_X:
////				if (IAmVox1) {
//////					pThisBond->TStrainSum1 = CurLocStrain.y + CurLocStrain.z;
////					pThisBond->CSArea1 = (1+CurLocStrain.y)*(1+CurLocStrain.z)*NominalSize*NominalSize;
////				}
////				else {
//////					pThisBond->TStrainSum2 = CurLocStrain.y + CurLocStrain.z;
////					pThisBond->CSArea2 = (1+CurLocStrain.y)*(1+CurLocStrain.z)*NominalSize*NominalSize;
////				}
////				break;
////			case AXIS_Y:
////				if (IAmVox1) {
//////					pThisBond->TStrainSum1 = CurLocStrain.x + CurLocStrain.z;
////					pThisBond->CSArea1 = (1+CurLocStrain.x)*(1+CurLocStrain.z)*NominalSize*NominalSize;
////				}
////				else {
//////					pThisBond->TStrainSum2 = CurLocStrain.x + CurLocStrain.z;
////					pThisBond->CSArea2 = (1+CurLocStrain.x)*(1+CurLocStrain.z)*NominalSize*NominalSize;
////				}
////				break;
////			case AXIS_Z:
////				if (IAmVox1) {
//////					pThisBond->TStrainSum1 = CurLocStrain.y + CurLocStrain.x;
////					pThisBond->CSArea1 = (1+CurLocStrain.y)*(1+CurLocStrain.x)*NominalSize*NominalSize;
////				}
////				else {
//////					pThisBond->TStrainSum2 = CurLocStrain.y + CurLocStrain.x;
////					pThisBond->CSArea2 = (1+CurLocStrain.y)*(1+CurLocStrain.x)*NominalSize*NominalSize;
////				}
////				break;
////			}
////		}
////	}
////	else { //volume effects off
//////		SizeCurrent = Vec3D<>(NominalSize, NominalSize, NominalSize);
////		//for (int i=0; i<NumLocBond; i++){
////		//	CVXS_Bond* pThisBond = GetBond(i);
////		for (int i=0; i<6; i++){ //update for collision bonds?
////			CVXS_Bond* pThisBond = InternalBondPointers[i];
////			if (pThisBond){
////				pThisBond->CSArea1 = pThisBond->CSArea2 = NominalSize*NominalSize;
////			}
////		}
////	}
////
////	if(pSim->pEnv->IsFloorEnabled()){
//	if(pSim->IsFeatureEnabled(VXSFEAT_FLOOR)){
//		TotalForce += CalcFloorEffect(Vec3D<vfloat>(TotalForce));
//		//else StaticFricFlag = false;
//	//	if (StaticFricFlag) {TotalForce.x = 0; TotalForce.y = 0;} //no lateral movement if static friction in effect
//	}
//
//	
//	CornerPosCur = (Vec3D<>(1,1,1)+poissonsstrain)*NominalSize/2;
//	CornerNegCur = -(Vec3D<>(1,1,1)+poissonsstrain)*NominalSize/2;
//
////	CornerPosCur = (Vec3D<>(1,1,1)+StrainPosDirsCur)*NominalSize/2;
////	CornerNegCur = -(Vec3D<>(1,1,1)+StrainNegDirsCur)*NominalSize/2;
//
//	//Enforce fixed degrees of freedom (put no force on them so they don't move)
////	if (IS_FIXED(DOF_X, DofFixed) && WithRestraint) TotalForce.x=0;
////	if (IS_FIXED(DOF_Y, DofFixed) && WithRestraint) TotalForce.y=0;
////	if (IS_FIXED(DOF_Z, DofFixed) && WithRestraint) TotalForce.z=0;
//
//	ForceCurrent=TotalForce;
//	return ForceCurrent;
return Vec3D<>(0,0,0);
}

Vec3D<> CVXS_Voxel::CalcTotalMoment(void)
{
//	Vec3D<> TotalMoment(0,0,0);
//	//permanent bonds
//	for (int i=0; i<6; i++){ //update for collision bonds?
//		CVXS_Bond* pThisBond = InternalBondPointers[i];
//		if (pThisBond){
//			if (IAmInternalVox2(i)){ TotalMoment -= InternalBondPointers[i]->GetMoment2(); } //if this is voxel 2		//add moments from bond
//			else { TotalMoment -= InternalBondPointers[i]->GetMoment1(); } //if this is voxel 1
//		}
//	}
//
//	//Moments from permanent bonds:
//	//for (int i=0; i<3; i++){
//	//	CVXS_Bond* pThisBondP = InternalBondPointers[2*i];
//	//	CVXS_Bond* pThisBondN = InternalBondPointers[2*i+1];
//
//	//	if (!pThisBondP && !pThisBondN) continue;
//	//	Vec3D<> tmpMoment;
//	//	
//	//	if (pThisBondP) tmpMoment += pThisBondP->GetMoment1();
//	//	if (pThisBondN) tmpMoment += pThisBondN->GetMoment2();
//
//	//	TotalMoment -= tmpMoment;
//	//	if (i==0) TotalMoment -= tmpMoment;
//	//}
//	TotalMoment = Angle.RotateVec3D(TotalMoment);
//
//
////	TotalMoment += -pSim->GetSlowDampZ() * Vel*_2xSqMxExS; 
////
////				vfloat AngMomFact = (1 - 10*pSim->GetSlowDampZ() * _inertiaInv *_2xSqIxExSxSxS*dt);
//
//
//
//	//EXTERNAL moments
//	TotalMoment += -pSim->GetSlowDampZ() * AngVel *_2xSqIxExSxSxS; 
//	TotalMoment += ExternalInputScale*ExternalTorque; //add in any external forces....
//
//	if (IS_FIXED(DOF_TX, DofFixed)) TotalMoment.x=0;
//	if (IS_FIXED(DOF_TY, DofFixed)) TotalMoment.y=0;
//	if (IS_FIXED(DOF_TZ, DofFixed)) TotalMoment.z=0;
//
//	return TotalMoment;
return Vec3D<>(0,0,0);

}

vfloat CVXS_Voxel::GetCurGroundPenetration() //how far into the ground penetrating (penetration is positive, no penetration is zero)
{
	vfloat Penetration = 0.5*Scale - Pos.z;
	return Penetration <= 0 ? 0 : Penetration;
}

inline bool CVXS_Voxel::IAmVox2Col(const int BondDirColIndex) const 
{
	return false;
//	return (this == ColBondPointers[BondDirColIndex]->GetpV2());
} //returns true if this voxel is Vox2 of the specified bond


Vec3D<> CVXS_Voxel::CalcFloorEffect(Vec3D<> TotalVoxForce) //calculates the object's interaction with a floor. should be calculated AFTER all other forces for static friction to work right...
{
	Vec3D<> FloorForce(0,0,0); //the force added by floor interactions...

	//StaticFricFlag = false; //assume not under static friction unless we decide otherwise
	vfloat CurPenetration = GetCurGroundPenetration();


	if (CurPenetration>0){ 
		vfloat LocA1 = GetLinearStiffness(); //p_Sim->LocalVXC.GetBaseMat(MatIndex)->GetElasticMod()*2*NominalSize; 
		vfloat LocUDynamic = _pMat->GetuDynamic();
		vfloat LocUStatic = _pMat->GetuStatic();

		vfloat NormalForce = LocA1 * CurPenetration; //positive for penetration...
		FloorForce.z += NormalForce; //force resisting penetration
	
		//do vertical damping here...
		FloorForce.z -= pSim->GetCollisionDampZ()*_2xSqMxExS*Vel.z;  //critically damp force for this bond to ground
//		FloorForce.z -= p_Sim->GetCollisionDampZ()*2*Mass*sqrt(LocA1/Mass)*Vel.z;  //critically damp force for this bond to ground
		
		//lateral friction
		vfloat SurfaceVel = sqrt(Vel.x*Vel.x + Vel.y*Vel.y); //velocity along the floor...
		vfloat SurfaceVelAngle = atan2(Vel.y, Vel.x); //angle of sliding along floor...
		vfloat SurfaceForce = sqrt(TotalVoxForce.x*TotalVoxForce.x + TotalVoxForce.y*TotalVoxForce.y);
	
		//alwyas acts in direction opposite to force in STATIC friction mode

		if (StaticFricFlag){
			if (SurfaceForce > LocUStatic*NormalForce) StaticFricFlag = false; //if we have enough to break static friction
		}

//		if (Vel.x == 0 && Vel.y == 0){ //STATIC FRICTION: if this point is stopped and in the static friction mode...
//			if (SurfaceForce < LocUStatic*NormalForce) StaticFricFlag = true; //if we don't have enough to break static friction
//		}
		else { //DYNAMIC FRICTION
//			if (dFrictionForce*pSim->dt < Mass*SurfaceVel){ //check momentum. if we are not possibly coming to a stop this timestep, add in the friction force
			vfloat dFrictionForce = LocUDynamic*NormalForce; 
			Vec3D<> FricForceToAdd = -Vec3D<>(cos(SurfaceVelAngle)*dFrictionForce, sin(SurfaceVelAngle)*dFrictionForce, 0); //always acts in direction opposed to velocity in DYNAMIC friction mode	
			FloorForce += FricForceToAdd;
//			}
//			else { //if we are coming to a stop, don't overshoot the stop. Set to zero and zero the momentum to get static friction to kick in.
//				StaticFricFlag = true;
//				LinMom.x = 0; //fully stop the voxel here! (caution...)
//				LinMom.y = 0;
//			}
		}
		
	}
	else {StaticFricFlag = false;}
	return FloorForce;

}

Vec3D<> CVXS_Voxel::CalcGndDampEffect() //damps everything to ground as qucik as possible...
{
	vfloat tmp = 0;
//	for (int i=0; i<6; i++) tmp+=sqrt(InternalBondPointers[i]->GetLinearStiffness()*Mass);

	return -pSim->GetSlowDampZ()*2*tmp*Vel;
}


vfloat CVXS_Voxel::GetMaxBondStrain(void) const
{
	vfloat MxSt = 0;
	for (int i=0; i<6; i++){
		if (InternalBondPointers[i] != NULL){
			vfloat TSt = 0; //InternalBondPointers[i]->GetEngStrain();
			if (TSt>MxSt) MxSt = TSt; 
		}
	}
	return MxSt;
}

vfloat CVXS_Voxel::GetMaxBondStrainE(void) const
{
	vfloat MxSt = 0;
	for (int i=0; i<6; i++){
		if (InternalBondPointers[i] != NULL){
			vfloat TSt = 0; //InternalBondPointers[i]->GetStrainEnergy();
			if (TSt>MxSt) MxSt = TSt; 
		}
	}
	return MxSt;
}

vfloat CVXS_Voxel::GetMaxBondStress(void) const
{
	vfloat MxSt = 0;
//	for (int i=0; i<NumLocalBonds; i++){
	for (int i=0; i<6; i++){
		if (InternalBondPointers[i] != NULL){
			vfloat TSt = 0; //InternalBondPointers[i]->GetEngStress();
			if (TSt>MxSt) MxSt = TSt; 
		}
	}
	return MxSt;

}

vfloat CVXS_Voxel::CalcVoxMatStress(const vfloat StrainIn, bool* const IsPastYielded, bool* const IsPastFail) const
{
		return _pMat->GetModelStress(StrainIn, IsPastYielded, IsPastFail);
}

Vec3D<> CVXS_Voxel::strain(bool tensionStrain)
{
	////if no connections in the positive and negative directions of a particular axis, strain is zero
	////if one connection in positive or negative direction of a particular axis, strain is that strain - ?? and force or constraint?
	////if connections in both the positive and negative directions of a particular axis, strain is the average. 
	//float intStrRet[3] = {0}; //intermediate strain return value. axes according to linkAxis enum
	//int numBondAxis[3] = {0}; //number of bonds in this axis (0,1,2). axes according to linkAxis enum
	//for (int i=0; i<6; i++){ //cycle through link directions
	//	if (InternalBondPointers[i]){
	//		int axis = toAxis((linkDirection)i);
	//		intStrRet[axis] += isNegative((linkDirection)i) ? InternalBondPointers[i]->GetStrainV2() : InternalBondPointers[i]->GetStrainV1();
	//		numBondAxis[axis]++;
	//	}
	//}
	//for (int i=0; i<3; i++){ //cycle through axes
	//	if (numBondAxis[i]==2) intStrRet[i]*= 0.5f; //average
	//	if (tensionStrain && numBondAxis[i]==1){ //if just one bond
	//		if (!IS_FIXED(DOF_X, DofFixed) && ExternalForce.x == 0) intStrRet[i]=0; //if no other external means of providing tension, zero out strain.
	//	}
	//}


	//return Vec3D<>(intStrRet[0], intStrRet[1], intStrRet[2]);
	return Vec3D<>(0,0,0);
}

void CVXS_Voxel::updatePoissonsstrain()
{
	if (pSim->IsFeatureEnabled(VXSFEAT_VOLUME_EFFECTS)){
		vfloat mu = GetPoisson();
		poissonsstrain = strain(true);
		Vec3D<> tmp(	poissonsstrain.x + pow(1+poissonsstrain.y + poissonsstrain.z, -mu)-1, 
						poissonsstrain.y + pow(1+poissonsstrain.x + poissonsstrain.z, -mu)-1, 
						poissonsstrain.z + pow(1+poissonsstrain.x + poissonsstrain.y, -mu)-1);
		poissonsstrain = tmp;
	}
	else poissonsstrain = strain(false);
}
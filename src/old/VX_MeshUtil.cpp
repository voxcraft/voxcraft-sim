/*******************************************************************************
Copyright (c) 2010, Jonathan Hiller (Cornell University)
If used in publication cite "J. Hiller and H. Lipson "Dynamic Simulation of Soft Heterogeneous Objects" In press. (2011)"

This file is part of Voxelyze.
Voxelyze is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Voxelyze is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
See <http://www.opensource.org/licenses/lgpl-3.0.html> for license details.
*******************************************************************************/

#include "VX_MeshUtil.h"
#include "VX_Object.h"
#include "VX_FEA.h"
#include "VX_Sim.h"
#include "VXS_SimGLView.h"
#include "MarchCube.h"

#include "VX_MeshRender.h"


CVX_MeshUtil::CVX_MeshUtil(void)
{
	pSim = NULL;
	pSimView = NULL;
	pFEA = NULL;
}

CVX_MeshUtil::~CVX_MeshUtil(void)
{
	Clear();
}

void CVX_MeshUtil::Clear(void)
{
	DefMesh.Clear();
	CalcVerts.clear();
	FacetToSIndex.clear();
	pSim = NULL;
	pFEA = NULL;
}


void CVX_MeshUtil::LinkSimExisting(CVX_Sim* pSimIn, CVXS_SimGLView* pSimViewIn, CMesh* pMeshIn) //links the simulation to an existing surface mesh so it can be deformed correctly
{
	if (pSimIn != NULL) pSim = pSimIn;
	if (pSimViewIn != NULL) pSimView = pSimViewIn;
	if (pMeshIn != NULL) DefMesh = *pMeshIn;
	if (!DefMesh.Exists()) return; //cna't inport if no mesh...

	vfloat FuzzDist = 1.5; //multiple of Lattice dimension to include voxels (radius)

	CVX_Object* pObj = pSimIn->pEnv->pObj; //temp, for convenience
	Vec3D<> WS = pObj->GetWorkSpace();
	Vec3D<> tmpVert, VoxCtr, Offset;
	int xi, yi, zi, Vi;
	vfloat dist;

	CalcVerts.clear();
	for (int i=0; i<(int)(DefMesh.Vertices.size()); i++){
		CVertexCalc tmpVertCalc;
		tmpVert = DefMesh.Vertices[i].v;

		//check for all voxel centerpoints within one voxel of the voxel this one would be contained in?
		//see which voxe this point is containted in...
		xi = (int)(tmpVert.x/WS.x*pObj->GetVXDim());
		yi = (int)(tmpVert.y/WS.y*pObj->GetVYDim());
		zi = (int)(tmpVert.z/WS.z*pObj->GetVZDim());
		int D2S = int(FuzzDist*(pObj->GetLatticeDim()/pObj->GetLatDimEnv().Min())); //distance to search... make sure we look far enough if some dimensions are smaller than lattice dimension
		for (int x=xi-D2S; x<=xi+D2S; x++){
			for (int y=yi-D2S; y<=yi+D2S; y++){
				for (int z=zi-D2S; z<=zi+D2S; z++){
					Vi = pObj->GetIndex(x, y, z);
					if (Vi != -1 && pObj->GetMat(Vi) != 0){ //if valid voxel location
						VoxCtr = pObj->GetXYZ(Vi);
						Offset = tmpVert-VoxCtr;
						dist = Offset.Length()/pObj->GetLatticeDim();
						if (dist<FuzzDist){ //FINALLY! make a connection!
							tmpVertCalc.ConVoxels.push_back(CVertexComp(Vi, Offset, (FuzzDist-dist)/FuzzDist));
						}
					}
				}
			}
		}
		CalcVerts.push_back(tmpVertCalc);
	}

	DefMesh.DrawSmooth = true;

	pSimView->defMeshVx2->generateMesh(); /////
//	pSimView->defMeshVx2->saveObj("test.obj");

}

void CVX_MeshUtil::LinkSimSmooth(CVX_Sim* pSimIn, CVXS_SimGLView* pSimViewIn) //creates a voxel mesh and links to simulation
{
	CMesh GeneratedSmoothMesh;
	CVX_Object* pObj = pSimIn->pEnv->pObj;
	int xD = pObj->GetVXDim();
	int yD = pObj->GetVYDim();
	int zD = pObj->GetVZDim();

	CArray3D<float> OccupancyArray;
	OccupancyArray.resize(Index3D(xD, yD, zD)); 
	int NumPossibleVox = pObj->GetStArraySize();
	for (int k=0; k<zD; k++){
		for (int j=0; j<yD; j++){
			for (int i=0; i<xD; i++){
				if (pObj->Structure.GetData(pObj->GetIndex(i,j,k))>0){
					if (j==19)
						int stophere=-3;
					OccupancyArray.addValue(Index3D(i,j,k), 1.0);
				}
			}
		}
	}
	CMarchCube::SingleMaterial(&GeneratedSmoothMesh, &OccupancyArray, 0.5, pSimIn->pEnv->pObj->GetLatticeDim());
	LinkSimExisting(pSimIn, pSimViewIn, &GeneratedSmoothMesh);
}

void CVX_MeshUtil::LinkSimVoxels(CVX_Sim* pSimIn, CVXS_SimGLView* pSimViewIn)
{
	if (pSimIn != NULL) pSim = pSimIn;
	if (pSimViewIn != NULL) pSimView = pSimViewIn;

	CMesh tmpMesh;

	CVX_Object* pDM = pSim->pEnv->pObj; //for convenience of access...

	//make all the potential lattice points:
//	std::vector<CVertexCalc> tmpVerts; //3D array [x][y][z]
	CalcVertsAll.clear();
	int tVX = pDM->GetVXDim()+1;
	int tVY = pDM->GetVYDim()+1;
	int tVZ = pDM->GetVZDim()+1;
	Vec3D<> Offset, tmp;

	FacetToSIndex.clear();
	CalcVertsAll.resize(tVX*tVY*tVZ); //each row in x, rolling around to consecutive Y's, then z.

	int NumVox = 0; /////////////// pSim->NumVox();
	for (int i=0; i<NumVox; i++){ //for every voxel in the simulation...
		//CVXS_Voxel& CurVox = pSim->VoxArray[i];
		//int CurXIndex = pSim->StoXIndexMap[i];
		//int cX, cY, cZ, tX, tY, tZ;
		//pDM->GetXYZNom(&cX, &cY, &cZ, CurXIndex); //get the nominal coords if this voxel...

		//for (int i=0; i<8; i++){
		//	int tmpVertIndex = D3IndCorner(cX, cY, cZ, tVX-1, tVY-1, i);
		//	CalcVertsAll[tmpVertIndex].ConVoxels.push_back(CVertexComp(CurXIndex, i));
		//}

		//	//look at all faces, add facets (no non-existant vertex indices.. but we'll shift down
		//	if (!pDM->IsVoxRelativeIndex(CurXIndex, 1, 0, 0)){ //if we have an open face +X
		//		tmpMesh.Facets.push_back(CFacet(Vec3D<>(1, 0, 0), D3Ind(cX+1, cY, cZ, tVX, tVY), D3Ind(cX+1, cY+1, cZ, tVX, tVY), D3Ind(cX+1, cY+1, cZ+1, tVX, tVY), CurXIndex));
		//		tmpMesh.Facets.push_back(CFacet(Vec3D<>(1, 0, 0), D3Ind(cX+1, cY, cZ, tVX, tVY), D3Ind(cX+1, cY+1, cZ+1, tVX, tVY), D3Ind(cX+1, cY, cZ+1, tVX, tVY), CurXIndex));
		//	}
		//	if (!pDM->IsVoxRelativeIndex(CurXIndex, -1, 0, 0)){ //if we have an open face -X
		//		tmpMesh.Facets.push_back(CFacet(Vec3D<>(-1, 0, 0), D3Ind(cX, cY, cZ, tVX, tVY), D3Ind(cX, cY+1, cZ+1, tVX, tVY), D3Ind(cX, cY+1, cZ, tVX, tVY), CurXIndex));
		//		tmpMesh.Facets.push_back(CFacet(Vec3D<>(-1, 0, 0), D3Ind(cX, cY, cZ, tVX, tVY), D3Ind(cX, cY, cZ+1, tVX, tVY), D3Ind(cX, cY+1, cZ+1, tVX, tVY), CurXIndex));

		//	}
		//	if (!pDM->IsVoxRelativeIndex(CurXIndex, 0, 1, 0)){ //if we have an open face + Y
		//		tmpMesh.Facets.push_back(CFacet(Vec3D<>(0, 1, 0), D3Ind(cX, cY+1, cZ, tVX, tVY), D3Ind(cX, cY+1, cZ+1, tVX, tVY), D3Ind(cX+1, cY+1, cZ+1, tVX, tVY), CurXIndex));
		//		tmpMesh.Facets.push_back(CFacet(Vec3D<>(0, 1, 0), D3Ind(cX, cY+1, cZ, tVX, tVY), D3Ind(cX+1, cY+1, cZ+1, tVX, tVY), D3Ind(cX+1, cY+1, cZ, tVX, tVY), CurXIndex));
		//	}
		//	if (!pDM->IsVoxRelativeIndex(CurXIndex, 0, -1, 0)){ //if we have an open face + Y
		//		tmpMesh.Facets.push_back(CFacet(Vec3D<>(0, -1, 0), D3Ind(cX, cY, cZ, tVX, tVY), D3Ind(cX+1, cY, cZ+1, tVX, tVY), D3Ind(cX, cY, cZ+1, tVX, tVY), CurXIndex));
		//		tmpMesh.Facets.push_back(CFacet(Vec3D<>(0, -1, 0), D3Ind(cX, cY, cZ, tVX, tVY), D3Ind(cX+1, cY, cZ, tVX, tVY), D3Ind(cX+1, cY, cZ+1, tVX, tVY), CurXIndex));
		//	}
		//	if (!pDM->IsVoxRelativeIndex(CurXIndex, 0, 0, 1)){ //if we have an open face + Y
		//		tmpMesh.Facets.push_back(CFacet(Vec3D<>(0, 0, 1), D3Ind(cX, cY, cZ+1, tVX, tVY), D3Ind(cX+1, cY, cZ+1, tVX, tVY), D3Ind(cX+1, cY+1, cZ+1, tVX, tVY), CurXIndex));
		//		tmpMesh.Facets.push_back(CFacet(Vec3D<>(0, 0, 1), D3Ind(cX, cY, cZ+1, tVX, tVY), D3Ind(cX+1, cY+1, cZ+1, tVX, tVY), D3Ind(cX, cY+1, cZ+1, tVX, tVY), CurXIndex));
		//	}
		//	if (!pDM->IsVoxRelativeIndex(CurXIndex, 0, 0, -1)){ //if we have an open face + Y
		//		tmpMesh.Facets.push_back(CFacet(Vec3D<>(0, 0, -1), D3Ind(cX, cY, cZ, tVX, tVY), D3Ind(cX+1, cY+1, cZ, tVX, tVY), D3Ind(cX+1, cY, cZ, tVX, tVY), CurXIndex));
		//		tmpMesh.Facets.push_back(CFacet(Vec3D<>(0, 0, -1), D3Ind(cX, cY, cZ, tVX, tVY), D3Ind(cX, cY+1, cZ, tVX, tVY), D3Ind(cX+1, cY+1, cZ, tVX, tVY), CurXIndex));
		//	}


		//	while (FacetToSIndex.size() < tmpMesh.Facets.size()) //update our linking list with all the facets we just added
		//		FacetToSIndex.push_back(pSim->XtoSIndexMap[CurXIndex]);


		//	//look at all 12 edges... (shared lines added twice, but oh well... we will delete them later
		//	//vertical lines
		//	if (!pDM->IsVoxRelativeIndex(CurXIndex, 1, 0, 0) || !pDM->IsVoxRelativeIndex(CurXIndex, 0, 1, 0))
		//		tmpMesh.Lines.push_back(CLine(D3Ind(cX+1, cY+1, cZ, tVX, tVY), D3Ind(cX+1, cY+1, cZ+1, tVX, tVY)));
		//	if (!pDM->IsVoxRelativeIndex(CurXIndex, -1, 0, 0) || !pDM->IsVoxRelativeIndex(CurXIndex, 0, 1, 0))
		//		tmpMesh.Lines.push_back(CLine(D3Ind(cX, cY+1, cZ, tVX, tVY), D3Ind(cX, cY+1, cZ+1, tVX, tVY)));
		//	if (!pDM->IsVoxRelativeIndex(CurXIndex, -1, 0, 0) || !pDM->IsVoxRelativeIndex(CurXIndex, 0, -1, 0))
		//		tmpMesh.Lines.push_back(CLine(D3Ind(cX, cY, cZ, tVX, tVY), D3Ind(cX, cY, cZ+1, tVX, tVY)));
		//	if (!pDM->IsVoxRelativeIndex(CurXIndex, 1, 0, 0) || !pDM->IsVoxRelativeIndex(CurXIndex, 0, -1, 0))
		//		tmpMesh.Lines.push_back(CLine(D3Ind(cX+1, cY, cZ, tVX, tVY), D3Ind(cX+1, cY, cZ+1, tVX, tVY)));

		//	//top lines
		//	if (!pDM->IsVoxRelativeIndex(CurXIndex, 1, 0, 0) || !pDM->IsVoxRelativeIndex(CurXIndex, 0, 0, 1))
		//		tmpMesh.Lines.push_back(CLine(D3Ind(cX+1, cY, cZ+1, tVX, tVY), D3Ind(cX+1, cY+1, cZ+1, tVX, tVY)));
		//	if (!pDM->IsVoxRelativeIndex(CurXIndex, 0, 1, 0) || !pDM->IsVoxRelativeIndex(CurXIndex, 0, 0, 1))
		//		tmpMesh.Lines.push_back(CLine(D3Ind(cX, cY+1, cZ+1, tVX, tVY), D3Ind(cX+1, cY+1, cZ+1, tVX, tVY)));
		//	if (!pDM->IsVoxRelativeIndex(CurXIndex, -1, 0, 0) || !pDM->IsVoxRelativeIndex(CurXIndex, 0, 0, 1))
		//		tmpMesh.Lines.push_back(CLine(D3Ind(cX, cY, cZ+1, tVX, tVY), D3Ind(cX, cY+1, cZ+1, tVX, tVY)));
		//	if (!pDM->IsVoxRelativeIndex(CurXIndex, 0, -1, 0) || !pDM->IsVoxRelativeIndex(CurXIndex, 0, 0, 1))
		//		tmpMesh.Lines.push_back(CLine(D3Ind(cX, cY, cZ+1, tVX, tVY), D3Ind(cX+1, cY, cZ+1, tVX, tVY)));

		//	//bottom lines
		//	if (!pDM->IsVoxRelativeIndex(CurXIndex, 1, 0, 0) || !pDM->IsVoxRelativeIndex(CurXIndex, 0, 0, -1))
		//		tmpMesh.Lines.push_back(CLine(D3Ind(cX+1, cY, cZ, tVX, tVY), D3Ind(cX+1, cY+1, cZ, tVX, tVY)));
		//	if (!pDM->IsVoxRelativeIndex(CurXIndex, 0, 1, 0) || !pDM->IsVoxRelativeIndex(CurXIndex, 0, 0, -1))
		//		tmpMesh.Lines.push_back(CLine(D3Ind(cX, cY+1, cZ, tVX, tVY), D3Ind(cX+1, cY+1, cZ, tVX, tVY)));
		//	if (!pDM->IsVoxRelativeIndex(CurXIndex, -1, 0, 0) || !pDM->IsVoxRelativeIndex(CurXIndex, 0, 0, -1))
		//		tmpMesh.Lines.push_back(CLine(D3Ind(cX, cY, cZ, tVX, tVY), D3Ind(cX, cY+1, cZ, tVX, tVY)));
		//	if (!pDM->IsVoxRelativeIndex(CurXIndex, 0, -1, 0) || !pDM->IsVoxRelativeIndex(CurXIndex, 0, 0, -1))
		//		tmpMesh.Lines.push_back(CLine(D3Ind(cX, cY, cZ, tVX, tVY), D3Ind(cX+1, cY, cZ, tVX, tVY)));
	}


	Vec3D<> tmpOffset = pDM->GetLatDimEnv()/2;
	CVX_Object tmpDM = *pDM;
	tmpDM.Resize(tVX, tVY, tVZ);
	//add all vertices..
	for (int k=0; k<tVZ; k++){
		for (int j=0; j<tVY; j++){
			for (int i=0; i<tVX; i++){
				Vec3D<> pos = tmpDM.GetXYZ(i, j, k)-tmpOffset;
				tmpMesh.Vertices.push_back(CVertex(pos));
			}
		}
	}

	tmpMesh.RemoveDupLines();

	//tmpMesh is now complete, but with way too many vertices, most unused.

	//commit to DefMesh, minus all the extra vertices!!!
	DefMesh.Clear();
	CalcVerts.clear();

	int* Map = new int[tVX*tVY*tVZ];
	int NewInd = 0;
	for (int i=0; i<tVX*tVY*tVZ; i++){
		if (CalcVertsAll[i].ConVoxels.size() != 0 && CalcVertsAll[i].ConVoxels.size() != 8){
			DefMesh.Vertices.push_back(tmpMesh.Vertices[i]);
			CalcVerts.push_back(CalcVertsAll[i]);
			Map[i] = NewInd++;

		}
		else Map[i] = -1;
	}

	DefMesh.Facets = tmpMesh.Facets;
	for (int i=0; i<(int)DefMesh.Facets.size(); i++){
		for (int j=0; j<3; j++){
			DefMesh.Facets[i].vi[j] = Map[tmpMesh.Facets[i].vi[j]];
		}
	}

	DefMesh.Lines = tmpMesh.Lines;
	for (int i=0; i<(int)DefMesh.Lines.size(); i++){
		for (int j=0; j<2; j++){
			DefMesh.Lines[i].vi[j] = Map[tmpMesh.Lines[i].vi[j]];
		}
	}
}

void CVX_MeshUtil::UpdateMesh(int CurSel) //updates mesh based on linked FEA/Relaxation
{

	if (pSim){/*
		//update vertex positions
		Vec3D avgPos, ThisPos;
		vfloat TotWeight, Ra, Ga, Ba;
		for (int i=0; i<(int)CalcVerts.size(); i++){
			avgPos = Vec3D(0,0,0);
			TotWeight = 0.0;
#ifdef USE_OPEN_GL
			CColor Tmp;
			Ra = 0.0;
			Ga = 0.0;
			Ba = 0.0;
#endif

			for (int j=0; j<(int)CalcVerts[i].ConVoxels.size(); j++){
				CVXS_Voxel& CurVox = pSim->VoxArray[pSim->XtoSIndexMap[CalcVerts[i].ConVoxels[j].XIndex]];

				//temp, until scale is updated by temp!
				vfloat TempFact = 1.0;
				if(pSim->pEnv->TempEnabled) TempFact = (1+(pSim->pEnv->CurTemp-pSim->pEnv->TempBase)*CurVox.GetCTE()); //pSim->LocalVXC.GetBaseMat(CurVox.MatIndex)->GetCTE());

				ThisPos = CurVox.GetCurPos() + (CurVox.Angle*Quat3D(TempFact*CalcVerts[i].ConVoxels[j].Off)*CurVox.Angle.Conjugate()).ToVec(); //todo:scale!
				avgPos += CalcVerts[i].ConVoxels[j].Weight*ThisPos;

#ifdef USE_OPEN_GL
				Tmp = pSim->GetCurColor(pSim->XtoSIndexMap[CalcVerts[i].ConVoxels[j].XIndex], CurSel);
				Ra += CalcVerts[i].ConVoxels[j].Weight*Tmp.r;
				Ga += CalcVerts[i].ConVoxels[j].Weight*Tmp.g;
				Ba += CalcVerts[i].ConVoxels[j].Weight*Tmp.b;
#endif

				TotWeight += CalcVerts[i].ConVoxels[j].Weight;
			}

			DefMesh.Vertices[i].DrawOffset = avgPos/TotWeight - DefMesh.Vertices[i].v; //yes, we're subtracting this out just to add it back.

#ifdef USE_OPEN_GL
			DefMesh.Vertices[i].VColor = CColor(Ra/TotWeight, Ga/TotWeight, Ba/TotWeight, 1.0);
#endif
		}*/

		Vec3D<> NewPos;
		for (int i=0; i<(int)CalcVerts.size(); i++){
			GetCurVLoc(CalcVerts[i], &NewPos);
			DefMesh.Vertices[i].DrawOffset = NewPos - DefMesh.Vertices[i].v; //yes, we're subtracting this out just to add it back.

#ifdef USE_OPEN_GL
			GetCurVCol(CalcVerts[i], &DefMesh.Vertices[i].VColor, CurSel);
#endif
		}


		//update colors that aren't by vertex!
#ifdef USE_OPEN_GL

		if (FacetToSIndex.size() != 0){ //if pulling color for full fact
			for (int i=0; i<(int)DefMesh.Facets.size(); i++){
				DefMesh.Facets[i].FColor = pSimView->GetCurVoxColor(FacetToSIndex[i], CurSel);
				//DefMesh.Facets[i].FColor = pSim->GetCurColor(i, CurSel);
			}
		}
		//else { //color by vertex!!!...

		//}
#endif

		//update normals!
		if (!DefMesh.DrawSmooth) //if drawing faces...
			DefMesh.CalcFaceNormals();
		else
			DefMesh.CalcVertNormals();
		

	}
	else if (pFEA){
		//todo...
	}
}

void CVX_MeshUtil::GetCurVLoc(CVertexCalc& VertCalc, Vec3D<>* pLocOut)
{
//	Vec3D<> avgPos(0,0,0);
//	Vec3D<> ThisPos;
//	vfloat TotWeight = 0;
//
//	for (int j=0; j<(int)VertCalc.ConVoxels.size(); j++){
//		CVertexComp CurComp = VertCalc.ConVoxels[j];
//		CVXS_Voxel& CurVox = pSim->VoxArray[pSim->XtoSIndexMap[CurComp.XIndex]];
//
//		Vec3D<> ThisOffset; //offset to this point
//		if (CurComp.Corner != -1){
//			switch (CurComp.Corner){
//			case cNNN: ThisOffset = CurVox.GetCornerNeg(); break;
//			case cNNP: ThisOffset = Vec3D<>(CurVox.GetCornerNeg().x, CurVox.GetCornerNeg().y, CurVox.GetCornerPos().z); break;
//			case cNPN: ThisOffset = Vec3D<>(CurVox.GetCornerNeg().x, CurVox.GetCornerPos().y, CurVox.GetCornerNeg().z); break;
//			case cNPP: ThisOffset = Vec3D<>(CurVox.GetCornerNeg().x, CurVox.GetCornerPos().y, CurVox.GetCornerPos().z); break;
//			case cPNN: ThisOffset = Vec3D<>(CurVox.GetCornerPos().x, CurVox.GetCornerNeg().y, CurVox.GetCornerNeg().z); break;
//			case cPNP: ThisOffset = Vec3D<>(CurVox.GetCornerPos().x, CurVox.GetCornerNeg().y, CurVox.GetCornerPos().z); break;
//			case cPPN: ThisOffset = Vec3D<>(CurVox.GetCornerPos().x, CurVox.GetCornerPos().y, CurVox.GetCornerNeg().z); break;
//			case cPPP: ThisOffset = CurVox.GetCornerPos(); break;
//			}
//
//
//		}
//		else {
//			vfloat ScaleFact = CurVox.GetCurScale() / CurVox.GetNominalSize(); //Assumes square, isotropic expansion
//			ThisOffset = ScaleFact*VertCalc.ConVoxels[j].Off;
//
//		//ThisPos = CurVox.GetCurPos() + (CurVox.GetCurAngle()*Quat3D<>(ScaleFact*VertCalc.ConVoxels[j].Off)*CurVox.GetCurAngle().Conjugate()).ToVec(); //todo:scale!
//		}
//		
////		ThisPos = CurVox.GetCurPos() + (CurVox.GetCurAngle()*Quat3D<>(ThisOffset)*CurVox.GetCurAngle().Conjugate()).ToVec(); //todo:scale!
//		ThisPos = CurVox.GetCurPos() + CurVox.GetCurAngle().RotateVec3D(ThisOffset);
//
//		avgPos += VertCalc.ConVoxels[j].Weight*ThisPos;
//		TotWeight += VertCalc.ConVoxels[j].Weight;
//	}
//
//	*pLocOut = avgPos/TotWeight; 
}

void CVX_MeshUtil::GetCurVCol(CVertexCalc& VertCalc, CColor* pColOut, int CurSel)
{
//#ifdef USE_OPEN_GL
//	vfloat TotWeight = 0.0;
//	vfloat Ra = 0.0f;
//	vfloat Ga = 0.0f;
//	vfloat Ba = 0.0f;
//	CColor Tmp;
//
//	for (int j=0; j<(int)VertCalc.ConVoxels.size(); j++){
//		Tmp = pSimView->GetCurVoxColor(pSim->XtoSIndexMap[VertCalc.ConVoxels[j].XIndex], CurSel);
//		Ra += VertCalc.ConVoxels[j].Weight*Tmp.r;
//		Ga += VertCalc.ConVoxels[j].Weight*Tmp.g;
//		Ba += VertCalc.ConVoxels[j].Weight*Tmp.b;
//
//		TotWeight += VertCalc.ConVoxels[j].Weight;
//	}
//
//
//	*pColOut = CColor(Ra/TotWeight, Ga/TotWeight, Ba/TotWeight, 1.0);
//#endif
}

int CVX_MeshUtil::D3IndCorner(int XInd, int YInd, int ZInd, const int XSize, const int YSize, int const Corner) //returns the vertex index in the [x+1, y+1, z+1] CalcVertsAll array corresponding to the specified corner of this voxel.
{
	switch (Corner){
//	case cNNN: *pOut = Vec3D(-1,-1,-1); break;
	case cNNP: ZInd++; break;
	case cNPN: YInd++; break;
	case cNPP: YInd++; ZInd++; break;
	case cPNN: XInd++; break;
	case cPNP: XInd++; ZInd++; break;
	case cPPN: XInd++; YInd++; break;
	case cPPP: XInd++; YInd++; ZInd++; break;
	}
	return D3Ind(XInd, YInd, ZInd, XSize+1, YSize+1);
}

void CVX_MeshUtil::GetCornerDir(int Corner, Vec3D<>* pOut) //returns vector with component +/- 1 dependeind on which corner
{
	switch (Corner){
	case cNNN: *pOut = Vec3D<>(-1,-1,-1); break;
	case cNNP: *pOut = Vec3D<>(-1,-1, 1); break;
	case cNPN: *pOut = Vec3D<>(-1, 1,-1); break;
	case cNPP: *pOut = Vec3D<>(-1, 1, 1); break;
	case cPNN: *pOut = Vec3D<>( 1,-1,-1); break;
	case cPNP: *pOut = Vec3D<>( 1,-1, 1); break;
	case cPPN: *pOut = Vec3D<>( 1, 1,-1); break;
	case cPPP: *pOut = Vec3D<>( 1, 1, 1); break;
	default: *pOut = Vec3D<>(0,0,0);
	}
}

void CVX_MeshUtil::Draw(void)
{
#ifdef USE_OPEN_GL
	DefMesh.Draw();
#endif
}


//Misc
bool CVX_MeshUtil::ToStl(std::string BasePath, CVX_Object* pObj, bool WantDefMes)
{
	if (BasePath == "") return false;

	CMesh tmpMesh;
	int XSize = pObj->GetVXDim();
	int YSize = pObj->GetVYDim();

	switch (pObj->Voxel.GetVoxName()){
	case VS_BOX:
		for (int m=1; m<=(int)pObj->Palette.size(); m++){ //for each material (skip null material)
			if (pObj->GetNumVox(m) == 0) continue; //see if there's material here, continue if not

			tmpMesh.Clear();
			int CheckIndex = 0;
			Vec3D<> Center;
			int NomX, NomY, NomZ;
			Vec3D<> LatDims = pObj->GetLatDimEnv();
			//vertices of the cube in order of (---, --+, -+-, -++, +--, +-+, ++-, +++)
			Vec3D<> V[8] = {Vec3D<>(0,0,0)};

			for (int i=0; i<pObj->GetStArraySize(); i++){ //for every element of the master array
				if (pObj->GetMat(i) != m) continue;  //if there's not a voxel here of the correct material, keep going

				pObj->GetXYZ(&Center, i); //loads the center variable with coordinates of center of cube
				pObj->GetXYZNom(&NomX, &NomY, &NomZ, i);

				if (WantDefMes){
					GetCurVLoc(CalcVertsAll[D3IndCorner(NomX, NomY, NomZ, XSize, YSize, cNNN)], &V[cNNN]);
					GetCurVLoc(CalcVertsAll[D3IndCorner(NomX, NomY, NomZ, XSize, YSize, cNNP)], &V[cNNP]);
					GetCurVLoc(CalcVertsAll[D3IndCorner(NomX, NomY, NomZ, XSize, YSize, cNPN)], &V[cNPN]);
					GetCurVLoc(CalcVertsAll[D3IndCorner(NomX, NomY, NomZ, XSize, YSize, cNPP)], &V[cNPP]);
					GetCurVLoc(CalcVertsAll[D3IndCorner(NomX, NomY, NomZ, XSize, YSize, cPNN)], &V[cPNN]);
					GetCurVLoc(CalcVertsAll[D3IndCorner(NomX, NomY, NomZ, XSize, YSize, cPNP)], &V[cPNP]);
					GetCurVLoc(CalcVertsAll[D3IndCorner(NomX, NomY, NomZ, XSize, YSize, cPPN)], &V[cPPN]);
					GetCurVLoc(CalcVertsAll[D3IndCorner(NomX, NomY, NomZ, XSize, YSize, cPPP)], &V[cPPP]);

					//V[cNNN] = 
				}
				else{ //not deformed mesh
				V[cNNN] = Center+Vec3D<>(-LatDims.x/2,-LatDims.y/2,-LatDims.z/2);
				V[cNNP] = Center+Vec3D<>(-LatDims.x/2,-LatDims.y/2, LatDims.z/2);
				V[cNPN] = Center+Vec3D<>(-LatDims.x/2, LatDims.y/2,-LatDims.z/2);
				V[cNPP] = Center+Vec3D<>(-LatDims.x/2, LatDims.y/2, LatDims.z/2);
				V[cPNN] = Center+Vec3D<>( LatDims.x/2,-LatDims.y/2,-LatDims.z/2);
				V[cPNP] = Center+Vec3D<>( LatDims.x/2,-LatDims.y/2, LatDims.z/2);
				V[cPPN] = Center+Vec3D<>( LatDims.x/2, LatDims.y/2,-LatDims.z/2);
				V[cPPP] = Center+Vec3D<>( LatDims.x/2, LatDims.y/2, LatDims.z/2);
				}

				for (int j=0; j<6; j++){ //for each face of the cube
					switch (j){
						case 0: //top! (+Z)
							CheckIndex = pObj->GetIndex(NomX, NomY, NomZ+1);
							if (CheckIndex == -1 || pObj->GetMat(CheckIndex) != m) //if its outside range or there's no cube there...
								tmpMesh.AddQuadFacet(V[cPPP], V[cNPP], V[cNNP], V[cPNP]);
						break;
						case 1: //bottom!
							CheckIndex = pObj->GetIndex(NomX, NomY, NomZ-1);
							if (CheckIndex == -1 || pObj->GetMat(CheckIndex) != m)
								tmpMesh.AddQuadFacet(V[cPNN], V[cNNN], V[cNPN], V[cPPN]);
						break;
						case 2: //right!
							CheckIndex = pObj->GetIndex(NomX+1, NomY, NomZ);
							if (CheckIndex == -1 || pObj->GetMat(CheckIndex) != m)
								tmpMesh.AddQuadFacet(V[cPNN], V[cPPN], V[cPPP], V[cPNP]);
						break;

						case 3: //left!
							CheckIndex = pObj->GetIndex(NomX-1, NomY, NomZ);
							if (CheckIndex == -1 || pObj->GetMat(CheckIndex) != m)
								tmpMesh.AddQuadFacet(V[cNNN], V[cNNP], V[cNPP], V[cNPN]);
						break;

						case 4: //front!
							CheckIndex = pObj->GetIndex(NomX, NomY+1, NomZ);
							if (CheckIndex == -1 || pObj->GetMat(CheckIndex) != m)
								tmpMesh.AddQuadFacet(V[cNPN], V[cNPP], V[cPPP], V[cPPN]);
						break;

						case 5: //back!
							CheckIndex = pObj->GetIndex(NomX, NomY-1, NomZ);
							if (CheckIndex == -1 || pObj->GetMat(CheckIndex) != m)
								tmpMesh.AddQuadFacet(V[cPNN], V[cPNP], V[cNNP], V[cNNN]);
						break;
					}
				}
			}
			std::string tmp = "_" + pObj->Palette[m].GetName();
			std::string ThisPath = BasePath;
			ThisPath.insert(ThisPath.size()-4, tmp);
			tmpMesh.SaveSTL(ThisPath);
		}
		
		break;
	case VS_SPHERE:
		for (int m=1; m<=(int)pObj->Palette.size(); m++){ //for each material (skip null material)
			if (pObj->GetNumVox(m) == 0) continue; //see if there's material here, continue if not

			tmpMesh.Clear();
			Vec3D<> Center, p1, p2, p3, p4;
			vfloat a1, a2, x1, x2, Rad1, Rad2;

			for (int i=0; i<pObj->GetStArraySize(); i++){ //for every element of the master array
				if (pObj->GetMat(i) != m) continue;  //if there's not a voxel here of the correct material, keep going
				pObj->GetXYZ(&Center, i); //loads the center variable with coordinates of center of cube
				
				//adds a sphere to the mesh:
				vfloat pi = atan(1.0)*4.0;
				vfloat da = pi/12;
				vfloat Scale = pObj->GetLatticeDim()/2.01;

				for (a2=0; a2<pi; a2+=da) {
					x1 = cos(a2);
					Rad1 = sin(a2);
					x2 = cos(a2+da);
					Rad2 = sin(a2+da);
					for (a1=0; a1<pi*2.01; a1+=da) {
						p1 = Vec3D<>(x1, Rad1*cos(a1), Rad1*sin(a1));
						p2 = Vec3D<>(x2, Rad2*cos(a1), Rad2*sin(a1));
						p3 = Vec3D<>(x1, Rad1*cos(a1+da), Rad1*sin(a1+da));
						p4 = Vec3D<>(x2, Rad2*cos(a1+da), Rad2*sin(a1+da));

						p1 = p1*Scale + Center;
						p2 = p2*Scale + Center;
						p3 = p3*Scale + Center;
						p4 = p4*Scale + Center;


						if (Rad1<1e-6)
							tmpMesh.AddFacet(p1, p2, p4);
						else if (Rad2<1e-6)
							tmpMesh.AddFacet(p1, p2, p3);
						else {
							tmpMesh.AddFacet(p1, p2, p3);
							tmpMesh.AddFacet(p2, p4, p3);
						}

					}
				}

			}


			std::string tmp = "_" + pObj->Palette[m].GetName();
			std::string ThisPath = BasePath;
			ThisPath.insert(ThisPath.size()-4, tmp);
			tmpMesh.SaveSTL(ThisPath);
		}
		break;

	default:
		return false;
	}
	return true;
}



bool CVX_MeshUtil::FromStl(CMesh* pMeshIn, CVX_Object* pObj, int MatIndex)
{
	Vec3D<> Loc;
	for (int i=0; i<pObj->GetStArraySize(); i++){
		Loc = pObj->GetXYZ(i);
		if(pMeshIn->IsInside(&Loc)) pObj->SetMat(i, MatIndex);
	}

	return true;

	/*
//	if (pObj->Structure.X_Voxels == 0)
//		pObj->InitializeMatter(); //initialize defaults if haven't loaded yet...

	Vec3D Point;
	//pObj->Cl
	int TotNumFacets = (int)pMeshIn->Facets.size();
	int* pZList = new int[TotNumFacets]; //array of facet indices that touch current Z plane...
	for (int i=0; i<TotNumFacets; i++) pZList[i] = -1; //initialize all to -1
	
	vfloat* pIntrscts = new vfloat[1000]; //the most intersections we expect to see...

	int Xv, Yv, Zv;
	int CurNumInt; //current number of intersections...

	pObj->GetVDim(&Xv, &Yv, &Zv);
	for (int k=0; k<Zv; k++){ //iterate through z
		pObj->GetXYZ(&Point, 0, 0, k); //this height
		vfloat ZHeight = Point.z;

		//clear the previous list:
		int iter = 0;
		while (pZList[iter] != -1 && iter < TotNumFacets)
			pZList[iter++] = -1;

		//add any Facets whose Z coordinates are not all above or all below this Z plane
		int NumZFacets = 0;
		bool IsAbove, IsBelow;

		//TRACE("ZHeight: %g\n", ZHeight);
		for (int m=0; m<TotNumFacets; m++){
			IsAbove = true; IsBelow = true;

			for (int n=0; n<3; n++){
				if(pMeshIn->Vertices[pMeshIn->Facets[m].vi[n]].v.z > ZHeight) IsBelow = false;
				if(pMeshIn->Vertices[pMeshIn->Facets[m].vi[n]].v.z < ZHeight) IsAbove = false;
			}
			if (!IsAbove && !IsBelow){ //if this facet is not fully above or fully below our ZPlane
				pZList[NumZFacets++] = m; //add to out current ZList...
				//TRACE("V1.z = %g, V2.z = %g, V3.z = %g\n", pModel->Vertices[pModel->Facets[m].vi[0]].v.z, pModel->Vertices[pModel->Facets[m].vi[1]].v.z, pModel->Vertices[pModel->Facets[m].vi[2]].v.z);
			}
		}

		
		for (int j=0; j<Yv; j++){ //iterate through y
			pObj->GetXYZ(&Point, 0, j, k); //this location
			CurNumInt = pMeshIn->GetXIntersections(Point.z, Point.y, pIntrscts, NumZFacets, pZList); //fill in pIntrscts

			int IntrInd = 0;
			for (int i=0; i<Xv; i++){ //iterate through x
				pObj->GetXYZ(&Point, i, j, k); //this location

				while (pIntrscts[IntrInd] < Point.x && IntrInd < CurNumInt) //step through our array of intersections until
					IntrInd++;

				if (IntrInd%2 == 1) pObj->SetVoxel(i, j, k, MatIndex); //if we've gone through an odd number of intersections, put material here
				//else pObj->SetVoxel(i, j, k, 0); //otherwise empty...

			}
		}
	}

	delete [] pZList;
	pZList = NULL;
	delete [] pIntrscts;
	pIntrscts = NULL;

	return true;
*/
}
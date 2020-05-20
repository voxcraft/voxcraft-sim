#ifdef _0
#include "TI_Collision.h"

#define ENVELOPE_RADIUS 0.625f

CUDA_DEVICE TI_Collision::TI_Collision(TI_Voxel* v1, TI_Voxel* v2)
{
	pV1 = v1;
    pV2 = v2;
    v1->material();
	penetrationStiff = 2.0f/(1.0f/v1->material()->penetrationStiffness()+1.0f/pV2->material()->penetrationStiffness());
	dampingC = 0.5f*(v1->material()->collisionDampingTranslateC() + v2->material()->collisionDampingTranslateC()); //average
}

CUDA_DEVICE TI_Collision& TI_Collision::operator=(const TI_Collision& col)
{
	pV1 = col.pV1;
	pV2 = col.pV2;
	penetrationStiff = col.penetrationStiff;
	force = col.force;
	return *this;
}

CUDA_DEVICE TI_Vec3D<float> const TI_Collision::contactForce(TI_Voxel* pVoxel)
{
	if (pVoxel == pV1) return force;
	else if (pVoxel == pV2) return -force;
	else return TI_Vec3D<float>(0,0,0);
}


CUDA_DEVICE void TI_Collision::updateContactForce() 
{
	//just basic sphere envelope, repel with the stiffness of the material... (assumes UpdateConstants has been called)
	TI_Vec3D<float> offset = (TI_Vec3D<float>)(pV2->position() - pV1->position());
	float NomDist = (float)((pV1->baseSizeAverage() + pV2->baseSizeAverage())*ENVELOPE_RADIUS); //effective diameter of 1.5 voxels... (todo: remove length2!!
	float RelDist = NomDist -offset.Length(); //negative for overlap!

	if (RelDist > 0){
		TI_Vec3D<float> unit = offset.Normalized(); //unit vector from voxel 1 in the direction of voxel 2
		float relativeVelocity = pV1->velocity().Dot((TI_Vec3D<double>)unit) - pV2->velocity().Dot((TI_Vec3D<double>)unit); //negative is moving towards each other

		force = unit * (penetrationStiff * RelDist + dampingC*relativeVelocity); //kx + cV - if we're overlapping
	}
	else force = TI_Vec3D<float>(0,0,0);
}

#endif
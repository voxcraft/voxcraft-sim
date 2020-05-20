#include "VX3_Collision.h"

__device__ VX3_Collision::VX3_Collision(VX3_Voxel* v1, VX3_Voxel* v2) {
    pV1 = v1;
    pV2 = v2;
	penetrationStiff = 2.0f/(1.0f/v1->material()->penetrationStiffness()+1.0f/v2->material()->penetrationStiffness());
	dampingC = 0.5f*(v1->material()->collisionDampingTranslateC() + v2->material()->collisionDampingTranslateC()); //average
}

__device__ VX3_Vec3D<double> const VX3_Collision::contactForce(VX3_Voxel* pVoxel)
{
	if (pVoxel == pV1) return force;
	else if (pVoxel == pV2) return -force;
	else return VX3_Vec3D<double>(0,0,0);
}

__device__ void VX3_Collision::updateContactForce()
{
	//just basic sphere envelope, repel with the stiffness of the material... (assumes UpdateConstants has been called)
	VX3_Vec3D<double> offset = (VX3_Vec3D<double>)(pV2->position() - pV1->position());
	double NomDist = (double)((pV1->baseSizeAverage() + pV2->baseSizeAverage())*COLLISION_ENVELOPE_RADIUS); //effective diameter of 1.5 voxels... (todo: remove length2!!
	double RelDist = NomDist -offset.Length(); //negative for overlap!

	if (RelDist > 0){ //Why RelDist > 0 and there's penetration?
		VX3_Vec3D<double> unit = offset.Normalized(); //unit vector from voxel 1 in the direction of voxel 2
		double relativeVelocity = pV1->velocity().Dot((VX3_Vec3D<double>)unit) - pV2->velocity().Dot((VX3_Vec3D<double>)unit); //negative is moving towards each other
		// printf("%f, %f, %f;  %f * %f + %f * %f\n", offset.Length(), pV1->baseSizeAverage(), NomDist, penetrationStiff, RelDist, dampingC, relativeVelocity);
		force = unit * (penetrationStiff * RelDist + dampingC*relativeVelocity); //kx + cV - if we're overlapping
	}
	else force = VX3_Vec3D<double>(0,0,0);
}

#if !defined(TI_COLLISION_H)
#define TI_COLLISION_H

#include "TI_Utils.h"

#include "TI_Voxel.h"

class TI_Collision
{
public:
	CUDA_DEVICE TI_Collision(TI_Voxel* v1, TI_Voxel* v2); //!< Constructor taking the two voxels to watch for collision between. The order is irrelevant. @param[in] v1 One voxel @param[in] v2 The other voxel
	CUDA_DEVICE TI_Collision& operator=(const TI_Collision& col); //!< Overload "=" operator.
	CUDA_DEVICE TI_Collision(const TI_Collision& col) {*this = col;} //!< copy constructor.

	CUDA_DEVICE TI_Vec3D<float> const contactForce(TI_Voxel* pVoxel); //!< Returns the repelling force acting on pVoxel from the penetration of the other voxel if their collision boundaries overlap. Otherwise returns a zero force. This force will only be accurate if updateContactForce() has been called since the voxels last moved. @param[in] pVoxel The voxel in question. This should be voxel1() or voxel2() to have any meaning. Otherwise a zero force is returned.
	CUDA_DEVICE void updateContactForce(); //!< Updates the state this collision based on the current positions and properties of voxel1() and voxel2(). Call contactForce() with either voxel as the argument to obtain the repelling penetration force (if any) that exisits.

	CUDA_DEVICE TI_Voxel* voxel1() const {return pV1;} //!<One voxel of this potential collision pair.
	CUDA_DEVICE TI_Voxel* voxel2() const {return pV2;} //!<The other voxel of this potential collision pair.

/* data */
	//static float envelopeRadius; //!<The collision envelope radius that these two voxels collide at. Even though voxels are cubic, a spherical collision envelope is used for computation efficiency. Values are multiplied by the length of an edge of the voxel to determine the actual collision radius. Values less than 0.5 or greater than 0.866 are probably of limited use. Prefer around 0.625 (default).

	TI_Voxel *pV1, *pV2;
	float penetrationStiff; //in N/m for these two voxels
	float dampingC; //damping factor for these two voxels
	TI_Vec3D<float> force;
};

#endif // TI_COLLISION_H

#if !defined(VX3_MATERIAL_LINK_H)
#define VX3_MATERIAL_LINK_H

#include "VX3.cuh"
#include "VX_MaterialLink.h" // for CVX_MaterialLink

#include "VX3_MaterialVoxel.h"
class VX3_VoxelyzeKernel;

class VX3_MaterialLink : public VX3_MaterialVoxel {
public:
	VX3_MaterialLink(CVX_MaterialLink* p, VX3_VoxelyzeKernel* k);
	~VX3_MaterialLink();

	__device__ VX3_MaterialLink(VX3_MaterialVoxel* mat1, VX3_MaterialVoxel* mat2); //!< Creates a link material from the two specified voxel materials. The order is unimportant. @param[in] mat1 voxel material on one side of the link. @param[in] mat2 voxel material on the other side of the link.
	__device__ VX3_MaterialLink(const VX3_MaterialLink& VIn) {*this = VIn;} //!< Copy constructor

	__device__ VX3_MaterialLink& operator=(const VX3_MaterialLink& VIn); //!< Equals operator

	__device__ virtual bool updateAll(); //!< Updates and recalculates eveything possible (used by inherited classed when material properties have changed)
	__device__ virtual bool updateDerived(); //!< Updates all the derived quantities cached as member variables for this and derived classes. (Especially if density, size or elastic modulus changes.)

/* data */
	VX3_MaterialVoxel *vox1Mat = NULL; //!< Constituent material 1 from one voxel
	VX3_MaterialVoxel *vox2Mat = NULL; //!< Constituent material 2 from the other voxel

	float _a1; //!< Cached a1 beam constant.
	float _a2; //!< Cached a2 beam constant.
	float _b1; //!< Cached b1 beam constant.
	float _b2; //!< Cached b2 beam constant.
	float _b3; //!< Cached b3 beam constant.
	float _sqA1; //!< Cached sqrt(a1) constant for damping calculations.
	float _sqA2xIp; //!< Cached sqrt(a2*L*L/6) constant for damping calculations.
	float _sqB1; //!< Cached sqrt(b1) constant for damping calculations.
	float _sqB2xFMp; //!< Cached sqrt(b2*L/2) constant for damping calculations.
	float _sqB3xIp; //!< Cached sqrt(b3*L*L/6) constant for damping calculations.
	
};

#endif // VX3_MATERIAL_LINK_H

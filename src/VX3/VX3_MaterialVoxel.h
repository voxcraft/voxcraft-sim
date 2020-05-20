#if !defined(VX3_MATERIALVOXEL_H)
#define VX3_MATERIALVOXEL_H

#include "VX3.cuh"
#include "VX_MaterialVoxel.h" // for CVX_MaterialVoxel
#include "VX3_Material.h"
class VX3_VoxelyzeKernel;

class VX3_MaterialVoxel : public VX3_Material {
public:
	VX3_MaterialVoxel( CVX_MaterialVoxel *p, VX3_VoxelyzeKernel* k);

	__device__ VX3_MaterialVoxel(float youngsModulus=1e6f, float density=1e3f, double nominalSize=0.001); //!< Default Constructor. @param[in] youngsModulus The Young's Modulus (stiffness) of this material in Pascals. @param[in] density The density of this material in Kg/m^3. @param[in] nominalSize The nominal voxel size in meters.
	__device__ VX3_MaterialVoxel(const VX3_Material& mat, double nominalSize=0.001); //!< Constructs from an existing material. @param[in] mat Material to construct from. @param[in] nominalSize The nominal voxel size in meters
	//virtual ~VX3_MaterialVoxel(void); //!< Destructor. Virtual so we can just keep track of VX3_Material pointers.
	__device__ VX3_MaterialVoxel(const VX3_MaterialVoxel& vIn) {*this = vIn;} //!< Copy constructor
	__device__ VX3_MaterialVoxel& operator=(const VX3_MaterialVoxel& vIn); //!< Equals operator

	//size and scaling
	__device__ bool setNominalSize(double size); //!< Sets the nominal cubic voxel size in order to calculate mass, moments of inertia, etc of this material. In ordinary circumstances this should be controlled by the overall simulation and never called directly. Use setExternalScaleFactor() if you wish to change the size of voxels of this material. @param[in] size The size of a voxel as measured by its linear outer dimension. (units: meters)
	__device__ double nominalSize(){return nomSize;} //!< Returns the nominal cubic voxel size in meters.
	__device__ VX3_Vec3D<double> size() {return nomSize*extScale;} //!< Returns the current nominal size (in meters) of any voxels of this material-including external scaling factors. The size is calculated according to baseSize()*externalScaleFactor(). This represents the nominal size for voxels of this material, and every instantiated voxel may have a different actual size based on the forces acting upon it in context. This value does not include thermal expansions and contractions which may also change the nominal size of a given voxel depending on its CTE and current temperature.

	//mass and inertia
	__device__ float mass(){return _mass;} //!<Returns the mass of a voxel of this material in Kg. Mass cannot be specified directly. Mass is indirectly calculated according to setDensity() and setBaseSize().
	__device__ float momentInertia(){return _momentInertia;} //!<Returns the first moment of inertia of a voxel of this material in kg*m^2. This quantity is indirectly calculated according to setDensity() and setBaseSize().

	//damping convenience functions
	__device__ float internalDampingTranslateC() const {return zetaInternal*_2xSqMxExS;} //!< Returns the internal material damping coefficient (translation).
	__device__ float internalDampingRotateC() const {return zetaInternal*_2xSqIxExSxSxS;} //!< Returns the internal material damping coefficient (rotation).
	__device__ float globalDampingTranslateC() const {return zetaGlobal*_2xSqMxExS;} //!< Returns the global material damping coefficient (translation)
	__device__ float globalDampingRotateC() const {return zetaGlobal*_2xSqIxExSxSxS;} //!< Returns the global material damping coefficient (rotation)
	__device__ float collisionDampingTranslateC() const {return zetaCollision*_2xSqMxExS;} //!< Returns the global material damping coefficient (translation)
	__device__ float collisionDampingRotateC() const {return zetaCollision*_2xSqIxExSxSxS;} //!< Returns the global material damping coefficient (rotation)

	//stiffness
	__device__ float penetrationStiffness() const {return (float)(2*E*nomSize);} //!< returns the stiffness with which this voxel will resist penetration. This is calculated according to E*A/L with L = voxelSize/2.

	__device__ void initialize(double nominalSize); //!< Initializes this voxel material with the specified voxel size. @param[in] nominalSize The nominal voxel size in meters.
	
	//only the main simulation should update gravity
	__device__ void setGravityMultiplier(float gravityMultiplier){gravMult = gravityMultiplier;} //!< Sets the multiple of gravity for this material. In normal circumstances this should only be called by the parent voxelyze object. @param[in] gravityMultiplier Gravity multiplier (1 = 1G gravity).
	__device__ float gravityMuliplier() {return gravMult;} //!< Returns the current gravity multiplier.
	__device__ float gravityForce() {return -_mass*9.80665f*gravMult;} //!< Returns the current gravitational force on this voxel according to F=ma.

	__device__ virtual bool updateAll() {return false;} //!< Updates and recalculates eveything possible (used by inherited classed when material properties have changed)
	__device__ virtual bool updateDerived(); //!< Updates all the derived quantities cached as member variables for this and derived classes. (Especially if density, size or elastic modulus changes.)

/* data */
	double nomSize; //!< Nominal size (i.e. lattice dimension) (m)
	float gravMult; //!< Multiplier for how strongly gravity should affect this block in g (1.0 = -9.81m/s^2)
	float _mass; //!< Cached mass of this voxel (kg)
	float _massInverse; //!< Cached 1/Mass (1/kg)
	float _sqrtMass; //!< Cached sqrt(mass). (sqrt(Kg))
	float _firstMoment; //!< Cached 1st moment "inertia" (needed for certain calculations) (kg*m)
	float _momentInertia; //!< Cached mass moment of inertia (i.e. rotational "mass") (kg*m^2)
	float _momentInertiaInverse; //!< Cached 1/Inertia (1/(kg*m^2))
	float _2xSqMxExS; //!< Cached value needed for quick damping calculations (Kg*m/s)
	float _2xSqIxExSxSxS; //!< Cached value needed for quick rotational damping calculations (Kg*m^2/s)
};


#endif // VX3_MATERIALVOXEL_H

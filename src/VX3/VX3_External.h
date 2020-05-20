#if !defined(VX3_EXTERNAL_H)
#define VX3_EXTERNAL_H

#include "VX3.cuh"
#include "VX_External.h"

class VX3_External {
public:
	VX3_External( CVX_External* p );

	__device__ inline void dofSet(dofObject& obj, dofComponent dof, bool set) {set ? obj|=dof : obj&=~dof;}
	__device__ inline void dofSetAll(dofObject& obj, bool set) {set ? obj|=0x3F : obj&=~0x3F;}
	__device__ inline bool dofIsSet(dofObject obj, dofComponent dof){return (dof&obj)?true:false;}
	__device__ inline bool dofIsAllSet(dofObject obj){return (obj&0x3F)==0x3F;}
	__device__ inline bool dofIsNoneSet(dofObject obj){return !(obj&0x3F);}
	__device__ inline dofObject dof(bool tx, bool ty, bool tz, bool rx, bool ry, bool rz) {dofObject ret=0; dofSet(ret, X_TRANSLATE, tx); dofSet(ret, Y_TRANSLATE, ty); dofSet(ret, Z_TRANSLATE, tz); dofSet(ret, X_ROTATE, rx); dofSet(ret, Y_ROTATE, ry); dofSet(ret, Z_ROTATE, rz); return ret;}

	__device__ VX3_External();
	__device__ VX3_External(const VX3_External& eIn) {*this = eIn;} //!< Copy constructor
	__device__ VX3_External& operator=(const VX3_External& eIn); //!< Equals operator
	__device__ inline bool operator==(const VX3_External& b) {return dofFixed==b.dofFixed && extForce==b.extForce && extMoment==b.extMoment && extTranslation==b.extTranslation && extRotation==b.extRotation;} //!< comparison operator

	__device__ void reset(); //!< Resets this external to defaults - i.e., no effect on a voxel  (forces, fixed, displacements, etc) 
	__device__ bool isEmpty() {return (dofFixed == 0 && extForce==VX3_Vec3D<float>() && extMoment==VX3_Vec3D<float>());} //!< returns true if this external is empty - i.e is exerting no effect on a voxel

	__device__ bool isFixed(dofComponent dof) {return dofIsSet(dofFixed, dof);}  //!< Returns true if the specified degree of freedom is fixed for this voxel. @param[in] dof Degree of freedom to query according to the dofComponent enum.
	__device__ bool isFixedAll() {return dofIsAllSet(dofFixed);} //!< Returns true if all 6 degrees of freedom are fixed for this voxel.
	__device__ bool isFixedAllTranslation() {return dofIsSet(dofFixed, X_TRANSLATE) && dofIsSet(dofFixed, Y_TRANSLATE) && dofIsSet(dofFixed, Z_TRANSLATE);} //!< Returns true if all translational degrees of freedom are fixed.
	__device__ bool isFixedAllRotation() {return dofIsSet(dofFixed, X_ROTATE) && dofIsSet(dofFixed, Y_ROTATE) && dofIsSet(dofFixed, Z_ROTATE);} //!< Returns true if all rotationsl degrees of freedom are fixed.
	
	__device__ bool isFixedAny() {return (dofFixed != 0);} //!< Returns true if any of the 6 degrees of freedom are fixed for this voxel.
	__device__ bool isFixedAnyTranslation() {return dofIsSet(dofFixed, X_TRANSLATE) || dofIsSet(dofFixed, Y_TRANSLATE) || dofIsSet(dofFixed, Z_TRANSLATE);} //!< Returns true if any of the three translational degrees of freedom are fixed.
	__device__ bool isFixedAnyRotation() {return dofIsSet(dofFixed, X_ROTATE) || dofIsSet(dofFixed, Y_ROTATE) || dofIsSet(dofFixed, Z_ROTATE);} //!< Returns true if any of the three rotational degrees of freedom are fixed.

	__device__ VX3_Vec3D<double> translation() {return extTranslation;} //!< Returns any external translation that has been applied to this external.
	__device__ VX3_Vec3D<double> rotation() {return extRotation;} //!< Returns any external rotation that has been applied to this external as a rotation vector.
	__device__ VX3_Quat3D<double> rotationQuat() {return _extRotationQ;} //!< Returns any external rotation that has been applied to this external as a quaternion.


	__device__ void setFixed(bool xTranslate, bool yTranslate, bool zTranslate, bool xRotate, bool yRotate, bool zRotate); //!< Sets any of the degrees of freedom specified as "true" to fixed for this voxel. (GCS) @param[in] xTranslate Translation in the X direction  @param[in] yTranslate Translation in the Y direction @param[in] zTranslate Translation in the Z direction @param[in] xRotate Rotation about the X axis @param[in] yRotate Rotation about the Y axis @param[in] zRotate Rotation about the Z axis
	__device__ void setFixed(dofComponent dof, bool fixed=true) {fixed?setDisplacement(dof):clearDisplacement(dof);} //!< Sets the specified degree of freedom to either fixed or free, depending on the value of fixed. Either way, sets the translational or rotational displacement of this degree of freedom to zero. @param[in] dof the degree of freedom in question @param[in] fixed Whether this degree of freedom should be fixed (true) or free (false).
	__device__ void setFixedAll(bool fixed=true) {fixed?setDisplacementAll():clearDisplacementAll();} //!< Sets all 6 degrees of freedom to either fixed or free depending on the value of fixed. Either way, sets all displacements to zero. @param[in] fixed Whether all degrees of freedom should be fixed (true) or free (false).

	__device__ void setDisplacement(dofComponent dof, double displacement=0.0); //!< Fixes the specified degree of freedom and applies the prescribed displacement if specified. @param[in] dof the degree of freedom in question. @param[in] displacement The displacement in meters (translational dofs) or radians (rotational dofs) to apply. Large fixed displacements may cause instability.
	__device__ void setDisplacementAll(const VX3_Vec3D<double>& translation = VX3_Vec3D<double>(0,0,0), const VX3_Vec3D<double>& rotation = VX3_Vec3D<double>(0,0,0)); //!< Fixes the all degrees of freedom and applies the specified translation and rotation. @param[in] translation The translation in meters from this voxel's nominal position to fix it at. @param[in] rotation The rotation (in the form of a rotation vector) from this voxel's nominal orientation to fix it at.

	__device__ void clearDisplacement(dofComponent dof); //!< Clears any prescribed displacement from this degree of freedom and unfixes it, too. @param[in] dof the degree of freedom in question.
	__device__ void clearDisplacementAll(); //!< Clears all prescribed displacement from this voxel and completely unfixes it, too.

	__device__ VX3_Vec3D<float> force() {return extForce;} //!< Returns the current applied external force in newtons.
	__device__ VX3_Vec3D<float> moment() {return extMoment;} //!< Returns the current applied external moment in N-m.

	__device__ void setForce(const float xForce, const float yForce, const float zForce) {extForce = VX3_Vec3D<float>(xForce, yForce, zForce);} //!< Applies forces to this voxel in the global coordinate system. Has no effect in any fixed degrees of freedom. @param xForce Force in the X direction in newtons.  @param yForce Force in the Y direction in newtons.  @param zForce Force in the Z direction in newtons. 
	__device__ void setForce(const VX3_Vec3D<float>& force) {extForce = force;} //!< Convenience function for setExternalForce(float, float, float).
	__device__ void setMoment(const float xMoment, const float yMoment, const float zMoment) {extMoment = VX3_Vec3D<float>(xMoment, yMoment, zMoment);}  //!< Applies moments to this voxel in the global coordinate system. All rotations according to the right-hand rule. Has no effect in any fixed degrees of freedom. @param xMoment Moment in the X axis rotation in newton-meters. @param yMoment Moment in the Y axis rotation in newton-meters. @param zMoment Moment in the Z axis rotation in newton-meters. 
	__device__ void setMoment(const VX3_Vec3D<float>& moment) {extMoment = moment;} //!< Convenience function for setExternalMoment(float, float, float).

	__device__ void addForce(const float xForce, const float yForce, const float zForce) {extForce += VX3_Vec3D<float>(xForce, yForce, zForce);} //!< Applies forces to this voxel in the global coordinate system. Has no effect in any fixed degrees of freedom. @param xForce Force in the X direction in newtons.  @param yForce Force in the Y direction in newtons.  @param zForce Force in the Z direction in newtons. 
	__device__ void addForce(const VX3_Vec3D<float>& force) {extForce += force;} //!< Convenience function for setExternalForce(float, float, float).
	__device__ void addMoment(const float xMoment, const float yMoment, const float zMoment) {extMoment += VX3_Vec3D<float>(xMoment, yMoment, zMoment);}  //!< Applies moments to this voxel in the global coordinate system. All rotations according to the right-hand rule. Has no effect in any fixed degrees of freedom. @param xMoment Moment in the X axis rotation in newton-meters. @param yMoment Moment in the Y axis rotation in newton-meters. @param zMoment Moment in the Z axis rotation in newton-meters. 
	__device__ void addMoment(const VX3_Vec3D<float>& moment) {extMoment += moment;} //!< Convenience function for setExternalMoment(float, float, float).

	__device__ void clearForce(){extForce = VX3_Vec3D<float>();} //!< Clears all applied forces from this voxel.
	__device__ void clearMoment(){extMoment = VX3_Vec3D<float>();} //!< Clears all applied moments from this voxel.

	__device__ void rotationChanged(); //called to keep cached quaternion rotation in sync

/* data */	
	dofObject dofFixed;
	
	VX3_Vec3D<float> extForce, extMoment; //External force, moment applied to these voxels (N, N-m) if relevant DOF are unfixed
	VX3_Vec3D<double> extTranslation, extRotation;
	VX3_Quat3D<double> _extRotationQ; //cached quaternion rotation (pointer to only create if needed)

};

#endif // VX3_EXTERNAL_H

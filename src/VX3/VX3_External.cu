#include "VX3_External.h"

VX3_External::VX3_External( CVX_External* p ) :
dofFixed(p->dofFixed), extForce(p->extForce), extMoment(p->extMoment),
extTranslation(p->extTranslation), extRotation(p->extRotation), 
_extRotationQ(p->_extRotationQ) {

}

__device__ VX3_External::VX3_External() 
{
	reset();
}

__device__ VX3_External& VX3_External::operator=(const VX3_External& eIn)
{
	dofFixed = eIn.dofFixed;
	extForce = eIn.extForce;
	extMoment = eIn.extMoment;
	extTranslation = eIn.extTranslation;
	extRotation = eIn.extRotation;
	rotationChanged();
	return *this;
}

__device__ void VX3_External::reset()
{
	dofFixed=0;
	extForce = extMoment = VX3_Vec3D<float>();
	extTranslation = VX3_Vec3D<double>();
	extRotation = VX3_Vec3D<double>();
	rotationChanged();
}


__device__ void VX3_External::setFixed(bool xTranslate, bool yTranslate, bool zTranslate, bool xRotate, bool yRotate, bool zRotate)
{
	dofFixed = dof(xTranslate, yTranslate, zTranslate, xRotate, yRotate, zRotate);
	extTranslation = extRotation = VX3_Vec3D<double>(); //clear displacements
}

__device__ void VX3_External::setDisplacement(dofComponent dof, double displacement)
{
	dofSet(dofFixed, dof, true);
	if (displacement != 0.0f){
		if (dof & X_TRANSLATE) extTranslation.x = displacement;
		if (dof & Y_TRANSLATE) extTranslation.y = displacement;
		if (dof & Z_TRANSLATE) extTranslation.z = displacement;
		if (dof & X_ROTATE) extRotation.x = displacement;
		if (dof & Y_ROTATE)	extRotation.y = displacement;
		if (dof & Z_ROTATE) extRotation.z = displacement;
	}

	rotationChanged();
}

__device__ void VX3_External::setDisplacementAll(const VX3_Vec3D<double>& translation, const VX3_Vec3D<double>& rotation)
{
	dofSetAll(dofFixed, true);
	extTranslation = translation;
	extRotation = rotation;

	rotationChanged();
}

__device__ void VX3_External::clearDisplacement(dofComponent dof)
{
	dofSet(dofFixed, dof, false);

	if (dof & X_TRANSLATE) extTranslation.x = 0.0;
	if (dof & Y_TRANSLATE) extTranslation.y = 0.0;
	if (dof & Z_TRANSLATE) extTranslation.z = 0.0;
	if (dof & X_ROTATE) extRotation.x = 0.0;
	if (dof & Y_ROTATE)	extRotation.y = 0.0;
	if (dof & Z_ROTATE) extRotation.z = 0.0;

	rotationChanged();
}

__device__ void VX3_External::clearDisplacementAll()
{
	dofSetAll(dofFixed, false);
	extTranslation = VX3_Vec3D<double>();
	extRotation = VX3_Vec3D<double>();

	rotationChanged();
}

__device__ void VX3_External::rotationChanged()
{
	if (extRotation != VX3_Vec3D<double>()){
		_extRotationQ = VX3_Quat3D<double>(extRotation);
	}
	else { //rotation is zero in all axes
		_extRotationQ = VX3_Quat3D<double>();
	}
}





#include "VX3_MaterialVoxel.h"
#include "VX3_VoxelyzeKernel.cuh"

VX3_MaterialVoxel::VX3_MaterialVoxel( CVX_MaterialVoxel *p, VX3_VoxelyzeKernel* k ):
VX3_Material( (CVX_Material*) p, k ),
nomSize(p->nomSize), gravMult(p->gravMult),_mass(p->_mass),
_massInverse(p->_massInverse), _sqrtMass(p->_sqrtMass), _firstMoment(p->_firstMoment),
_momentInertia(p->_momentInertia), _momentInertiaInverse(p->_momentInertiaInverse),
_2xSqMxExS(p->_2xSqMxExS), _2xSqIxExSxSxS(p->_2xSqIxExSxSxS) {

}

__device__ VX3_MaterialVoxel::VX3_MaterialVoxel(float youngsModulus, float density, double nominalSize) : VX3_Material(youngsModulus, density)
{
	initialize(nominalSize);
}

__device__ VX3_MaterialVoxel::VX3_MaterialVoxel(const VX3_Material& mat, double nominalSize) : VX3_Material(mat)
{
	initialize(nominalSize);
}

__device__ void VX3_MaterialVoxel::initialize(double nominalSize)
{
	nomSize = nominalSize;
	gravMult = 0.0f;
	updateDerived();
}

__device__ VX3_MaterialVoxel& VX3_MaterialVoxel::operator=(const VX3_MaterialVoxel& vIn)
{
	VX3_Material::operator=(vIn); //set base VX3_Material class variables equal

	nomSize=vIn.nomSize;
	gravMult=vIn.gravMult;
	_eHat = vIn._eHat;
	_mass=vIn._mass;
	_massInverse=vIn._massInverse;
	_sqrtMass=-vIn._sqrtMass;
	_firstMoment=vIn._firstMoment;
	_momentInertia=vIn._momentInertia;
	_momentInertiaInverse=vIn._momentInertiaInverse;
	_2xSqMxExS=vIn._2xSqMxExS;
	_2xSqIxExSxSxS=vIn._2xSqIxExSxSxS;

	return *this;
}

__device__ bool VX3_MaterialVoxel::updateDerived() 
{
	VX3_Material::updateDerived(); //update base VX3_Material class derived variables

	double volume = nomSize*nomSize*nomSize;
	_mass = (float)(volume*rho); 
	_momentInertia = (float)(_mass*nomSize*nomSize / 6.0f); //simple 1D approx
	_firstMoment = (float)(_mass*nomSize / 2.0f);

	if (volume==0 || _mass==0 || _momentInertia==0){
		_massInverse = _sqrtMass = _momentInertiaInverse = _2xSqMxExS = _2xSqIxExSxSxS = 0.0f; //zero everything out
		return false;
	}


	_massInverse = 1.0f / _mass;
	_sqrtMass = sqrt(_mass);
	_momentInertiaInverse = 1.0f / _momentInertia;
	_2xSqMxExS = (float)(2.0f*sqrt(_mass*E*nomSize));
	_2xSqIxExSxSxS = (float)(2.0f*sqrt(_momentInertia*E*nomSize*nomSize*nomSize));

	return true;
}


__device__ bool VX3_MaterialVoxel::setNominalSize(double size)
{
	if (size <= 0) size = FLT_MIN;
	nomSize=size;
	return updateDerived(); //update derived quantities
}
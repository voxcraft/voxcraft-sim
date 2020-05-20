#if !defined(VX3_MATERIAL_H)
#define VX3_MATERIAL_H

#include <string>
#include <vector>
#include "VX3_vector.cuh"

#include "VX3.cuh"
#include "VX_Material.h" // for CVX_Material

class VX3_VoxelyzeKernel;

class VX3_Material {
public:

	VX3_Material( CVX_Material* p, VX3_VoxelyzeKernel* k );

	__device__ VX3_Material(float youngsModulus=1e6f, float density=1e3f); //!< Default Constructor. @param[in] youngsModulus The Young's Modulus (stiffness) of this material in Pascals. @param[in] density The density of this material in Kg/m^3
	//__device__ virtual ~VX3_Material(void) {}; //!< Destructor. Specified as virtual so we can just keep track of generic material pointers for voxel and link materials.
	__device__ VX3_Material(const VX3_Material& vIn) {*this = vIn;} //!< Copy constructor
	__device__ VX3_Material& operator=(const VX3_Material& vIn); //!< Equals operator

	__device__ void clear(); //!< Resets all material information to default.
	//__device__ const char* lastError() const {return error.c_str();} //!< Returns the last error encountered for this object.

	//__device__ void setName(const char* name) {myName = std::string(name);} //!< Adds an optional name to the material. @param[in] name Desired name. 
	//__device__ const char* name() const {return myName.c_str();} //!< Returns the optional material name if one was specifed. Otherwise returns an empty string.

	__device__ float stress(float strain, float transverseStrainSum=0.0f, bool forceLinear = false); //!<returns the stress of the material model accounting for volumetric strain effects. @param [in] strain The strain to query. The resulting stress in this direction will be returned. @param [in] transverseStrainSum The sum of the two principle normal strains in the plane perpendicular to strain. @param [in] forceLinear If true, the result will be calculated according to the elastic modulus of the material regardless of non-linearities in the model.
	__device__ float modulus(float strain); //!<returns the modulus (slope of the stress/strain curve) of the material model at the specified strain. @param [in] strain The strain to query.
	__device__ bool isYielded(float strain) {return epsilonYield != -1.0f && strain>epsilonYield;} //!< Returns true if the specified strain is past the yield point (if one is specified). @param [in] strain The strain to query.
	__device__ bool isFailed(float strain) {return epsilonFail != -1.0f && strain>epsilonFail;} //!< Returns true if the specified strain is past the failure point (if one is specified). @param [in] strain The strain to query.

	//color
	__device__ void setColor(int red, int green, int blue, int alpha=255); //!< Sets the material color. Values from [0,255]. @param [in] red Red channel @param [in] green Green channel @param [in] blue Blue channel @param [in] alpha Alpha channel
	__device__ void setRed(int red); //!< Sets the red channel of the material color. @param [in] red Red channel [0,255]
	__device__ void setGreen(int green); //!< Sets the green channel of the material color. @param [in] green Green channel [0,255]
	__device__ void setBlue(int blue); //!< Sets the blue channel of the material color. @param [in] blue Blue channel [0,255]
	__device__ void setAlpha(int alpha); //!< Sets the alpha channel of the material color. @param [in] alpha Alpha channel [0,255]
	__device__ int red() const {return r;} //!< Returns the red channel of the material color [0,255] or -1 if unspecified.
	__device__ int green() const {return g;} //!< Returns the green channel of the material color [0,255] or -1 if unspecified.
	__device__ int blue() const {return b;} //!< Returns the blue channel of the material color [0,255] or -1 if unspecified.
	__device__ int alpha() const {return a;} //!< Returns the alpha channel of the material color [0,255] or -1 if unspecified.

	//Material Model
	__device__ bool setModel(int dataPointCount, float* pStrainValues, float* pStressValues); //!< Defines the physical material behavior with a series of true stress/strain data points. @param [in] dataPointCount The expected number of data points. @param [in] pStrainValues pointer to the first strain value data point in a contiguous array (Units: Pa). @param [in] pStressValues pointer to the first stress value data point in a contiguous array (Units: Pa).
	__device__ bool setModelLinear(float youngsModulus, float failureStress=-1); //!< Convenience function to quickly define a linear material. @param [in] youngsModulus Young's modulus (Units: Pa). @param [in] failureStress Optional failure stress (Units: Pa). -1 indicates failure is not an option.
	__device__ bool setModelBilinear(float youngsModulus, float plasticModulus, float yieldStress, float failureStress=-1); //!< Convenience function to quickly define a bilinear material. @param [in] youngsModulus Young's modulus (Units: Pa). @param [in] plasticModulus Plastic modulus (Units: Pa). @param [in] yieldStress Yield stress. @param [in] failureStress Optional failure stress (Units: Pa). -1 indicates failure is not an option.
	__device__ bool isModelLinear() const {return linear;} //!< Returns true if the material model is a simple linear behavior.
	
	__device__ float youngsModulus() const {return E;} //!< Returns Youngs modulus in Pa.
	__device__ float yieldStress() const {return sigmaYield;} //!<Returns the yield stress in Pa or -1 if unspecified.
	__device__ float failureStress() const {return sigmaFail;} //!<Returns the failure stress in Pa or -1 if unspecified.
	__device__ int modelDataPoints() {return d_strainData.size();} //!< Returns the number of data points in the current material model data arrays.
	__device__ const float* modelDataStrain() {return &(d_strainData[0]);} //!< Returns a pointer to the first strain value data point in a continuous array. The number of values can be determined from modelDataPoints(). The assumed first value of 0 is included.
	__device__ const float* modelDataStress() {return &(d_stressData[0]);} //!< Returns a pointer to the first stress value data point in a continuous array. The number of values can be determined from modelDataPoints(). The assumed first value of 0 is included.

	__device__ void setPoissonsRatio(float poissonsRatio); //!< Defines Poisson's ratio for the material. @param [in] poissonsRatio Desired Poisson's ratio [0, 0.5).
	__device__ float poissonsRatio() const {return nu;} //!< Returns the current Poissons ratio
	__device__ float bulkModulus() const {return E/(3*(1-2*nu));} //!< Calculates the bulk modulus from Young's modulus and Poisson's ratio.
	__device__ float lamesFirstParameter() const {return (E*nu)/((1+nu)*(1-2*nu));} //!< Calculates Lame's first parameter from Young's modulus and Poisson's ratio.
	__device__ float shearModulus() const {return E/(2*(1+nu));} //!< Calculates the shear modulus from Young's modulus and Poisson's ratio.
	__device__ bool isXyzIndependent() const {return nu==0.0f;} //!< Returns true if poisson's ratio is zero - i.e. deformations in each dimension are independent of those in other dimensions.

	//other material properties
	__device__ void setDensity(float density); //!< Defines the density for the material in Kg/m^3. @param [in] density Desired density (0, INF)
	__device__ float density() const {return rho;} //!< Returns the current density.
	__device__ void setStaticFriction(float staticFrictionCoefficient); //!< Defines the coefficient of static friction. @param [in] staticFrictionCoefficient Coefficient of static friction [0, INF).
	__device__ float staticFriction() const {return muStatic;} //!< Returns the current coefficient of static friction.
	__device__ void setKineticFriction(float kineticFrictionCoefficient); //!< Defines the coefficient of kinetic friction. @param [in] kineticFrictionCoefficient Coefficient of kinetc friction [0, INF).
	__device__ float kineticFriction() const {return muKinetic;} //!< Returns the current coefficient of kinetic friction.

	//damping
	//http://www.roush.com/Portals/1/Downloads/Articles/Insight.pdf
	__device__ void setInternalDamping(float zeta); //!<Defines the internal material damping ratio. The effect is to damp out vibrations within a structure. zeta = mu/2 (mu = loss factor) = 1/(2Q) (Q = amplification factor). High values of zeta may lead to simulation instability. Recommended value: 1.0.  @param [in] zeta damping ratio [0, INF). (unitless)
	__device__ float internalDamping() const {return zetaInternal;} //!< Returns the internal material damping ratio.
	__device__ void setGlobalDamping(float zeta); //!<Defines the viscous damping of any voxels using this material relative to ground (no motion). Translation C (damping coefficient) is calculated according to zeta*2*sqrt(m*k) where k=E*nomSize. Rotational damping coefficient is similarly calculated High values relative to 1.0 may cause simulation instability. @param [in] zeta damping ratio [0, INF). (unitless)
	__device__ float globalDamping() const {return zetaGlobal;} //!< Returns the global material damping ratio.
	__device__ void setCollisionDamping(float zeta); //!<Defines the material damping ratio for when this material collides with something. This gives some control over the elasticity of a collision. A value of zero results in a completely elastic collision. @param [in] zeta damping ratio [0, INF). (unitless)
	__device__ float collisionDamping() const {return zetaCollision;} //!< Returns the collision material damping ratio.


	//size and scaling
	__device__ void setExternalScaleFactor(VX3_Vec3D<double> factor); //!< Scales all voxels of this material by a specified factor in each dimension (1.0 is no scaling). This allows enables volumetric displacement-based actuation within a structure. As such, mass is unchanged when the external scale factor changes. Actual size is obtained by multiplying nominal size by the provided factor. @param[in] factor Multiplication factor (0, INF) for the size of all voxels of this material in its local x, y, and z axes. (unitless) 
	__device__ void setExternalScaleFactor(double factor) {setExternalScaleFactor(VX3_Vec3D<double>(factor, factor, factor));} //!< Convenience function to specify isotropic external scaling factor. See setExternalScaleFactor(VX3_Vec3D<> factor). @param[in] factor external scaling factor (0, INF).
	__device__ VX3_Vec3D<double> externalScaleFactor() {return extScale;} //!< Returns the current external scaling factor (unitless). See description of setExternalScaleFactor().

	//thermal expansion
	__device__ void setCte(float cte) {alphaCTE=cte;} //!< Defines the coefficient of thermal expansion. @param [in] cte Desired coefficient of thermal expansion per degree C (-INF, INF)
	__device__ float cte() const {return alphaCTE;} //!< Returns the current coefficient of thermal expansion per degree C.

	//derived quantities to cache
	__device__ virtual bool updateAll() {return false;} //!< Updates and recalculates eveything possible (used by inherited classed when material properties have changed)
	__device__ virtual bool updateDerived(); //!< Updates all the derived quantities cached as member variables for this and derived classes. (Especially if density, size or elastic modulus changes.)

	__device__ bool setYieldFromData(float percentStrainOffset=0.2); //!< Sets sigmaYield and epsilonYield assuming strainData, stressData, E, and failStrain are set correctly.
	__device__ float strain(float stress); //!< Returns a simple reverse lookup of the first strain that yields this stress from data point lookup.

	//VX3_vector
	__device__ void syncVectors();



/* data */
	//std::string error; //!< The last error encountered
	//std::string myName; //!< The name of this material. Default is "".
	int r; //!< Red color value of this material from 0-255. Default is -1 (invalid/not set).
	int g; //!< Green color value of this material from 0-255. Default is -1 (invalid/not set).
	int b; //!< Blue color value of this material from 0-255. Default is -1 (invalid/not set).
	int a; //!< Alpha value of this material from 0-255. Default is -1 (invalid/not set).

	int matid=0;

	//material model
	bool fixed = false;
	bool sticky = false;
	double Cilia = 0;
	bool linear; //!< Set to true if this material is specified as linear.
	float E; //!< Young's modulus (stiffness) in Pa.
	float sigmaYield; //!< Yield stress in Pa.
	float sigmaFail; //!< Failure stress in Pa
	float epsilonYield; //!< Yield strain
	float epsilonFail; //!< Failure strain
	VX3_hdVector<float> hd_strainData;
	VX3_dVector<float> d_strainData;
	VX3_hdVector<float> hd_stressData;
	VX3_dVector<float> d_stressData;

	float nu; //!< Poissonss Ratio
	float rho; //!< Density in Kg/m^3
	float alphaCTE; //!< Coefficient of thermal expansion (CTE)
	float muStatic; //!< Static coefficient of friction
	float muKinetic; //!< Kinetic coefficient of friction
	float zetaInternal; //!< Internal damping ratio
	float zetaGlobal; //!< Global damping ratio
	float zetaCollision; //!< Collision damping ratio

	VX3_Vec3D<double> extScale; //!< A prescribed scaling factor. default of (1,1,1) is no scaling.

	float _eHat; //!< Cached effective elastic modulus for materials with non-zero Poisson's ratio.

	VX3_dVector<VX3_Material*> d_dependentMaterials; //!< Any materials in this list will have updateDerived() called whenever it's called for this material. For example, in Voxelyze this is used for updatng link materials when one or both voxel materials change
	VX3_hdVector<VX3_Material*>hd_dependentMaterials;
	//Future parameters:
	//piezo?
	//compressive strength? (/compressive data)
	//heat conduction
	// dependentMaterials only for combinedMaterial, not used.

	bool isTarget = false;
	bool isPaceMaker = false;
	double PaceMakerPeriod=0;
	bool isElectricalActive=false;

	double RemoveFromSimulationAfterThisManySeconds = 0.0;
	double TurnOnThermalExpansionAfterThisManySeconds = 0.0;
	double TurnOnCiliaAfterThisManySeconds = 0.0;

	double signalValueDecay = 0.9; // ratio from [0,1]
	double signalTimeDelay = 0.0; // in sec
	double inactivePeriod = 0.05; //in sec
	int isMeasured=1;


	// for Secondary Experiment
    bool removed = false;
};




#endif // VX3_MATERIAL_H

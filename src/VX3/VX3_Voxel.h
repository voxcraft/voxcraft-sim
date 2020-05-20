#if !defined(VX3_VOXEL_H)
#define VX3_VOXEL_H
#include "VX3.cuh"

#include "VX_Voxel.h" // for CVX_Voxel

#include "VX3_External.h"
#include "VX3_Link.h"
#include "VX3_MaterialVoxel.h"
#include "VX3_Signal.h"
#include "VX3_queue.cuh"

class VX3_Collision;
class VX3_VoxelyzeKernel;

class VX3_Voxel {
  public:
    VX3_Voxel() = default;
    VX3_Voxel(CVX_Voxel *p, VX3_VoxelyzeKernel *k);
    ~VX3_Voxel();

    __device__ VX3_Link *link(linkDirection direction) const {
        return links[direction];
    } //!< Returns a pointer to the link object in the specified direction if it exists. Returns null if a link does not exist in this
      //!< direction.
    __device__ int linkCount() const {
        int retVal = 0;
        for (int i = 0; i < 6; i++)
            if (links[i])
                retVal++;
        return retVal;
    } //!< Returns the number of links present for this voxel out of a total 6 possible.
    __device__ VX3_Voxel *adjacentVoxel(linkDirection direction)
        const; //!< Returns a pointer to the voxel in the specified direction if one exists, or NULL otherwise. @param[in] direction
               //!< Positive or negative X, Y, or Z direction according to the linkDirection enum.
    __device__ short indexX() { return ix; } //!< Returns the global X index of this voxel.
    __device__ short indexY() { return iy; } //!< Returns the global Y index of this voxel.
    __device__ short indexZ() { return iz; } //!< Returns the global Z index of this voxel.

    __device__ VX3_MaterialVoxel *material() {
        return mat;
    } //!< Returns the linked material object containing the physical properties of this voxel.

    __device__ bool externalExists() {
        return ext ? true : false;
    } //!< Returns true if this voxel has had its VX3_External object created. This does not mecessarily imply that this external object
      //!< actually contains any fixes or forces.
    __device__ VX3_External *external() {
        if (!ext)
            ext = new VX3_External();
        return ext;
    } //!< Returns a pointer to this voxel's unique external object that contains fixes, forces, and/or displacements. Allocates a new empty
      //!< one if it doesn't already exist. Use externalExists() to determine if external() has been previously called at any time.

    __device__ void
    timeStep(double dt, double currentTime,
             VX3_VoxelyzeKernel *k); //!< Advances this voxel's state according to all forces and moments acting on it. Large timesteps will
                                     //!< cause instability. Use CVoxelyze::recommendedTimeStep() to get the recommended largest stable
                                     //!< timestep. @param[in] dt Timestep (in second) to advance.

    // physical location
    __device__ VX3_Vec3D<double> position() const {
        return pos;
    } //!< Returns the center position of this voxel in meters (GCS). This is the origin of the local coordinate system (LCS).
    __device__ VX3_Vec3D<double> originalPosition() const {
        double s = mat->nominalSize();
        return VX3_Vec3D<double>(ix * s, iy * s, iz * s);
    } //!< Returns the initial (nominal) position of this voxel.
    __device__ VX3_Vec3D<double> displacement() const {
        return (pos - originalPosition());
    } //!< Returns the 3D displacement of this voxel from its original location in meters (GCS)/
    __device__ VX3_Vec3D<float> size() const {
        return cornerOffset(PPP) - cornerOffset(NNN);
    } //!< Returns the current deformed size of this voxel in the local voxel coordinates system (LCS). If asymmetric forces are acting on
      //!< this voxel, the voxel may not be centered on position(). Use cornerNegative() and cornerPositive() to determine this information.
    __device__ VX3_Vec3D<float> cornerPosition(
        voxelCorner corner) const; //!< Returns the deformed location of the voxel corner in the specified corner in the global coordinate
                                   //!< system (GCS). Essentially cornerOffset() with the voxel's current global position/rotation applied.
    __device__ VX3_Vec3D<float> cornerOffset(voxelCorner corner)
        const; //!< Returns the deformed location of the voxel corner in the specified corner in the local voxel coordinate system (LCS).
               //!< Used to draw the deformed voxel in the correct position relative to the position().
    __device__ bool isInterior() const {
        return (boolStates & SURFACE) ? true : false;
    } //!< Returns true if the voxel is surrounded by other voxels on its 6 coordinate faces. Returns false if 1 or more faces are exposed.
    __device__ bool isSurface() const {
        return !isInterior();
    } //!< Convenience function to enhance code readibility. The inverse of isInterior(). Returns true 1 or more faces are exposed. Returns
      //!< false if the voxel is surrounded by other voxels on its 6 coordinate faces.

    __device__ VX3_Vec3D<double> baseSize() const {
        return mat->size() * (1 + tempe * mat->alphaCTE);
    } //!< Returns the nominal size of this voxel (LCS) accounting for any specified temperature and external actuation. Specifically,
      //!< returns the zero-stress size of the voxel if all forces/moments were removed.

    // Change to: Actuation due to Voltage!!
    __device__ double baseSize(linkAxis axis) const {
        return mat->size()[axis] * (1 + tempe * mat->alphaCTE);
        // return mat->size()[axis] * (1 + voltage * mat->alphaCTE);
    } //!< Returns the nominal size of this voxel in the specified axis accounting for any specified temperature and external actuation.
      //!< Specifically, returns the zero-stress dimension of the voxel if all forces/moments were removed.

    __device__ double baseSizeAverage() const {
        VX3_Vec3D<double> bSize = baseSize();
        return (bSize.x + bSize.y + bSize.z) / 3.0f;
    } //!< Returns the average nominal size of the voxel in a zero-stress (no force) state. (X+Y+Z/3)

    __device__ VX3_Quat3D<double> orientation() const {
        return orient;
    } //!< Returns the orientation of this voxel in quaternion form (GCS). This orientation defines the relative orientation of the local
      //!< coordinate system (LCS). The unit quaternion represents the original orientation of this voxel.
    __device__ float orientationAngle() const {
        return (float)orient.Angle();
    } //!< Use with orientationAxis() to get the orientation of this voxel in angle/axis form. Returns the angle in radians.
    __device__ VX3_Vec3D<double> orientationAxis() const {
        return orient.Axis();
    } //!< Use with orientationAngle() to get the orientation of this voxel in angle/axis form. Returns a unit vector in the global
      //!< coordinate system (GCS).

    __device__ float displacementMagnitude() const {
        return (float)displacement().Length();
    } //!< Returns the distance (magnitude of displacement) this voxel has moved from its initial nominal position. (GCS)
    __device__ float angularDisplacementMagnitude() const {
        return (float)orient.Angle();
    } //!< Returns the angle (magnitude of angular displacement) this voxel has rotated from its initial nominal orientation. (GCS)
    __device__ VX3_Vec3D<double> velocity() const {
        return linMom * mat->_massInverse;
    } //!< Returns the 3D velocity of this voxel in m/s (GCS)
    __device__ float velocityMagnitude() const {
        return (float)(linMom.Length() * mat->_massInverse);
    } //!< Returns the velocity of this voxel in m/s.
    __device__ VX3_Vec3D<double> angularVelocity() const {
        return angMom * mat->_momentInertiaInverse;
    } //!< Returns the 3D angular velocity of this voxel in rad/s (GCS)
    __device__ float angularVelocityMagnitude() const {
        return (float)(angMom.Length() * mat->_momentInertiaInverse);
    } //!< Returns the angular velocity of this voxel in rad/s.
    __device__ float kineticEnergy() const {
        return (float)(0.5 * (mat->_massInverse * linMom.Length2() + mat->_momentInertiaInverse * angMom.Length2()));
    } //!< Returms the kinetic energy of this voxel in Joules.
    __device__ float volumetricStrain() const {
        return (float)(strain(false).x + strain(false).y + strain(false).z);
    } //!< Returns the volumetric strain of the voxel according to the definition at
      //!< http://www.colorado.edu/engineering/CAS/courses.d/Structures.d/IAST.Lect05.d/IAST.Lect05.pdf
    __device__ float pressure() const {
        return -mat->youngsModulus() * volumetricStrain() / (3 * (1 - 2 * mat->poissonsRatio()));
    } //!< Returns the engineering internal "pressure" in Pa according to the definition at
      //!< http://www.colorado.edu/engineering/CAS/courses.d/Structures.d/IAST.Lect05.d/IAST.Lect05.pdf

    // material state
    __device__ bool
    isYielded() const; //!< Returns true if the stress in this voxel has ever exceeded the yield stress. Technically, this returns true if
                       //!< any of the connected links have yielded since the stress state of the voxel is never expressly calculated.
    __device__ bool
    isFailed() const; //!< Returns true if the stress in this voxel has ever exceeded the failure stress. Technically, this returns true if
                      //!< any of the connected links have failed since the stress state of the voxel is never expressly calculated.

    //@ voxel level for heat diffusion experiments later
    __device__ float temperature() { return tempe; } //!< Returns the current temperature of this voxel in degrees Celsius.
    __device__ void
    setTemperature(float temperature); //!< Specifies the temperature for this voxel. This adds (or subtracts) the correct amount of thermal
                                       //!< energy to leave the voxel at ths specified temperature, but this temperature will not be
                                       //!< maintaned without subsequent determines the amount of scaling from the temperature

    __device__ VX3_Vec3D<float> externalForce();  //!< Returns the current external force applied to this voxel in newtons. If the voxel is
                                                  //!< not fixed this will return any applied external forces. If fixed it will return the
                                                  //!< current reaction force necessary to enforce the zero-motion constraint.
    __device__ VX3_Vec3D<float> externalMoment(); //!< Returns the current external moment applied to this voxel in N-m. If the voxel is not
                                                  //!< fixed this will return any applied external moments. If fixed it will return the
                                                  //!< current reaction moment necessary to enforce the zero-motion constraint.

    __device__ void haltMotion() {
        linMom = angMom = VX3_Vec3D<>(0, 0, 0);
    } //!< Halts all momentum of this block. Unless fixed the voxel will continue to move in subsequent timesteps.

    __device__ void enableFloor(bool enabled) {
        enabled ? boolStates |= FLOOR_ENABLED : boolStates &= ~FLOOR_ENABLED;
    } //!< Enables this voxel interacting with the floor at Z=0. @param[in] enabled Enable interaction
    __host__ __device__ bool isFloorEnabled() const {
        return boolStates & FLOOR_ENABLED ? true : false;
    } //!< Returns true of this voxel will interact with the floor at Z=0.
    __device__ bool isFloorStaticFriction() const {
        return boolStates & FLOOR_STATIC_FRICTION ? true : false;
    } //!< Returns true if this voxel is in contact with the floor and stationary in the horizontal directions. This corresponds to that
      //!< voxel being in the mode of static friction (as opposed to kinetic) with the floor.
    __device__ float floorPenetration() const {
        return (float)(baseSizeAverage() / 2 - mat->nominalSize() / 2 - pos.z);
    } //!< Returns the interference (in meters) between the collision envelope of this voxel and the floor at Z=0. Positive numbers
      //!< correspond to interference. If the voxel is not touching the floor 0 is returned.

    __device__ VX3_Vec3D<double>
    force(); //!< Calculates and returns the sum of the current forces on this voxel. This would normally only be called internally, but can
             //!< be used to query the state of a voxel for visualization or debugging.
    __device__ VX3_Vec3D<double>
    moment(); //!< Calculates and returns the sum of the current moments on this voxel. This would normally only be called internally, but
              //!< can be used to query the state of a voxel for visualization or debugging.

    __device__ float
    transverseArea(linkAxis axis); //!< Returns the transverse area of this voxel with respect to the specified axis. This would normally be
                                   //!< called only internally, but can be used to calculate the correct relationship between force and
                                   //!< stress for this voxel if Poisson's ratio is non-zero.
    __device__ float
    transverseStrainSum(linkAxis axis); //!< Returns the sum of the current strain of this voxel in the two mutually perpindicular axes to
                                        //!< the specified axis. This would normally be called only internally, but can be used to correctly
                                        //!< calculate stress for this voxel if Poisson's ratio is non-zero.

    __device__ float dampingMultiplier() {
        return 2 * mat->_sqrtMass * mat->zetaInternal / previousDt;
    } //!< Returns the damping multiplier for this voxel. This would normally be called only internally for the internal damping
      //!< calculations.

    // a couple global convenience functions to have wherever the link enums are used
    __device__ static inline linkAxis toAxis(linkDirection direction) {
        return (linkAxis)((int)direction / 2);
    } //!< Returns the link axis of the specified link direction.
    __device__ static inline linkDirection toDirection(linkAxis axis, bool positiveDirection) {
        return (linkDirection)(2 * ((int)axis) + positiveDirection ? 0 : 1);
    } //!< Returns the link direction of the specified link axis and sign.
    __device__ static inline bool isNegative(linkDirection direction) {
        return direction % 2 == 1;
    } //!< Returns true if the specified link direction is negative.
    __device__ static inline bool isPositive(linkDirection direction) {
        return direction % 2 == 0;
    } //!< Returns true if the specified link direction is positive.
    __device__ static inline linkDirection toOpposite(linkDirection direction) {
        return (linkDirection)(direction - direction % 2 + (direction + 1) % 2);
    } //!< Returns the opposite (negated) link direction of the specified direction.

    __device__ void replaceMaterial(
        VX3_MaterialVoxel *newMaterial); //!< Replaces the material properties of this voxel (but not links) to this new CVX_Material. May
                                         //!< cause unexpected behavior if certain material properties are changed mid-simulation. @param
                                         //!< [in] newMaterial The new material properties for this voxel.

    __device__ void addLinkInfo(linkDirection direction,
                                VX3_Link *link); // adds the information about a link connected to this voxel in the specified direction
    __device__ void
    removeLinkInfo(linkDirection direction); // removes the information about a link connected to this voxel in the specified direction

    __device__ void setFloorStaticFriction(bool active) {
        active ? boolStates |= FLOOR_STATIC_FRICTION : boolStates &= ~FLOOR_STATIC_FRICTION;
    }

    __device__ void
    floorForce(float dt,
               VX3_Vec3D<double> *pTotalForce); // modifies pTotalForce to include the object's interaction with a floor. This should be
                                                // calculated as the last step of sumForce so that pTotalForce is complete.

    __device__ VX3_Vec3D<float>
    strain(bool poissonsStrain) const; // LCS returns voxel strain. if tensionStrain true and no actual tension in that
    __device__ VX3_Vec3D<float> poissonsStrain();

    __device__ void eulerStep(float dt); // execute an euler time step at the specified dt

    __device__ void updateSurface();
    __device__ void enableCollisions(bool enabled, float watchRadius = 0.0f); // watchRadius in voxel units
    __device__ bool isCollisionsEnabled() const { return boolStates & COLLISIONS_ENABLED ? true : false; }
    __device__ void generateNearby(int linkDepth, int gindex, bool surfaceOnly = true);

    __device__ void propagateSignal(double currentTime);
    __device__ void packMaker(double currentTime);
    __device__ void localSignalDecay(double currentTime);

    __device__ void receiveSignal(double value, double currentTime, bool force);
    __device__ void getSignal(double currentTime);

    __device__ void syncVectors();

    /* data */
    CVX_Voxel *_voxel;

    VX3_MaterialVoxel *mat = NULL;
    int matid;
    short ix, iy, iz;
    VX3_External *ext = NULL;

    VX3_Link *links[6]; // links in the 6 cardinal directions according to linkDirection enumeration

    // voxel state
    VX3_Vec3D<double> pos;     // current center position (meters) (GCS)
    VX3_Vec3D<double> linMom;  // current linear momentum (kg*m/s) (GCS)
    VX3_Quat3D<double> orient; // current orientation (GCS)
    VX3_Vec3D<double> angMom;  // current angular momentum (kg*m^2/s) (GCS)

    voxState boolStates; // single int to store many boolean state values as bit flags according to

    float tempe; // 0 is no expansion

    VX3_Vec3D<float> pStrain;   // cached poissons strain
    bool poissonsStrainInvalid; // flag for recomputing poissons strain.

    float previousDt; // remember the duration of the last timestep of this voxel

    VX3_Vec3D<float> lastColWatchPosition;

    double phaseOffset;
    bool isDetached; // true if the voxel is on main body, false if it fell on the ground.

    VX3_Vec3D<> contactForce;
    VX3_Vec3D<> baseCiliaForce;
    VX3_Vec3D<> shiftCiliaForce;
    VX3_Vec3D<> CiliaForce;

    bool enableAttach = true;

    VX3_dQueue<VX3_Signal *> d_signals;
    VX3_Signal d_signal;
    double localSignal = 0.0;
    double localSignaldt = 0.0;
    double packmakerNextPulse = 0.0;
    double inactiveUntil = 0.0;


    // for Secondary Experiment
    bool removed = false;
    bool isMeasured = true;
};

#endif // VX3_VOXEL_H

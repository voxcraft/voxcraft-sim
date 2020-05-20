#include "VX3_Voxel.h"
#include <vector>

#include "VX3_External.h"
#include "VX3_Link.h"
#include "VX3_MaterialVoxel.h"
#include "VX3_MemoryCleaner.h"
#include "VX3_VoxelyzeKernel.cuh"

VX3_Voxel::VX3_Voxel(CVX_Voxel *p, VX3_VoxelyzeKernel *k)
    : matid(p->matid), ix(p->ix), iy(p->iy), iz(p->iz), pos(p->pos), linMom(p->linMom), orient(p->orient), angMom(p->angMom),
      boolStates(p->boolStates), tempe(p->temp), pStrain(p->pStrain), poissonsStrainInvalid(p->poissonsStrainInvalid),
      previousDt(p->previousDt), phaseOffset(p->phaseOffset), isDetached(p->isDetached), baseCiliaForce(p->baseCiliaForce),
      shiftCiliaForce(p->shiftCiliaForce) {
    _voxel = p;

    for (int i = 0; i < k->num_d_voxelMats; i++) {
        if (k->h_voxelMats[i] == p->mat) {
            mat = &k->d_voxelMats[i];
            break;
        }
    }

    for (unsigned i = 0; i < 6; i++) {
        if (p->links[i]) {
            links[i] = k->h_lookup_links[p->links[i]];
            // for (int j=0;j<k->num_d_links;j++) {
            // 	if (p->links[i] == k->h_links[j]) {
            // 		links[i] = &k->d_links[j];
            // 		continue;
            // 	}
            // }
        } else {
            links[i] = NULL;
        }
    }

    // mat = new VX3_MaterialVoxel(p->mat);
    if (p->ext) {
        VcudaMalloc((void **)&ext, sizeof(VX3_External));
        VX3_External temp2(p->ext);
        VcudaMemcpy(ext, &temp2, sizeof(VX3_External), VcudaMemcpyHostToDevice);
    } else {
        ext = NULL;
    }
}

VX3_Voxel::~VX3_Voxel() {
    if (ext) {
        MycudaFree(ext);
        ext = NULL;
    }
}

__device__ void VX3_Voxel::syncVectors() {
    d_signal.value = 0;
    d_signal.activeTime = 0;

    d_signals.clear();
}
__device__ VX3_Voxel *VX3_Voxel::adjacentVoxel(linkDirection direction) const {
    VX3_Link *pL = links[(int)direction];
    if (pL)
        return pL->voxel(true) == this ? pL->voxel(false) : pL->voxel(true);
    else
        return NULL;
}

__device__ void VX3_Voxel::addLinkInfo(linkDirection direction, VX3_Link *link) {
    links[direction] = link;
    updateSurface();
}

__device__ void VX3_Voxel::removeLinkInfo(linkDirection direction) {
    links[direction] = NULL;
    updateSurface();
}

__device__ void VX3_Voxel::replaceMaterial(VX3_MaterialVoxel *newMaterial) {
    if (newMaterial != NULL) {

        linMom *= newMaterial->_mass / mat->_mass; // adjust momentums to keep velocity constant across material change
        angMom *= newMaterial->_momentInertia / mat->_momentInertia;
        setFloorStaticFriction(false);
        poissonsStrainInvalid = true;

        mat = newMaterial;
    }
}

__device__ bool VX3_Voxel::isYielded() const {
    for (int i = 0; i < 6; i++) {
        if (links[i] && links[i]->isYielded())
            return true;
    }
    return false;
}

__device__ bool VX3_Voxel::isFailed() const {
    for (int i = 0; i < 6; i++) {
        if (links[i] && links[i]->isFailed())
            return true;
    }
    return false;
}

__device__ void VX3_Voxel::setTemperature(float temperature) {
    tempe = temperature;
    for (int i = 0; i < 6; i++) {
        if (links[i] != NULL)
            links[i]->updateRestLength();
    }
}

__device__ VX3_Vec3D<float> VX3_Voxel::externalForce() {
    VX3_Vec3D<float> returnForce(ext->force());
    if (ext->isFixed(X_TRANSLATE) || ext->isFixed(Y_TRANSLATE) || ext->isFixed(Z_TRANSLATE)) {
        VX3_Vec3D<float> thisForce = (VX3_Vec3D<float>)-force();
        if (ext->isFixed(X_TRANSLATE))
            returnForce.x = thisForce.x;
        if (ext->isFixed(Y_TRANSLATE))
            returnForce.y = thisForce.y;
        if (ext->isFixed(Z_TRANSLATE))
            returnForce.z = thisForce.z;
    }
    return returnForce;
}

__device__ VX3_Vec3D<float> VX3_Voxel::externalMoment() {
    VX3_Vec3D<float> returnMoment(ext->moment());
    if (ext->isFixed(X_ROTATE) || ext->isFixed(Y_ROTATE) || ext->isFixed(Z_ROTATE)) {
        VX3_Vec3D<float> thisMoment = (VX3_Vec3D<float>)-moment();
        if (ext->isFixed(X_ROTATE))
            returnMoment.x = thisMoment.x;
        if (ext->isFixed(Y_ROTATE))
            returnMoment.y = thisMoment.y;
        if (ext->isFixed(Z_ROTATE))
            returnMoment.z = thisMoment.z;
    }
    return returnMoment;
}

__device__ VX3_Vec3D<float> VX3_Voxel::cornerPosition(voxelCorner corner) const {
    return (VX3_Vec3D<float>)pos + orient.RotateVec3D(cornerOffset(corner));
}

__device__ VX3_Vec3D<float> VX3_Voxel::cornerOffset(voxelCorner corner) const {
    VX3_Vec3D<> strains;
    for (int i = 0; i < 3; i++) {
        bool posLink = corner & (1 << (2 - i)) ? true : false;
        VX3_Link *pL = links[2 * i + (posLink ? 0 : 1)];
        if (pL && !pL->isFailed()) {
            strains[i] = (1 + pL->axialStrain(posLink)) * (posLink ? 1 : -1);
        } else
            strains[i] = posLink ? 1.0 : -1.0;
    }

    return (0.5 * baseSize()).Scale(strains);
}

// http://klas-physics.googlecode.com/svn/trunk/src/general/Integrator.cpp (reference)
__device__ void VX3_Voxel::timeStep(double dt, double currentTime, VX3_VoxelyzeKernel *k) {
    previousDt = dt;
    if (dt == 0.0f)
        return;

    if (ext && ext->isFixedAll()) {

        pos = originalPosition() + ext->translation();
        orient = ext->rotationQuat();
        haltMotion();
        return;
    }

    // Translation
    VX3_Vec3D<double> curForce = force();

    // Apply Force Field
    curForce.x += k->force_field.x_forcefield(pos.x, pos.y, pos.z, k->collisionCount, currentTime, k->recentAngle, k->targetCloseness,
                                              k->numClosePairs, k->num_d_voxels);
    curForce.y += k->force_field.y_forcefield(pos.x, pos.y, pos.z, k->collisionCount, currentTime, k->recentAngle, k->targetCloseness,
                                              k->numClosePairs, k->num_d_voxels);
    curForce.z += k->force_field.z_forcefield(pos.x, pos.y, pos.z, k->collisionCount, currentTime, k->recentAngle, k->targetCloseness,
                                              k->numClosePairs, k->num_d_voxels);

    VX3_Vec3D<double> fricForce = curForce;

    if (isFloorEnabled()) {
        floorForce(dt, &curForce); // floor force needs dt to calculate threshold to "stop" a slow voxel into static friction.
    }
    fricForce = curForce - fricForce;

    assert(!(curForce.x != curForce.x) || !(curForce.y != curForce.y) || !(curForce.z != curForce.z)); // assert non QNAN
    linMom += curForce * dt;

    VX3_Vec3D<double> translate(linMom * (dt * mat->_massInverse)); // movement of the voxel this timestep

    //	we need to check for friction conditions here (after calculating the translation) and stop things accordingly
    if (isFloorEnabled() &&
        floorPenetration() >=
            0) { // we must catch a slowing voxel here since it all boils down to needing access to the dt of this timestep.

        double work = fricForce.x * translate.x + fricForce.y * translate.y;                // F dot disp
        double hKe = 0.5 * mat->_massInverse * (linMom.x * linMom.x + linMom.y * linMom.y); // horizontal kinetic energy
        if (hKe + work <= 0)
            setFloorStaticFriction(true); // this checks for a change of direction according to the work-energy principle

        if (isFloorStaticFriction()) { // if we're in a state of static friction, zero out all horizontal motion
            linMom.x = linMom.y = 0;
            translate.x = translate.y = 0;
        }
    } else
        setFloorStaticFriction(false);

    pos += translate;

    // Rotation
    VX3_Vec3D<> curMoment = moment();
    angMom += curMoment * dt;

    orient = VX3_Quat3D<>(angMom * (dt * mat->_momentInertiaInverse)) * orient; // update the orientation
    if (ext) {
        double size = mat->nominalSize();
        if (ext->isFixed(X_TRANSLATE)) {
            pos.x = ix * size + ext->translation().x;
            linMom.x = 0;
        }
        if (ext->isFixed(Y_TRANSLATE)) {
            pos.y = iy * size + ext->translation().y;
            linMom.y = 0;
        }
        if (ext->isFixed(Z_TRANSLATE)) {
            pos.z = iz * size + ext->translation().z;
            linMom.z = 0;
        }
        if (ext->isFixedAnyRotation()) { // if any rotation fixed, all are fixed
            if (ext->isFixedAllRotation()) {
                orient = ext->rotationQuat();
                angMom = VX3_Vec3D<double>();
            } else { // partial fixes: slow!
                VX3_Vec3D<double> tmpRotVec = orient.ToRotationVector();
                if (ext->isFixed(X_ROTATE)) {
                    tmpRotVec.x = 0;
                    angMom.x = 0;
                }
                if (ext->isFixed(Y_ROTATE)) {
                    tmpRotVec.y = 0;
                    angMom.y = 0;
                }
                if (ext->isFixed(Z_ROTATE)) {
                    tmpRotVec.z = 0;
                    angMom.z = 0;
                }
                orient.FromRotationVector(tmpRotVec);
            }
        }
    }
    //	we need to check for friction conditions here (after calculating the translation) and stop things accordingly
    if (isFloorEnabled() && floorPenetration() >= 0) {
        // we must catch a slowing voxel here since it all boils down to needing access to the dt of this timestep.
        if (isFloorStaticFriction()) {
            angMom = VX3_Vec3D<>(0, 0, 0);
        }
    }

    poissonsStrainInvalid = true;

    // if (d_signals.size()==1) printf("< %p, %f, %s, %d, %d\n", this, currentTime, k->vxa_filename, d_signals.sizeof_chunk,
    // d_signals.size());

    if (k->EnableSignals) {
        // printf("%f) before propagateSignal. this=%p.\n",currentTime, this);
        propagateSignal(currentTime);
        packMaker(currentTime);
        localSignalDecay(currentTime);
    }
}

__device__ void VX3_Voxel::localSignalDecay(double currentTime) {
    if (localSignaldt > currentTime)
        return;
    if (localSignal < 0.1) {
        // lower than threshold, simply ignore.
        localSignal = 0;
    } else {
        localSignal = localSignal * 0.9;
        localSignaldt = currentTime + 0.01;
    }
}

__device__ void VX3_Voxel::packMaker(double currentTime) {
    if (!mat->isPaceMaker)
        return;
    if (packmakerNextPulse > currentTime)
        return;

    receiveSignal(100, currentTime, true);
    packmakerNextPulse = currentTime + mat->PaceMakerPeriod;
}

__device__ void VX3_Voxel::receiveSignal(double signalValue, double activeTime, bool force) {
    if (!force) {
        if (inactiveUntil > activeTime)
            return;
    }
    if (signalValue < 0.1) {
        // lower than threshold, simply ignore.
        return;
    }

    // if received a signal, this cell will activate at activeTime, and before that, no need to receive another signal.
    inactiveUntil = activeTime + mat->inactivePeriod;

    localSignal = signalValue;
    // VX3_Signal *s = new VX3_Signal();
    d_signal.value = signalValue * mat->signalValueDecay;
    if (d_signal.value < 0.1)
        d_signal.value = 0;
    d_signal.activeTime = activeTime;
    // d_signals.push_back(s);
    // InsertSignalQueue(signalValue, currentTime + mat->signalTimeDelay);
}
__device__ void VX3_Voxel::propagateSignal(double currentTime) {
    // first one in queue, check the time
    // if (inactiveUntil > currentTime)
    //     return;
    if (d_signal.activeTime > currentTime)
        return;
    if (d_signal.value < 0.1) {
        return;
    }
    for (int i = 0; i < 6; i++) {
        if (links[i]) {
            if (links[i]->pVNeg == this) {
                links[i]->pVPos->receiveSignal(d_signal.value, currentTime + mat->signalTimeDelay, false);
            } else {
                links[i]->pVNeg->receiveSignal(d_signal.value, currentTime + mat->signalTimeDelay, false);
            }
        }
    }

    d_signal.value=0;
    d_signal.activeTime=0;
    inactiveUntil = currentTime + 2 * mat->signalTimeDelay + mat->inactivePeriod;
    // if (s)
    //     delete s;
    //     // printf("%f) delete s. this=%p. d_signals.size() %d. \n",currentTime, this, d_signals.size() );
}

__device__ VX3_Vec3D<double> VX3_Voxel::force() {

    // forces from internal bonds
    VX3_Vec3D<double> totalForce(0, 0, 0);
    for (int i = 0; i < 6; i++) {
        if (links[i])
            totalForce += links[i]->force(isNegative((linkDirection)i)); // total force in LCS
    }

    totalForce = orient.RotateVec3D(totalForce); // from local to global coordinates
    assert(!(totalForce.x != totalForce.x) || !(totalForce.y != totalForce.y) || !(totalForce.z != totalForce.z)); // assert non QNAN

    // other forces
    if (externalExists())
        totalForce += external()->force();                     // external forces
    totalForce -= velocity() * mat->globalDampingTranslateC(); // global damping f-cv
    totalForce.z += mat->gravityForce();                       // gravity, according to f=mg

    // no collision yet
    // if (isCollisionsEnabled()){
    // for (int i=0;i<colWatch.size();i++){
    // }
    // }
    totalForce -= contactForce;
    contactForce.clear();

    totalForce += CiliaForce * mat->Cilia;
    CiliaForce.clear();

    return totalForce;
}

__device__ VX3_Vec3D<double> VX3_Voxel::moment() {
    // moments from internal bonds
    VX3_Vec3D<double> totalMoment(0, 0, 0);
    for (int i = 0; i < 6; i++) {
        if (links[i]) {
            totalMoment += links[i]->moment(isNegative((linkDirection)i)); // total force in LCS
        }
    }
    totalMoment = orient.RotateVec3D(totalMoment);

    // other moments
    if (externalExists())
        totalMoment += external()->moment();                        // external moments
    totalMoment -= angularVelocity() * mat->globalDampingRotateC(); // global damping
    return totalMoment;
}

__device__ void VX3_Voxel::floorForce(float dt, VX3_Vec3D<double> *pTotalForce) {
    float CurPenetration = floorPenetration(); // for now use the average.
    if (CurPenetration >= 0) {
        VX3_Vec3D<double> vel = velocity();
        VX3_Vec3D<double> horizontalVel(vel.x, vel.y, 0);

        float normalForce = mat->penetrationStiffness() * CurPenetration;

        pTotalForce->z += normalForce - mat->collisionDampingTranslateC() * vel.z; // in the z direction: k*x-C*v - spring and damping

        if (isFloorStaticFriction()) { // If this voxel is currently in static friction mode (no lateral motion)
            assert(horizontalVel.Length2() == 0);
            float surfaceForceSq =
                (float)(pTotalForce->x * pTotalForce->x + pTotalForce->y * pTotalForce->y); // use squares to avoid a square root
            float frictionForceSq = (mat->muStatic * normalForce) * (mat->muStatic * normalForce);

            if (surfaceForceSq > frictionForceSq)
                setFloorStaticFriction(false); // if we're breaking static friction, leave the forces as they currently have been calculated
                                               // to initiate motion this time step
        } else { // even if we just transitioned don't process here or else with a complete lack of momentum it'll just go back to static
                 // friction
            *pTotalForce -=
                mat->muKinetic * normalForce * horizontalVel.Normalized(); // add a friction force opposing velocity according to the normal
                                                                           // force and the kinetic coefficient of friction
        }
    } else
        setFloorStaticFriction(false);
}

__device__ VX3_Vec3D<float> VX3_Voxel::strain(bool poissonsStrain) const {
    // if no connections in the positive and negative directions of a particular axis, strain is zero
    // if one connection in positive or negative direction of a particular axis, strain is that strain - ?? and force or constraint?
    // if connections in both the positive and negative directions of a particular axis, strain is the average.

    VX3_Vec3D<float> intStrRet(0, 0, 0); // intermediate strain return value. axes according to linkAxis enum
    int numBondAxis[3] = {0};            // number of bonds in this axis (0,1,2). axes according to linkAxis enum
    bool tension[3] = {false};
    for (int i = 0; i < 6; i++) { // cycle through link directions
        if (links[i]) {
            int axis = toAxis((linkDirection)i);
            intStrRet[axis] += links[i]->axialStrain(isNegative((linkDirection)i));
            numBondAxis[axis]++;
        }
    }
    for (int i = 0; i < 3; i++) { // cycle through axes
        if (numBondAxis[i] == 2)
            intStrRet[i] *= 0.5f; // average
        if (poissonsStrain) {
            tension[i] = ((numBondAxis[i] == 2) ||
                          (ext && (numBondAxis[i] == 1 &&
                                   (ext->isFixed((dofComponent)(1 << i)) ||
                                    ext->force()[i] != 0)))); // if both sides pulling, or just one side and a fixed or forced voxel...
        }
    }

    if (poissonsStrain) {
        if (!(tension[0] && tension[1] && tension[2])) { // if at least one isn't in tension
            float add = 0;
            for (int i = 0; i < 3; i++)
                if (tension[i])
                    add += intStrRet[i];
            float value = pow(1.0f + add, -mat->poissonsRatio()) - 1.0f;
            for (int i = 0; i < 3; i++)
                if (!tension[i])
                    intStrRet[i] = value;
        }
    }

    return intStrRet;
}

__device__ VX3_Vec3D<float> VX3_Voxel::poissonsStrain() {
    if (poissonsStrainInvalid) {
        pStrain = strain(true);
        poissonsStrainInvalid = false;
    }
    return pStrain;
}

__device__ float VX3_Voxel::transverseStrainSum(linkAxis axis) {
    if (mat->poissonsRatio() == 0)
        return 0;

    VX3_Vec3D<float> psVec = poissonsStrain();

    switch (axis) {
    case X_AXIS:
        return psVec.y + psVec.z;
    case Y_AXIS:
        return psVec.x + psVec.z;
    case Z_AXIS:
        return psVec.x + psVec.y;
    default:
        return 0.0f;
    }
}

__device__ float VX3_Voxel::transverseArea(linkAxis axis) {
    float size = (float)mat->nominalSize();
    if (mat->poissonsRatio() == 0)
        return size * size;

    VX3_Vec3D<> psVec = poissonsStrain();

    switch (axis) {
    case X_AXIS:
        return (float)(size * size * (1 + psVec.y) * (1 + psVec.z));
    case Y_AXIS:
        return (float)(size * size * (1 + psVec.x) * (1 + psVec.z));
    case Z_AXIS:
        return (float)(size * size * (1 + psVec.x) * (1 + psVec.y));
    default:
        return size * size;
    }
}

__device__ void VX3_Voxel::updateSurface() {
    bool interior = true;
    for (int i = 0; i < 6; i++)
        if (!links[i]) {
            interior = false;
        } else if (links[i]->isDetached) {
            interior = false;
        }
    interior ? boolStates |= SURFACE : boolStates &= ~SURFACE;
}

__device__ void VX3_Voxel::enableCollisions(bool enabled, float watchRadius) {
    enabled ? boolStates |= COLLISIONS_ENABLED : boolStates &= ~COLLISIONS_ENABLED;
}

__device__ void VX3_Voxel::generateNearby(int linkDepth, int gindex, bool surfaceOnly) {
    assert(false); // not used. near by has logic flaws.
}

#include "VX3_Link.h"
#include "VX3_MaterialLink.h"
#include "VX3_VoxelyzeKernel.cuh"
#include <vector>

VX3_Link::VX3_Link(CVX_Link *p, VX3_VoxelyzeKernel *k)
    : forceNeg(p->forceNeg), forcePos(p->forcePos), momentNeg(p->momentNeg), momentPos(p->momentPos), strain(p->strain),
      maxStrain(p->maxStrain), strainOffset(p->strainOffset), boolStates(p->boolStates), axis(p->axis), strainRatio(p->strainRatio),
      pos2(p->pos2), angle1v(p->angle1v), angle2v(p->angle2v), angle1(p->angle1), angle2(p->angle2), smallAngle(p->smallAngle),
      currentRestLength(p->currentRestLength), currentTransverseArea(p->currentTransverseArea),
      currentTransverseStrainSum(p->currentTransverseStrainSum), _stress(p->_stress) {
    pVNeg = k->h_lookup_voxels[p->pVNeg];
    pVPos = k->h_lookup_voxels[p->pVPos];
    // for (int i=0;i<k->num_d_voxels;i++) {
    // 	if (k->h_voxels[i] == p->pVNeg) {
    // 		pVNeg = &k->d_voxels[i];
    // 		assert(pVNeg == k->h_lookup_voxels[p->pVNeg]);
    // 	}
    // 	if (k->h_voxels[i] == p->pVPos) {
    // 		pVPos = &k->d_voxels[i];
    // 	}
    // }

    for (int i = 0; i < k->num_d_linkMats; i++) {
        if (k->h_linkMats[i] == p->mat) {
            mat = &k->d_linkMats[i];
        }
    }
}
// VX3_Link can also be initialized in device
__device__ VX3_Link::VX3_Link(VX3_Voxel *voxel1, linkDirection dir1, VX3_Voxel *voxel2, linkDirection dir2, linkAxis link_axis,
                              VX3_VoxelyzeKernel *k) {
    // TODO: is this the right orientation? let's see...

    voxel1->links[dir1] = this;
    voxel2->links[dir2] = this;
    axis = link_axis;
    pVNeg = voxel2;
    pVPos = voxel1;
    // for (int i=0;i<6;i++) {
    // 	if (voxel1->links[i]==NULL) {
    // 		voxel1->links[i]=this;
    // 		break;
    // 	}
    // }
    // for (int i=0;i<6;i++) {
    // 	if (voxel2->links[i]==NULL) {
    // 		voxel2->links[i]=this;
    // 		break;
    // 	}
    // }
    //
    mat = k->combinedMaterial(voxel1->material(), voxel2->material());
    boolStates = 0;
    reset();
}

__device__ void VX3_Link::reset() {
    pos2 = angle1v = angle2v = VX3_Vec3D<double>();
    angle1 = angle2 = VX3_Quat3D<double>();
    forceNeg = forcePos = momentNeg = momentPos = VX3_Vec3D<double>();
    strain = maxStrain = strainOffset = _stress = 0.0f;
    strainRatio = pVPos->material()->E / pVNeg->material()->E;
    smallAngle = true;

    setBoolState(LOCAL_VELOCITY_VALID, false);

    updateRestLength();
    updateTransverseInfo();
}

__device__ float VX3_Link::axialStrain(bool positiveEnd) const {
    return positiveEnd ? 2.0f * strain * strainRatio / (1.0f + strainRatio) : 2.0f * strain / (1.0f + strainRatio);
}

__device__ bool VX3_Link::isYielded() const { return mat->isYielded(maxStrain); }

__device__ bool VX3_Link::isFailed() const { return mat->isFailed(maxStrain); }

__device__ void VX3_Link::updateRestLength() {
    // update rest length according to temperature of both end
    currentRestLength = 0.5 * (pVNeg->baseSize(axis) + pVPos->baseSize(axis));
}

__device__ void VX3_Link::updateTransverseInfo() {
    currentTransverseArea = 0.5f * (pVNeg->transverseArea(axis) + pVPos->transverseArea(axis));
    currentTransverseStrainSum = 0.5f * (pVNeg->transverseStrainSum(axis) + pVPos->transverseStrainSum(axis));
}

__device__ VX3_Quat3D<double> VX3_Link::orientLink() // updates pos2, angle1, angle2, and smallAngle
{
    VX3_Vec3D<> _pos2 = pVPos->position() - pVNeg->position();
    pos2 = toAxisX(_pos2); // digit truncation happens here...
    VX3_Quat3D<> _angle1 = pVNeg->orientation();
    angle1 = toAxisX(_angle1);
    VX3_Quat3D<> _angle2 = pVPos->orientation();
    angle2 = toAxisX(_angle2);

    VX3_Quat3D<double> totalRot = angle1.Conjugate(); // keep track of the total rotation of this bond
                                                      // (after toAxisX())
    pos2 = totalRot.RotateVec3D(pos2);
    angle2 = totalRot * angle2;
    angle1 = VX3_Quat3D<>(); // zero for now...

    // small angle approximation?
    float SmallTurn = (float)((abs(pos2.z) + abs(pos2.y)) / pos2.x);
    float ExtendPerc = (float)(abs(1 - pos2.x / currentRestLength));
    if (!smallAngle /*&& angle2.IsSmallAngle()*/ && SmallTurn < SA_BOND_BEND_RAD && ExtendPerc < SA_BOND_EXT_PERC) {
        smallAngle = true;
        setBoolState(LOCAL_VELOCITY_VALID, false);
    } else if (smallAngle && (/*!angle2.IsSmallishAngle() || */ SmallTurn > HYSTERESIS_FACTOR * SA_BOND_BEND_RAD ||
                              ExtendPerc > HYSTERESIS_FACTOR * SA_BOND_EXT_PERC)) {
        smallAngle = false;
        setBoolState(LOCAL_VELOCITY_VALID, false);
    }

    if (smallAngle) {                 // Align so Angle1 is all zeros
        pos2.x -= currentRestLength;  // only valid for small angles
    } else {                          // Large angle. Align so that Pos2.y, Pos2.z are zero.
        angle1.FromAngleToPosX(pos2); // get the angle to align Pos2 with the X axis
        totalRot = angle1 * totalRot; // update our total rotation to reflect this
        angle2 = angle1 * angle2;     // rotate angle2
        pos2 = VX3_Vec3D<>(pos2.Length() - currentRestLength, 0, 0);
    }

    angle1v = angle1.ToRotationVector();
    angle2v = angle2.ToRotationVector();

    assert(!(angle1v.x != angle1v.x) || !(angle1v.y != angle1v.y) || !(angle1v.z != angle1v.z)); // assert non QNAN
    assert(!(angle2v.x != angle2v.x) || !(angle2v.y != angle2v.y) || !(angle2v.z != angle2v.z)); // assert non QNAN

    return totalRot;
}

__device__ void VX3_Link::updateForces() {
    VX3_Vec3D<double> oldPos2 = pos2;
    VX3_Vec3D<double> oldAngle1v = angle1v;
    VX3_Vec3D<double> oldAngle2v = angle2v; // remember the positions/angles from last timestep to
                                            // calculate velocity

    orientLink();                                     // sets pos2, angle1, angle2
    VX3_Vec3D<double> dPos2 = 0.5 * (pos2 - oldPos2); // deltas for local damping. velocity at center
                                                      // is half the total velocity
    VX3_Vec3D<double> dAngle1 = 0.5 * (angle1v - oldAngle1v);
    VX3_Vec3D<double> dAngle2 = 0.5 * (angle2v - oldAngle2v);
    // if volume effects..
    if (!mat->isXyzIndependent() || currentTransverseStrainSum != 0) { // currentTransverseStrainSum != 0 catches when we disable
                                                                       // poissons mid-simulation
        // updateTransverseInfo();
    }
    _stress = updateStrain((float)(pos2.x / currentRestLength));
    if (isFailed()) {
        forceNeg = forcePos = momentNeg = momentPos = VX3_Vec3D<double>(0, 0, 0);
        return;
    }
    float b1 = mat->_b1, b2 = mat->_b2, b3 = mat->_b3,
          a2 = mat->_a2; // local copies
    // Beam equations. All relevant terms are here, even though some are zero
    // for small angle and others are zero for large angle (profiled as
    // negligible performance penalty)
    forceNeg = VX3_Vec3D<double>(_stress * currentTransverseArea, // currentA1*pos2.x,
                                 b1 * pos2.y - b2 * (angle1v.z + angle2v.z),
                                 b1 * pos2.z + b2 * (angle1v.y + angle2v.y)); // Use Curstress instead of -a1*Pos2.x
                                                                              // to account for non-linear deformation
    forcePos = -forceNeg;

    momentNeg = VX3_Vec3D<double>(a2 * (angle2v.x - angle1v.x), -b2 * pos2.z - b3 * (2 * angle1v.y + angle2v.y),
                                  b2 * pos2.y - b3 * (2 * angle1v.z + angle2v.z));
    momentPos = VX3_Vec3D<double>(a2 * (angle1v.x - angle2v.x), -b2 * pos2.z - b3 * (angle1v.y + 2 * angle2v.y),
                                  b2 * pos2.y - b3 * (angle1v.z + 2 * angle2v.z));
    // local damping:
    if (isLocalVelocityValid()) { // if we don't have the basis for a good
                                  // damping calculation, don't do any damping.

        float sqA1 = mat->_sqA1, sqA2xIp = mat->_sqA2xIp, sqB1 = mat->_sqB1, sqB2xFMp = mat->_sqB2xFMp, sqB3xIp = mat->_sqB3xIp;
        VX3_Vec3D<double> posCalc(sqA1 * dPos2.x, sqB1 * dPos2.y - sqB2xFMp * (dAngle1.z + dAngle2.z),
                                  sqB1 * dPos2.z + sqB2xFMp * (dAngle1.y + dAngle2.y));

        forceNeg += pVNeg->dampingMultiplier() * posCalc;
        forcePos -= pVPos->dampingMultiplier() * posCalc;

        momentNeg -= 0.5 * pVNeg->dampingMultiplier() *
                     VX3_Vec3D<>(-sqA2xIp * (dAngle2.x - dAngle1.x), sqB2xFMp * dPos2.z + sqB3xIp * (2 * dAngle1.y + dAngle2.y),
                                 -sqB2xFMp * dPos2.y + sqB3xIp * (2 * dAngle1.z + dAngle2.z));
        momentPos -= 0.5 * pVPos->dampingMultiplier() *
                     VX3_Vec3D<>(sqA2xIp * (dAngle2.x - dAngle1.x), sqB2xFMp * dPos2.z + sqB3xIp * (dAngle1.y + 2 * dAngle2.y),
                                 -sqB2xFMp * dPos2.y + sqB3xIp * (dAngle1.z + 2 * dAngle2.z));

    } else
        setBoolState(LOCAL_VELOCITY_VALID,
                     true); // we're good for next go-around unless something changes
    //	transform forces and moments to local voxel coordinates
    if (!smallAngle) {
        forceNeg = angle1.RotateVec3DInv(forceNeg);
        momentNeg = angle1.RotateVec3DInv(momentNeg);
    }
    forcePos = angle2.RotateVec3DInv(forcePos);
    momentPos = angle2.RotateVec3DInv(momentPos);

    toAxisOriginal(&forceNeg);
    toAxisOriginal(&forcePos);
    toAxisOriginal(&momentNeg);
    toAxisOriginal(&momentPos);
    // printf("momentNeg: (%f, %f, %f).\n", momentNeg.x, momentNeg.y,
    // momentNeg.z);
    if (isNewLink) {
        // for debug
        forceNeg = forceNeg * 0.01;
        forcePos = forcePos * 0.01;
        momentNeg = momentNeg * 0.01;
        momentPos = momentPos * 0.01;
        isNewLink -= 1;
    }
    // assert(!(forceNeg.x != forceNeg.x) || !(forceNeg.y != forceNeg.y) ||
    // !(forceNeg.z != forceNeg.z)); //assert non QNAN assert(!(forcePos.x !=
    // forcePos.x) || !(forcePos.y != forcePos.y) || !(forcePos.z !=
    // forcePos.z)); //assert non QNAN
}

__device__ float VX3_Link::updateStrain(float axialStrain) {
    // int di = 0;

    strain = axialStrain; // redundant?

    if (mat->linear) {

        if (axialStrain > maxStrain)
            maxStrain = axialStrain; // remember this maximum for easy reference

        return mat->stress(axialStrain, currentTransverseStrainSum);
    } else {
        float returnStress;

        if (axialStrain > maxStrain) { // if new territory on the stress/strain curve
            maxStrain = axialStrain;   // remember this maximum for easy reference
            returnStress = mat->stress(axialStrain, currentTransverseStrainSum);

            if (mat->nu != 0.0f)
                strainOffset = maxStrain - mat->stress(axialStrain) / (mat->_eHat * (1 - mat->nu)); // precalculate strain offset for when
                                                                                                    // we back off
            else
                strainOffset = maxStrain - returnStress / mat->E; // precalculate strain offset for
                                                                  // when we back off

        } else { // backed off a non-linear material, therefore in linear
                 // region.

            float relativeStrain = axialStrain - strainOffset; // treat the material as linear with
                                                               // a strain offset according to the
                                                               // maximum plastic deformation

            if (mat->nu != 0.0f)
                returnStress = mat->stress(relativeStrain, currentTransverseStrainSum, true);
            else
                returnStress = mat->E * relativeStrain;
        }

        return returnStress;
    }
}

__device__ float VX3_Link::strainEnergy() const {
    return forceNeg.x * forceNeg.x / (2.0f * mat->_a1) +                                                            // Tensile strain
           momentNeg.x * momentNeg.x / (2.0 * mat->_a2) +                                                           // Torsion strain
           (momentNeg.z * momentNeg.z - momentNeg.z * momentPos.z + momentPos.z * momentPos.z) / (3.0 * mat->_b3) + // Bending Z
           (momentNeg.y * momentNeg.y - momentNeg.y * momentPos.y + momentPos.y * momentPos.y) / (3.0 * mat->_b3);  // Bending Y
}

__device__ float VX3_Link::axialStiffness() {
    if (mat->isXyzIndependent())
        return mat->_a1;
    else {
        updateRestLength();
        updateTransverseInfo();

        return (float)(mat->_eHat * currentTransverseArea / ((strain + 1) * currentRestLength)); // _a1;
    }
}

__device__ float VX3_Link::a1() const { return mat->_a1; }
__device__ float VX3_Link::a2() const { return mat->_a2; }
__device__ float VX3_Link::b1() const { return mat->_b1; }
__device__ float VX3_Link::b2() const { return mat->_b2; }
__device__ float VX3_Link::b3() const { return mat->_b3; }

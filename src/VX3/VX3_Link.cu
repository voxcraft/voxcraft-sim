#include "math_constants.h"

#include "VX3_Link.h"
#include "VX3_MaterialLink.h"
#include "VX3_VoxelyzeKernel.cuh"
#include <vector>

VX3_Link::VX3_Link(CVX_Link *p, VX3_VoxelyzeKernel *k)
    : forceNeg(p->forceNeg), forcePos(p->forcePos), momentNeg(p->momentNeg), momentPos(p->momentPos), strain(p->strain),
      maxStrain(p->maxStrain), strainOffset(p->strainOffset), boolStates(p->boolStates), strainRatio(p->strainRatio),
      pos2(p->pos2), angle1v(p->angle1v), angle2v(p->angle2v), angle1(p->angle1), angle2(p->angle2), smallAngle(p->smallAngle),
      currentRestLength(p->currentRestLength), currentTransverseArea(p->currentTransverseArea),
      currentTransverseStrainSum(p->currentTransverseStrainSum), _stress(p->_stress) {
    pVNeg = k->h_lookup_voxels[p->pVNeg];
    pVPos = k->h_lookup_voxels[p->pVPos];
    for (int i = 0; i < k->num_d_linkMats; i++) {
        if (k->h_linkMats[i] == p->mat) {
            mat = &k->d_linkMats[i];
        }
    }
}
// VX3_Link can also be initialized in device
__device__ VX3_Link::VX3_Link(VX3_Voxel *voxelNeg, linkDirection dirNeg, VX3_Voxel *voxelPos, linkDirection dirPos, 
                                linkAxis link_axis, linkAxis link_axis_pos, VX3_VoxelyzeKernel *k) {
    syncVectors(k);
    voxelNeg->links[dirNeg] = this;
    voxelPos->links[dirPos] = this;
    pVNeg = voxelNeg;
    pVPos = voxelPos;
    linkdirNeg = dirNeg;
    linkdirPos = dirPos;
    
    mat = k->combinedMaterial(voxelNeg->material(), voxelPos->material());
    boolStates = 0;

    reset();
}

__device__ void VX3_Link::syncVectors(VX3_VoxelyzeKernel *k) {
    d_kernel = k;

    double num71 = cos(CUDART_PIO4);
    quat_linkDirection[(int)X_POS] = VX3_Quat3D<double>(1,0,0,0);
    quat_linkDirection[(int)X_NEG] = VX3_Quat3D<double>(0,0,0,-1);
    quat_linkDirection[(int)Y_POS] = VX3_Quat3D<double>(num71,0,0,num71);
    quat_linkDirection[(int)Y_NEG] = VX3_Quat3D<double>(num71,0,0,-num71);
    quat_linkDirection[(int)Z_POS] = VX3_Quat3D<double>(num71,0,-num71,0);
    quat_linkDirection[(int)Z_NEG] = VX3_Quat3D<double>(num71,0,num71,0);
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
    currentRestLength = 0.5 * (pVNeg->baseSize(toAxis(linkdirNeg)) + pVPos->baseSize(toAxis(linkdirPos)));
}

__device__ void VX3_Link::updateTransverseInfo() {
    currentTransverseArea = 0.5f * (pVNeg->transverseArea(toAxis(linkdirNeg)) + pVPos->transverseArea(toAxis(linkdirPos)));
    currentTransverseStrainSum = 0.5f * (pVNeg->transverseStrainSum(toAxis(linkdirNeg)) + pVPos->transverseStrainSum(toAxis(linkdirPos)));
}

/*
    Important Notes [Sida's Best Guess]:

    refer to: VX3_Link_orientLink.jpg

    In local coordinates, the link is always in the X_POS direction. pVNeg--->pVPos

    raw_pos2: position of pVPos relative to pVNeg in raw coordinates.
    local_pos2: position of pVPos relative to pVNeg in local coordinates.
    raw_angle1: orientation of pVNeg in raw coordinates.
    raw_angle2: orientation of pVPos in raw coordinates.
    totalRotation: a Quat rotate from real (raw) coordinates to local coordinates.
    pos2: position of pVPos in local corrdinates.
    angle1: pVNeg's orientation in local coordinates.
    angle2: pVPos's orientation in local coordinates
    final pos2: position of pVPos relative to rest place in local coordinates.

    About Quat3D: When the observer look from -inf to +inf along x axis, adding x term will causing clock-wise rotation.
 */
__device__ void VX3_Link::orientLink_new() // updates pos2, angle1, angle2, and smallAngle
{
    if (isNewLink) {
        printf("Debug.\n");
    }


    VX3_Quat3D<> totalRotation = quat_linkDirection[linkdirNeg] * pVNeg->orientation(); // (3) First Rotate according to pVNeg->orientation (on the right), then Rotate according to quat_linkDirection[linkdirNeg] (on the left).
    VX3_Vec3D<> raw_pos2 = pVPos->position() - pVNeg->position(); //(4)
    VX3_Quat3D<> raw_angle2 = quat_linkDirection[X_NEG] * quat_linkDirection[linkdirPos] * pVPos->orientation(); // (5,6,7)
    pos2 = totalRotation.RotateVec3DInv(raw_pos2); //(8)
    debug++;
    angle2 = totalRotation.Conjugate() * raw_angle2; // (9)
    angle1 = VX3_Quat3D<>(); // always (1,0,0,0)

    //SmallTurn: To what extent does the beam bent.
    //ExtendPerc: The extend percentage as a spring. 
    //smallAngle: The beam act like a spring in X-axis in local coordinates.
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
        angle2 = angle1 * angle2;     // rotate angle2
        pos2 = VX3_Vec3D<>(pos2.Length() - currentRestLength, 0, 0);
    }

    angle1v = angle1.ToRotationVector();
    angle2v = angle2.ToRotationVector();

    assert(!(angle1v.x != angle1v.x) || !(angle1v.y != angle1v.y) || !(angle1v.z != angle1v.z)); // assert non QNAN
    assert(!(angle2v.x != angle2v.x) || !(angle2v.y != angle2v.y) || !(angle2v.z != angle2v.z)); // assert non QNAN

}

__device__ void VX3_Link::orientLink() // updates pos2, angle1, angle2, and smallAngle
{
    VX3_Vec3D<> tmp_pos2, old_pos2;
    VX3_Quat3D<> tmp_angle2, old_angle2;
    if (true) {
        // 1. pVNeg->orientation
        // 2. linkdirNeg
        // 1+2 => 3. totalRotation
        // 4. raw pos2
        // 5. pVPos->orientation
        // 6. linkdirPos
        // 5+6 => 7. raw angle2
        // 3+4 => 8. pos2
        // 3+7 => 9. angle2
        // Imagine that pVNeg is placed on the origin, and the linkdirNeg is pointing towards +X direction.
        // (1) raw pos2: position of pVPos relative to pVNeg in real coordinates.
        // (2) axis pos2: Orthogonally rotate the position, makes it close to the X_POS direction of pVNeg.
        //                The reason not do rotation together is that orthogonal rotation can be performed as swapping axis, instead of using Quaternion, which will involve Pi and Cos. (More accurate)
        // (3) pos2: The position is in the imaginary coordinates. (Must be close to (1,0,0))
        // (4) raw angle 1: The orientation (X_POS) of pVNeg in real coordinates.
        // (5) angle 1: The orientation of linkdirNeg in real coordinates.
        VX3_Vec3D<> raw_pos2 = pVPos->position() - pVNeg->position(); //(1)
        VX3_Quat3D<> rotate_pos2 = quat_linkDirection[linkdirNeg].Conjugate() * pVNeg->orientation().Conjugate(); // (A)
        // (A) is the same as rotating twice:
        // pos2 = pVNeg->orientation().RotateVec3DInv(raw_pos2);
        // pos2 = quat_linkDirection[linkdirNeg].RotateVec3DInv(pos2);
        pos2 = rotate_pos2.RotateVec3D(raw_pos2);

        VX3_Quat3D<> raw_angle2 = pVPos->orientation(); // (5,6,7)
        angle2 = quat_linkDirection[linkdirPos] * raw_angle2;
        angle2 = rotate_pos2 * angle2;
        // VX3_Quat3D<> totalRotation = quat_linkDirection[linkdirNeg] * pVNeg->orientation(); // (3) First Rotate according to pVNeg->orientation (on the right), then Rotate according to quat_linkDirection[linkdirNeg] (on the left).
        // totalRotation = totalRotation.Conjugate();
        // totalRotation = totalRotation * quat_linkDirection[linkdirNeg];
        // pos2 = totalRotation.RotateVec3D(raw_pos2); //(8)
        // debug++;
        // angle2 = totalRotation * raw_angle2; // (9)
        // angle1 = VX3_Quat3D<>(); // always (1,0,0,0)
        old_angle2 = raw_angle2;
        tmp_angle2 = angle2;
    }
    if (true) {
        VX3_Vec3D<> _pos2 = pVPos->position() - pVNeg->position();
        pos2 = toAxisX(_pos2, toAxis(linkdirNeg)); // digit truncation happens here...
        VX3_Quat3D<> _angle1 = pVNeg->orientation();
        angle1 = toAxisX(_angle1, toAxis(linkdirNeg));
        VX3_Quat3D<> _angle2 = pVPos->orientation();
        angle2 = toAxisX(_angle2, toAxis((linkDirection)oppositeDir(linkdirPos)));

        VX3_Quat3D<double> totalRot = angle1.Conjugate(); // keep track of the total rotation of this bond
                                                        // (after toAxisX())
        pos2 = totalRot.RotateVec3D(pos2);
        angle2 = totalRot * angle2;
        angle1 = VX3_Quat3D<>(); // zero for now...
        if (d_kernel->currentTime>=0.08) {
            printf("%d) angle2: %e, %e, %e, %e \n \t %e, %e, %e, %e \n\t\t\t\t\t (old) %e, %e, %e, %e.\n", linkdirNeg,
            tmp_angle2.w, tmp_angle2.x, tmp_angle2.y, tmp_angle2.z,
            angle2.w, angle2.x, angle2.y, angle2.z,
            old_angle2.w, old_angle2.x, old_angle2.y, old_angle2.z
        );
        }      
    }

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
        angle2 = angle1 * angle2;     // rotate angle2
        pos2 = VX3_Vec3D<>(pos2.Length() - currentRestLength, 0, 0);
    }

    angle1v = angle1.ToRotationVector();
    angle2v = angle2.ToRotationVector();

    assert(!(angle1v.x != angle1v.x) || !(angle1v.y != angle1v.y) || !(angle1v.z != angle1v.z)); // assert non QNAN
    assert(!(angle2v.x != angle2v.x) || !(angle2v.y != angle2v.y) || !(angle2v.z != angle2v.z)); // assert non QNAN
}

__device__ void VX3_Link::updateForces() {
    if (isNewLink) {
        isNewLink -= 1;
    }
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
        float neg_damp = pVNeg->dampingMultiplier();
        float pos_damp = pVPos->dampingMultiplier();

        float sqA1 = mat->_sqA1, sqA2xIp = mat->_sqA2xIp, sqB1 = mat->_sqB1, sqB2xFMp = mat->_sqB2xFMp, sqB3xIp = mat->_sqB3xIp;
        VX3_Vec3D<double> posCalc(sqA1 * dPos2.x, sqB1 * dPos2.y - sqB2xFMp * (dAngle1.z + dAngle2.z),
                                  sqB1 * dPos2.z + sqB2xFMp * (dAngle1.y + dAngle2.y));

        forceNeg += neg_damp * posCalc;
        forcePos -= pos_damp * posCalc;

        momentNeg -= 0.5 * neg_damp *
                     VX3_Vec3D<>(-sqA2xIp * (dAngle2.x - dAngle1.x), sqB2xFMp * dPos2.z + sqB3xIp * (2 * dAngle1.y + dAngle2.y),
                                 -sqB2xFMp * dPos2.y + sqB3xIp * (2 * dAngle1.z + dAngle2.z));
        momentPos -= 0.5 * pos_damp *
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

    // Rewrite Rotation back, so linkNeg and linkPos don't need to be paired. (For arbitrary attachment.)
    if (true) {
        // new method
        forceNeg = quat_linkDirection[linkdirNeg].RotateVec3D(forceNeg);
        forcePos = quat_linkDirection[oppositeDir(linkdirPos)].RotateVec3D(forcePos);
        momentNeg = quat_linkDirection[linkdirNeg].RotateVec3D(momentNeg);
        momentPos = quat_linkDirection[oppositeDir(linkdirPos)].RotateVec3D(momentPos);
    } 
    if (false) {
        // old method
        toAxisOriginal(&forceNeg, toAxis(linkdirNeg));
        toAxisOriginal(&forcePos, toAxis((linkDirection)oppositeDir(linkdirPos)));
        toAxisOriginal(&momentNeg, toAxis(linkdirNeg));
        toAxisOriginal(&momentPos, toAxis((linkDirection)oppositeDir(linkdirPos)));
    }

    assert(!(forceNeg.x != forceNeg.x) || !(forceNeg.y != forceNeg.y) || !(forceNeg.z != forceNeg.z)); //assert non QNAN
    assert(!(forcePos.x != forcePos.x) || !(forcePos.y != forcePos.y) || !(forcePos.z != forcePos.z)); //assert non QNAN
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

__device__ __host__ int VX3_Link::oppositeDir(int linkdir) {
    // X_NEG for X_POS, etc.
    int residual = linkdir % 2;
    return linkdir - residual + !(linkdir % 2);
}

// Because using quaternion to rotate exact 90 degree is impossible, instead, we use this method to do rotation along axis
__device__ VX3_Vec3D<double> VX3_Link::pseudoRotation(linkDirection linkdir, VX3_Vec3D<double> pos, bool inverse){

}

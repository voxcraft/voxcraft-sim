/*******************************************************************************
Copyright (c) 2010, Jonathan Hiller (Cornell University)
If used in publication cite "J. Hiller and H. Lipson "Dynamic Simulation of Soft Heterogeneous Objects" In press. (2011)"

This file is part of Voxelyze.
Voxelyze is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Voxelyze is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
See <http://www.opensource.org/licenses/lgpl-3.0.html> for license details.
*******************************************************************************/

#ifndef VX_ENUMS_H
#define VX_ENUMS_H

//Enums for Voxelyze

////SIMULATION
////! Defines the type of view
//enum ViewMode {
//	RVM_NONE, //!< Disable drawing for maximum simulation speed
//	RVM_VOXELS, //!< Draws the actual voxels of the simulation
//	RVM_BONDS //!< Draws the bonds underlying the simulation
//};
////! Defines the color to display each voxel
//enum ViewColor {
//	RVC_TYPE, //!< Defined material color for voxels, type of bond for bonds
//	RVC_KINETIC_EN, //!< Kinetic energy
//	RVC_DISP, //!< Color coded displacement from inital position.
//	RVC_STATE, //Color coded information about the state of each voxel (reserved for future implementation)
//	RVC_STRAIN_EN, //!< Strain energy
//	RVC_STRAIN, //!< Color coded engineering srain, or internal deformation percent.
//	RVC_STRESS //!< Color coded engineerign stress, or internal force
//};
////! Defines the way  to display each voxel
//enum ViewVoxel {
//	RVV_DISCRETE, //!< Draws discrete voxels with position and orientation
//	RVV_DEFORMED, //!< Draws deformed voxels
//	RVV_SMOOTH //!< Draws a smooth mesh, whether provided or generated with marching cubes
//};

//! The collision detection scheme to use if enabled.
enum ColSystem{
	COL_BASIC, //!< Basic n2 collision detection. Extremely slow.
	COL_SURFACE, //!< Basic n2 collision detection, but only considering voxels on the surface of the object. Very slow.
	COL_BASIC_HORIZON, //!< Hierarchical collision detection between all voxels. (One level) Updates potential collision list only when aggregated motion requires it.
	COL_SURFACE_HORIZON //!< Hierarchical collision detection between surface voxels. (One level) This is the fastest collision scheme implemented, and suffice for most cases. Updates potential collision list only when aggregated motion requires it.
};

enum MatBlendModel{
	MB_LINEAR, //blended materials combine stiffness linearly (x)
	MB_EXPONENTIAL, //blended materials combine stiffness exponentially (2^x-1) (-1 puts in Y range of 0 to 1)
	MB_POLYNOMIAL //blended materials combine stiffness polynomially (x^n) with n stored seperately
};

//!Determines optional condition for the simulation to stop.
enum StopCondition {
	SC_NONE, //!<Runs indefinitely
	SC_MAX_TIME_STEPS, //!<Runs to a set maximum number of timesteps
	SC_MAX_SIM_TIME, //!<Runs for a set number of simulation seconds
	SC_TEMP_CYCLES, //!<IF temperature is varying, runs for a set number of cycles. Otherwise runs indefinitely.
	SC_CONST_MAXENERGY, //!<runs until kinetic+potential energy stabilizes to threshhold. Begins checking after 50 time steps, energy must stay within StopConditionValue threshold for 10 consecutive readings at 50 simulation steps apart.
	SC_MIN_KE, //!<runs until kinetic energy is below a threshhold. Begins checking after 10 time steps, energy must stay below StopConditionValue threshold for 10 consecutive readings at 50 simulation steps apart.
	SC_MIN_MAXMOVE //!<runs until maximum voxel displacement/timestep (of any voxel in the simulation) is below a threshhold (in mm?)
};

enum Axis {  //which axis do we refer to?
	AXIS_NONE,
	AXIS_X,
	AXIS_Y,
	AXIS_Z
};


//What statistics to calculate
enum SimStat : int {
	CALCSTAT_NONE = 0,
	CALCSTAT_ALL = 0xffff,
	CALCSTAT_COM = 1<<0,
	CALCSTAT_DISP = 1<<1,
	CALCSTAT_VEL = 1<<2,
	CALCSTAT_KINE = 1<<3,
	CALCSTAT_STRAINE = 1<<4,
	CALCSTAT_ENGSTRAIN = 1<<5,
	CALCSTAT_ENGSTRESS = 1<<6,
	CALCSTAT_PRESSURE = 1<<7

};

//Simulation features that can be turned on and off
enum SimFeature : int {
	VXSFEAT_NONE = 0,
	VXSFEAT_COLLISIONS = 1<<0,
	VXSFEAT_GRAVITY = 1<<1,
	VXSFEAT_FLOOR = 1<<2,
	VXSFEAT_TEMPERATURE = 1<<3,
	VXSFEAT_TEMPERATURE_VARY = 1<<4,
	VXSFEAT_PLASTICITY = 1<<5,
	VXSFEAT_FAILURE = 1<<6,
	VXSFEAT_BLENDING = 1<<7,
	VXSFEAT_VOLUME_EFFECTS = 1<<8,
	VXSFEAT_MAX_VELOCITY = 1<<9,
	VXSFEAT_EQUILIBRIUM_MODE = 1<<10,
	VXSFEAT_IGNORE_ANG_DAMP = 1<<11

};

//VOXELS

//old
enum BondDir {	//Defines the direction of a link
	BD_PX=0,			//Positive X direction
	BD_NX=1,			//Negative X direction
	BD_PY=2,			//Positive Y direction
	BD_NY=3,			//Negative Y direction
	BD_PZ=4,			//Positive Z direction
	BD_NZ=5				//Negative Z direction
}; 
#define NO_BOND -1 //if there is no bond present on a specified direction


//New
//
//enum linkDirection {	//Defines the direction of a link
//	X_POS=0,			//Positive X direction
//	X_NEG=1,			//Negative X direction
//	Y_POS=2,			//Positive Y direction
//	Y_NEG=3,			//Negative Y direction
//	Z_POS=4,			//Positive Z direction
//	Z_NEG=5				//Negative Z direction
//}; 
//enum linkAxis {
//	X_AXIS = 0,
//	Y_AXIS = 1,
//	Z_AXIS = 2
//};
////a couple global convenience functions to have wherever the link enums are used
//inline linkAxis toAxis(linkDirection direction) {return (linkAxis)((int)direction/2);}
//inline linkDirection toDirection(linkAxis axis, bool positiveDirection) {return (linkDirection)(2*((int)axis) + positiveDirection?0:1);}
//inline bool isNegative(linkDirection direction) {return direction%2==1;}
//inline bool isPositive(linkDirection direction) {return direction%2==0;}
//inline linkDirection toOpposite(linkDirection direction) {return (linkDirection)(direction-direction%2 + (direction+1)%2);}
//
//typedef unsigned char dofObject;  //Bits: 0 0 Tz, Ty, Tz, Z, Y, X. 0 if free, 1 if fixed
//enum dofComponent {
//	X_TRANSLATE=1<<0,
//	Y_TRANSLATE=1<<1,
//	Z_TRANSLATE=1<<2,
//	X_ROTATE=1<<3,
//	Y_ROTATE=1<<4,
//	Z_ROTATE=1<<5
//};
//inline void dofSet(dofObject& obj, dofComponent dof, bool set) {set ? obj|=dof : obj&=~dof;}
//inline void dofSetAll(dofObject& obj, bool set) {set ? obj|=0x3F : obj&=~0x3F;}
//inline bool dofIsSet(dofObject obj, dofComponent dof){return (dof&obj)?true:false;}
//inline bool dofIsAllSet(dofObject obj){return (obj&0x3F)==0x3F;}
//inline bool dofIsNoneSet(dofObject obj){return !(obj&0x3F);}
//inline dofObject dof(bool tx, bool ty, bool tz, bool rx, bool ry, bool rz) {dofObject ret=0; dofSet(ret, X_TRANSLATE, tx); dofSet(ret, Y_TRANSLATE, ty); dofSet(ret, Z_TRANSLATE, tz); dofSet(ret, X_ROTATE, rx); dofSet(ret, Y_ROTATE, ry); dofSet(ret, Z_ROTATE, rz); return ret;}
//

#endif //VX_ENUMS_H

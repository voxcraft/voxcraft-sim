#if !defined(TYPE_H)
#define TYPE_H

#include "debug_host_dev.h"

enum linkFlags { //default of each should be zero for easy clearing
    LOCAL_VELOCITY_VALID = 1<<0 //has something changes to render local velocity calculations (for damping) invalid? (this happens when small angle or global base size has changed)
};

typedef int linkState;

//! Defines the direction of a link relative to a given voxel.
enum linkDirection {	
    X_POS=0,			//!< Positive X direction
    X_NEG=1,			//!< Negative X direction
    Y_POS=2,			//!< Positive Y direction
    Y_NEG=3,			//!< Negative Y direction
    Z_POS=4,			//!< Positive Z direction
    Z_NEG=5				//!< Negative Z direction
}; 
//! Defines each of 8 corners of a voxel.
enum voxelCorner {
    NNN = 0, //0b000
    NNP = 1, //0b001
    NPN = 2, //0b010
    NPP = 3, //0b011
    PNN = 4, //0b100
    PNP = 5, //0b101
    PPN = 6, //0b110
    PPP = 7  //0b111
}; 

typedef int voxState;
enum voxFlags { //default of each should be zero for easy clearing
    SURFACE = 1<<1, //on the surface?
    FLOOR_ENABLED = 1<<2, //interact with a floor at z=0?
    FLOOR_STATIC_FRICTION = 1<<3, //is the voxel in a state of static friction with the floor?
    COLLISIONS_ENABLED = 1<<5
};

enum linkAxis {			
    X_AXIS = 0,			//!< X Axis
    Y_AXIS = 1,			//!< Y Axis
    Z_AXIS = 2			//!< Z Axis
};

typedef unsigned char dofObject;  //Bits: 0 0 Tz, Ty, Tz, Z, Y, X. 0 if free, 1 if fixed
enum dofComponent {
	X_TRANSLATE=1<<0,
	Y_TRANSLATE=1<<1,
	Z_TRANSLATE=1<<2,
	X_ROTATE=1<<3,
	Y_ROTATE=1<<4,
	Z_ROTATE=1<<5
};

static const float HYSTERESIS_FACTOR = 1.2f; //Amount for small angle bond calculations
static const float SA_BOND_BEND_RAD = 0.05f; //Amount for small angle bond calculations
static const float SA_BOND_EXT_PERC = 0.50f; //Amount for small angle bond calculations

#endif // TYPE_H

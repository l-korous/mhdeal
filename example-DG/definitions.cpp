#include "common.h"

const ui DG_ORDER = 0;
const ui INIT_REF_NUM = 0;
const ui COMPONENT_COUNT = 5;

// boundary id
const unsigned int BOUNDARY_FRONT = 1;
const unsigned int BOUNDARY_RIGHT = 2;
const unsigned int BOUNDARY_BACK = 3;
const unsigned int BOUNDARY_LEFT = 4;
const unsigned int BOUNDARY_BOTTOM = 5;
const unsigned int BOUNDARY_TOP = 6;

const dealii::Point<DIM> p1(0., 0., 0.);
const dealii::Point<DIM> p2(.3, 1., 1.);
const dealii::Point<DIM> p3(.3, .1, 0.);
const dealii::Point<DIM> p4(1.5, 1., 1.);

const d T_FINAL = 1000.0;
const d DELTA_T = 0.001;

const bool PRINT_ALGEBRA = false;

const d GAMMA = 5.0/3.0;  // monoatomic gas in 3D
const d ETA = 1.0e-8;
// ideal gas constant
const d R = 287.14;
// specific heat capacity at constant volume
const d C_V = 717.5;
const d KAPPA = 1.0 + (R / C_V);

// Exterior pressure (dimensionless).
const double P_EXT = 2.5;
// Inlet density (dimensionless).
const double RHO_EXT = 1.0;
// Inlet x-velocity (dimensionless).
const double V1_EXT = 1.;
// Inlet y-velocity (dimensionless).
const double V2_EXT = 0.0;
// Inlet z-velocity (dimensionless).
const double V3_EXT = 0.0;
// Inlet energy.
const double E_EXT = P_EXT / (KAPPA - 1.0) + RHO_EXT * (V1_EXT*V1_EXT + V2_EXT*V2_EXT + V3_EXT*V3_EXT) / 2.0;
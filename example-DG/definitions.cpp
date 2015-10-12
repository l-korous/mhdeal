#include "definitions.h"

const ui DG_ORDER = 0;
const ui INIT_REF_NUM = 4;
const ui COMPONENT_COUNT = 5;

// boundary id
const unsigned int BOUNDARY_FRONT = 1;
const unsigned int BOUNDARY_RIGHT = 2;
const unsigned int BOUNDARY_BACK = 3;
const unsigned int BOUNDARY_LEFT = 4;
const unsigned int BOUNDARY_BOTTOM = 5;
const unsigned int BOUNDARY_TOP = 6;

bool BC_IS_IN_WEAKFORM(const unsigned int bnd_marker)
{
  if (bnd_marker == BOUNDARY_LEFT || bnd_marker == BOUNDARY_TOP)
    return false;
  else
    return true;
}

bool BC_IS_OUTFLOW(const unsigned int bnd_marker)
{
  if (bnd_marker == BOUNDARY_RIGHT)
    return true;
  else
    return false;
}

const dealii::Point<DIM> p1(0., 0., 0.);
//const dealii::Point<DIM> p2(.3, 1., 1.);
//const dealii::Point<DIM> p3(.3, .1, 0.);
const dealii::Point<DIM> p4(4.1, 1., 1.);

const d T_FINAL = 1000.0;
const d DELTA_T = 0.001;

const bool PRINT_ALGEBRA = false;
const bool DELETE_VTK = true;

const d GAMMA = 5.0 / 3.0;  // monoatomic gas in 3D
const d ETA = 1.0e-8;
// ideal gas constant
const d R = 287.14;
// specific heat capacity at constant volume
const d C_V = 717.5;
const d KAPPA = 1.0 + (R / C_V);

// Left
const double RHO_IN_LEFT = 1.0;
const double V1_IN_LEFT = 2.9;
const double V2_IN_LEFT = 0.;
const double V3_IN_LEFT = 0.;
const double P_IN_LEFT = 0.714286;

// Top
const double RHO_IN_TOP = 1.7;
const double V1_IN_TOP = 2.619334;
const double V2_IN_TOP = -0.5063;
const double V3_IN_TOP = 0.;
const double P_IN_TOP = 1.52819;

// Init
const double RHO_INIT = RHO_IN_LEFT;
const double V1_INIT = V1_IN_LEFT;
const double V2_INIT = V2_IN_LEFT;
const double V3_INIT = V3_IN_LEFT;
const double P_INIT = P_IN_LEFT;

// Inlet energy.
const double E_IN_LEFT = P_IN_LEFT / (KAPPA - 1.0) + RHO_IN_LEFT * (V1_IN_LEFT*V1_IN_LEFT + V2_IN_LEFT*V2_IN_LEFT + V3_IN_LEFT*V3_IN_LEFT) / 2.0;

const double E_IN_TOP = P_IN_TOP / (KAPPA - 1.0) + RHO_IN_TOP * (V1_IN_TOP*V1_IN_TOP + V2_IN_TOP*V2_IN_TOP + V3_IN_TOP*V3_IN_TOP) / 2.0;

const double E_INIT = P_INIT / (KAPPA - 1.0) + RHO_INIT * (V1_INIT*V1_INIT + V2_INIT*V2_INIT + V3_INIT*V3_INIT) / 2.0;

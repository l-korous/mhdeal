#include "definitions.h"

// Rad, 0 = konecne objemy
const ui DG_ORDER = 0;
// Zjemneni site
const ui INIT_REF_NUM_X = 4;
const ui INIT_REF_NUM_Y = 4;
const ui INIT_REF_NUM_Z = 4;

bool BC_INFLOW_OUTFLOW(const unsigned int bnd_marker)
{
  if (bnd_marker == BOUNDARY_LEFT || bnd_marker == BOUNDARY_TOP || bnd_marker == BOUNDARY_RIGHT)
    return true;
  else
    return false;
}

bool BC_SOLID_WALL(const unsigned int bnd_marker)
{
  if (bnd_marker == BOUNDARY_BOTTOM)
    return true;
  else
    return false;
}

bool BC_IS_SYMMETRIC(const unsigned int bnd_marker)
{
  if ((bnd_marker == BOUNDARY_RIGHT) || (bnd_marker == BOUNDARY_LEFT))
    return true;
  else
    return false;
}

const dealii::Point<DIM> p1(0., 0., 0.);
const dealii::Point<DIM> p4(4.1, 1., 1.);

const d T_FINAL = 1000.0;
const d DELTA_T = 0.003;

const bool PRINT_ALGEBRA = false;
// Delete VTK on start
const bool DELETE_VTK = true;

const d GAMMA = 5.0 / 3.0;  // monoatomic gas in 3D
const d ETA = 1.0e-8;
// ideal gas constant
const d R = 287.14;
// specific heat capacity at constant volume
const d C_V = 717.5;
// This is somehow wrong -- ?
// const d KAPPA = 1.0 + (R / C_V);

const d KAPPA = 1.4;

// External quantities
const double RHO_EXT = 1.0;
const double V1_EXT = 0.;
const double V2_EXT = 0.;
const double V3_EXT = 0.;
const double P_EXT = 0.714286;
// External energy.
const double E_EXT = P_EXT / (KAPPA - 1.0) + RHO_EXT * (V1_EXT*V1_EXT + V2_EXT*V2_EXT + V3_EXT*V3_EXT) / 2.0;
const double B1_EXT = 0.;
const double B2_EXT = 0.;
const double B3_EXT = 0.;
const double J1_EXT = 0.;
const double J2_EXT = 0.;
const double J3_EXT = 0.;
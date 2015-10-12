#include "numericalFlux.h"

using namespace dealii;

void NumFluxCentral::calculate(vec U_L, vec U_R, Point<DIM> quadPoint, Point<DIM> normal, vec& result)
{
  d x = quadPoint(0), y = quadPoint(1), z = quadPoint(2);
  d nx = normal(0), ny = normal(1), nz = normal(2);

  d r, p1, p2, p3, E;

  r = U_L[0];
  p1 = U_L[1];
  p2 = U_L[2];
  p3 = U_L[3];
  E = U_L[4];

  result[0] = f_1_1 * nx + f_2_1 * ny + f_3_1 * nz;
  result[1] = f_1_2 * nx + f_2_2 * ny + f_3_2 * nz;
  result[2] = f_1_3 * nx + f_2_3 * ny + f_3_3 * nz;
  result[3] = f_1_4 * nx + f_2_4 * ny + f_3_4 * nz;
  result[4] = f_1_5 * nx + f_2_5 * ny + f_3_5 * nz;

  r = U_R[0];
  p1 = U_R[1];
  p2 = U_R[2];
  p3 = U_R[3];
  E = U_R[4];

  result[0] += f_1_1 * nx + f_2_1 * ny + f_3_1 * nz;
  result[1] += f_1_2 * nx + f_2_2 * ny + f_3_2 * nz;
  result[2] += f_1_3 * nx + f_2_3 * ny + f_3_3 * nz;
  result[3] += f_1_4 * nx + f_2_4 * ny + f_3_4 * nz;
  result[4] += f_1_5 * nx + f_2_5 * ny + f_3_5 * nz;

  result[0] = result[0] / 2.;
  result[1] = result[1] / 2.;
  result[2] = result[2] / 2.;
  result[3] = result[3] / 2.;
  result[4] = result[4] / 2.;
}

void NumFluxUpwind::calculate(vec U_L, vec U_R, Point<DIM> quadPoint, Point<DIM> normal, vec& result)
{
  d x = quadPoint(0), y = quadPoint(1), z = quadPoint(2);
  d nx = normal(0), ny = normal(1), nz = normal(2);

  // prumerna rychlost pro rozhodnuti kam tece medium
  vec U_C(DIM);
  U_C[0] = (U_L[1] / U_L[0] + U_R[1] / U_R[0]) / 2.0;
  U_C[1] = (U_L[2] / U_L[0] + U_R[2] / U_R[0]) / 2.0;
  U_C[2] = (U_L[3] / U_L[0] + U_R[3] / U_R[0]) / 2.0;

  d v_cdot_n = U_C[0] * nx + U_C[1] * ny + U_C[2] * nz;

  d r, p1, p2, p3, E;

  if (v_cdot_n >= 0.)
  {
    r = U_L[0];
    p1 = U_L[1];
    p2 = U_L[2];
    p3 = U_L[3];
    E = U_L[4];
  }
  else
  {
    r = U_R[0];
    p1 = U_R[1];
    p2 = U_R[2];
    p3 = U_R[3];
    E = U_R[4];
  }

  result[0] = f_1_1 * nx + f_2_1 * ny + f_3_1 * nz;
  result[1] = f_1_2 * nx + f_2_2 * ny + f_3_2 * nz;
  result[2] = f_1_3 * nx + f_2_3 * ny + f_3_3 * nz;
  result[3] = f_1_4 * nx + f_2_4 * ny + f_3_4 * nz;
  result[4] = f_1_5 * nx + f_2_5 * ny + f_3_5 * nz;
}

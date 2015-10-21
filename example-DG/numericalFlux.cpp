#include "numericalFlux.h"

using namespace dealii;

static double calc_energy(double rho, double rho_v_x, double rho_v_y, double rho_v_z, double pressure, double KAPPA)
{
  double to_return = pressure / (KAPPA - 1.0) + (rho_v_x*rho_v_x + rho_v_y*rho_v_y + rho_v_z*rho_v_z) / (2.0*rho);
  if (std::abs(to_return) < 1E-12 || to_return < 0.0)
    return 1E-12;
  return to_return;
}

// Calculates pressure from other quantities.
static double calc_pressure(double rho, double rho_v_x, double rho_v_y, double rho_v_z, double energy, double KAPPA)
{
  double to_return = (KAPPA - 1.0) * (energy - (rho_v_x*rho_v_x + rho_v_y*rho_v_y + rho_v_z*rho_v_z) / (2.0*rho));
  if (std::abs(to_return) < 1E-12 || to_return < 0.0)
    return 1E-12;
  return to_return;
}

// Calculates speed of sound.
static double calc_sound_speed(double rho, double rho_v_x, double rho_v_y, double rho_v_z, double energy, double KAPPA)
{
  double to_return = std::sqrt(KAPPA * calc_pressure(rho, rho_v_x, rho_v_y, rho_v_z, energy, KAPPA) / rho);
  if (std::abs(to_return) < 1E-12 || to_return < 0.0)
    return 1E-12;
  return to_return;
}

void NumFlux::Q(double result[5], double state_vector[5], double nx, double ny, double nz)
{
  result[0] = state_vector[0];
  double temp_result_1 = nx * state_vector[1] + ny * state_vector[2] + nz * state_vector[3];
  double temp_result_2 = -ny * state_vector[1] + nx * state_vector[2];
  double temp_result_3 = -nz * state_vector[1] + nx * state_vector[3];
  result[1] = temp_result_1;
  result[2] = temp_result_2;
  result[3] = temp_result_3;
  result[4] = state_vector[4];
}

void NumFlux::Q_inv(double result[5], double state_vector[5], double nx, double ny, double nz)
{
  result[0] = state_vector[0];
  double temp_result_1 = nx * state_vector[1] - ny * state_vector[2] - nz * state_vector[3];
  double temp_result_2 = ny * state_vector[1] + nx * state_vector[2];
  double temp_result_3 = nz * state_vector[1] + nx * state_vector[3];
  result[1] = temp_result_1;
  result[2] = temp_result_2;
  result[3] = temp_result_3;
  result[4] = state_vector[4];
}

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

void NumFluxVijayasundaram::calculate(vec U_L, vec U_R, Point<DIM> quadPoint, Point<DIM> normal, vec& result)
{
  d x = quadPoint(0), y = quadPoint(1), z = quadPoint(2);
  d nx = normal(0), ny = normal(1), nz = normal(2);

  d r, p1, p2, p3, E;

  d result_arr[5];
  d U_L_arr[5];
  d U_R_arr[5];
  for (unsigned int i = 0; i < COMPONENT_COUNT; i++)
  {
    U_L_arr[i] = U_L[i];
    U_R_arr[i] = U_R[i];
  }

  d result_temp[5];
  r = U_L[0];
  p1 = U_L[1];
  p2 = U_L[2];
  p3 = U_L[3];
  E = U_L[4];

  result_temp[0] = f_1_1 * nx + f_2_1 * ny + f_3_1 * nz;
  result_temp[1] = f_1_2 * nx + f_2_2 * ny + f_3_2 * nz;
  result_temp[2] = f_1_3 * nx + f_2_3 * ny + f_3_3 * nz;
  result_temp[3] = f_1_4 * nx + f_2_4 * ny + f_3_4 * nz;
  result_temp[4] = f_1_5 * nx + f_2_5 * ny + f_3_5 * nz;

  //////////////
  r = U_R[0];
  p1 = U_R[1];
  p2 = U_R[2];
  p3 = U_R[3];
  E = U_R[4];

  result_temp[0] += f_1_1 * nx + f_2_1 * ny + f_3_1 * nz;
  result_temp[1] += f_1_2 * nx + f_2_2 * ny + f_3_2 * nz;
  result_temp[2] += f_1_3 * nx + f_2_3 * ny + f_3_3 * nz;
  result_temp[3] += f_1_4 * nx + f_2_4 * ny + f_3_4 * nz;
  result_temp[4] += f_1_5 * nx + f_2_5 * ny + f_3_5 * nz;

  for (unsigned int i = 0; i < 5; i++)
    result[i] = result_temp[i] / 2.;

  d w_mean[5];
  w_mean[0] = U_L[0] + U_R[0];
  w_mean[1] = U_L[1] + U_R[1];
  w_mean[2] = U_L[2] + U_R[2];
  w_mean[3] = U_L[3] + U_R[3];
  w_mean[4] = U_L[4] + U_R[4];

  P_plus(result_temp, w_mean, U_L_arr, nx, ny, nz);
  P_minus(result_arr, w_mean, U_R_arr, nx, ny, nz);

  for (unsigned int i = 0; i < COMPONENT_COUNT; i++)
    result[i] = result_arr[i] + result_temp[i];
}

void NumFluxVijayasundaram::P_plus(double* result, double w[5], double param[5], double nx, double ny, d nz)
{
  Q(q, w, nx, ny, nz);

  // Initialize the matrices.
  double T[5][5];
  for (unsigned int i = 0; i < 5; i++)
    for (unsigned int j = 0; j < 5; j++)
      T[i][j] = 0.0;

  double T_inv[5][5];
  for (unsigned int i = 0; i < 5; i++)
    for (unsigned int j = 0; j < 5; j++)
      T_inv[i][j] = 0.0;

  // Calculate Lambda^+.
  Lambda_plus(result);

  // Calculate the necessary rows / columns of T(T_inv).
  if (result[0] > 0) {
    T_1(T);
    T_inv_1(T_inv);
  }
  if (result[1] > 0) {
    T_2(T);
    T_inv_2(T_inv);
  }
  if (result[2] > 0) {
    T_3(T);
    T_inv_3(T_inv);
  }
  if (result[3] > 0) {
    T_4(T);
    T_inv_4(T_inv);
  }
  if (result[4] > 0) {
    T_5(T);
    T_inv_5(T_inv);
  }

  // The matrix T * Lambda * T^{-1}
  double diag_inv[5][5];
  double A_1[5][5];
  for (unsigned int i = 0; i < 5; i++)
    for (unsigned int j = 0; j < 5; j++)
      diag_inv[i][j] = result[i] * T_inv[i][j];
  for (unsigned int i = 0; i < 5; i++)
    for (unsigned int j = 0; j < 5; j++) {
      A_1[i][j] = 0;
      for (unsigned int k = 0; k < 5; k++)
        A_1[i][j] += T[i][k] * diag_inv[k][j];
    }

  // Finale.
  Q(param, param, nx, ny, nz);
  for (unsigned int i = 0; i < 5; i++) {
    result[i] = 0;
    for (unsigned int j = 0; j < 5; j++)
      result[i] += A_1[i][j] * param[j];
  }
  Q_inv(result, result, nx, ny, nz);
}

void NumFluxVijayasundaram::P_minus(double* result, double w[5], double param[5], double nx, double ny, d nz)
{
  Q(q, w, nx, ny, nz);

  // Initialize the matrices.
  double T[5][5];
  for (unsigned int i = 0; i < 5; i++)
    for (unsigned int j = 0; j < 5; j++)
      T[i][j] = 0.0;

  double T_inv[5][5];
  for (unsigned int i = 0; i < 5; i++)
    for (unsigned int j = 0; j < 5; j++)
      T_inv[i][j] = 0.0;

  // Calculate Lambda^-.
  Lambda_minus(result);

  // Calculate the necessary rows / columns of T(T_inv).
  if (result[0] < 0) {
    T_1(T);
    T_inv_1(T_inv);
  }
  if (result[1] < 0) {
    T_2(T);
    T_inv_2(T_inv);
  }
  if (result[2] < 0) {
    T_3(T);
    T_inv_3(T_inv);
  }
  if (result[3] < 0) {
    T_4(T);
    T_inv_4(T_inv);
  }
  if (result[4] < 0) {
    T_5(T);
    T_inv_5(T_inv);
  }

  // The matrix T * Lambda * T^{-1}
  double diag_inv[5][5];
  double A_1[5][5];
  for (unsigned int i = 0; i < 5; i++)
    for (unsigned int j = 0; j < 5; j++)
      diag_inv[i][j] = result[i] * T_inv[i][j];
  for (unsigned int i = 0; i < 5; i++)
    for (unsigned int j = 0; j < 5; j++) {
      A_1[i][j] = 0;
      for (unsigned int k = 0; k < 5; k++)
        A_1[i][j] += T[i][k] * diag_inv[k][j];
    }

  // Finale.
  Q(param, param, nx, ny, nz);
  for (unsigned int i = 0; i < 5; i++) {
    result[i] = 0;
    for (unsigned int j = 0; j < 5; j++)
      result[i] += A_1[i][j] * param[j];
  }
  Q_inv(result, result, nx, ny, nz);
}

void NumFluxVijayasundaram::Lambda_plus(double result[5])
{
  a = calc_sound_speed(q[0], q[1], q[2], q[3], q[4], KAPPA);
  u = q[1] / q[0];
  v = q[2] / q[0];
  w = q[3] / q[0];
  V = u*u + v*v + w*w;
  result[0] = u - a < 0 ? 0 : u - a;
  result[1] = u < 0 ? 0 : u;
  result[2] = u < 0 ? 0 : u;
  result[3] = u < 0 ? 0 : u;
  result[4] = u + a < 0 ? 0 : u + a;
}

void NumFluxVijayasundaram::Lambda_minus(double result[5])
{
  a = calc_sound_speed(q[0], q[1], q[2], q[3], q[4], KAPPA);
  u = q[1] / q[0];
  v = q[2] / q[0];
  w = q[3] / q[0];
  V = u*u + v*v + w*w;
  result[0] = u - a < 0 ? u - a : 0;
  result[1] = u < 0 ? u : 0;
  result[2] = u < 0 ? u : 0;
  result[3] = u < 0 ? u : 0;
  result[4] = u + a < 0 ? u + a : 0;
}

void NumFluxVijayasundaram::Lambda(double result[5])
{
  a = calc_sound_speed(q[0], q[1], q[2], q[3], q[4], KAPPA);
  u = q[1] / q[0];
  v = q[2] / q[0];
  w = q[3] / q[0];
  V = u*u + v*v + w*w;
  result[0] = u - a;
  result[1] = u;
  result[2] = u;
  result[3] = u;
  result[4] = u + a;
}

void NumFluxVijayasundaram::T_1(double result[5][5])
{
  result[0][0] = 1.0;
  result[1][0] = u - a;
  result[2][0] = v;
  result[3][0] = w;
  result[4][0] = (V / 2.0) + (a*a / (KAPPA - 1.0)) - (u * a);
}
void NumFluxVijayasundaram::T_2(double result[5][5])
{
  result[0][1] = 1.0;
  result[1][1] = u;
  result[2][1] = v;
  result[3][1] = w;
  result[4][1] = V / 2.0;
}
void NumFluxVijayasundaram::T_3(double result[5][5])
{
  result[0][2] = 1.0;
  result[1][2] = u;
  result[2][2] = v - a;
  result[3][2] = w;
  result[4][2] = (V / 2.0) - v * a;
}
void NumFluxVijayasundaram::T_4(double result[5][5])
{
  result[0][3] = 1.0;
  result[1][3] = u;
  result[2][3] = v;
  result[3][3] = w - a;
  result[4][3] = (V / 2.0) - w * a;
}
void NumFluxVijayasundaram::T_5(double result[5][5])
{
  result[0][3] = 1.0;
  result[1][3] = u + a;
  result[2][3] = v;
  result[3][3] = w;
  result[4][3] = (V / 2.0) + (a * a / (KAPPA - 1.0)) + (u * a);
}

void NumFluxVijayasundaram::T_inv_1(double result[5][5])
{
  result[0][0] = (1.0 / (a * a)) * (0.5 * (((KAPPA - 1) * V / 2.0) + u * a));
  result[0][1] = (1.0 / (a * a)) * (-(a + u * (KAPPA - 1.0)) / 2.0);
  result[0][2] = (1.0 / (a * a)) * (-(v * (KAPPA - 1.0)) / 2.0);
  result[0][3] = (1.0 / (a * a)) * (-(w * (KAPPA - 1.0)) / 2.0);
  result[0][4] = (1.0 / (a * a)) * (KAPPA - 1.0) / 2.0;
}
void NumFluxVijayasundaram::T_inv_2(double result[5][5])
{
  result[1][0] = (1.0 / (a * a)) * (a * a - (v + w) * a - (KAPPA - 1.0) * (V / 2.0));
  result[1][1] = (1.0 / (a * a)) * u * (KAPPA - 1.0);
  result[1][2] = (1.0 / (a * a)) * (a + v * (KAPPA - 1.0));
  result[1][3] = (1.0 / (a * a)) * (a + w * (KAPPA - 1.0));
  result[1][4] = (1.0 / (a * a)) * (1.0 - KAPPA);
}
void NumFluxVijayasundaram::T_inv_3(double result[5][5])
{
  result[2][0] = (1.0 / (a * a)) * v * a;
  result[2][1] = (1.0 / (a * a)) * 0.0;
  result[2][2] = (1.0 / (a * a)) * (-a);
  result[2][3] = (1.0 / (a * a)) * 0.0;
  result[2][4] = (1.0 / (a * a)) * 0.0;
}
void NumFluxVijayasundaram::T_inv_4(double result[5][5])
{
  result[3][0] = (1.0 / (a * a)) * w * a;
  result[3][1] = (1.0 / (a * a)) * 0.0;
  result[3][2] = (1.0 / (a * a)) * 0.;
  result[3][3] = (1.0 / (a * a)) * (-a);
  result[3][4] = (1.0 / (a * a)) * 0.0;
}
void NumFluxVijayasundaram::T_inv_5(double result[5][5])
{
  result[4][0] = (1.0 / (a * a)) * (0.5 * (((KAPPA - 1.0) * V / 2.0) - u * a));
  result[4][1] = (1.0 / (a * a)) * (a - u * (KAPPA - 1.0)) / 2.0;
  result[4][2] = (1.0 / (a * a)) * (-(v * (KAPPA - 1.0)) / 2.0);
  result[4][3] = (1.0 / (a * a)) * (-(w * (KAPPA - 1.0)) / 2.0);
  result[4][4] = (1.0 / (a * a)) * (KAPPA - 1.0) / 2.0;
}



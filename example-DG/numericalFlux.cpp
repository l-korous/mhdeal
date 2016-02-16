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

template<typename arr>
void NumFlux::Q(arr result, arr state_vector, double nx, double ny, double nz)
{
  // Density
  result[0] = state_vector[0];
  // Momentum (x3)
  double temp_result_1 = nx * state_vector[1] + ny * state_vector[2] + nz * state_vector[3];
  double temp_result_2 = -ny * state_vector[1] + nx * state_vector[2];
  double temp_result_3 = -nz * state_vector[1] + nx * state_vector[3];
  result[1] = temp_result_1;
  result[2] = temp_result_2;
  result[3] = temp_result_3;
  // Energy
  result[4] = state_vector[4];

  // B (x3)
  double temp_result_5 = nx * state_vector[5] + ny * state_vector[6] + nz * state_vector[7];
  double temp_result_6 = -ny * state_vector[5] + nx * state_vector[6];
  double temp_result_7 = -nz * state_vector[5] + nx * state_vector[7];
  result[5] = temp_result_5;
  result[6] = temp_result_6;
  result[7] = temp_result_7;
}
template void NumFlux::Q(double[COMPONENT_COUNT_T],double[COMPONENT_COUNT_T],double,double,double);
template void NumFlux::Q(vec,vec,double,double,double);

template<typename arr>
void NumFlux::Q_inv(arr result, arr state_vector, double nx, double ny, double nz)
{
  // Density
  result[0] = state_vector[0];
  // Momentum (x3)
  double temp_result_1 = nx * state_vector[1] - ny * state_vector[2] - nz * state_vector[3];
  double temp_result_2 = ny * state_vector[1] + nx * state_vector[2];
  double temp_result_3 = nz * state_vector[1] + nx * state_vector[3];
  result[1] = temp_result_1;
  result[2] = temp_result_2;
  result[3] = temp_result_3;
  // Energy
  result[4] = state_vector[4];

  // B (x3)
  double temp_result_5 = nx * state_vector[5] - ny * state_vector[6] - nz * state_vector[7];
  double temp_result_6 = ny * state_vector[5] + nx * state_vector[6];
  double temp_result_7 = nz * state_vector[5] + nx * state_vector[7];
  result[4] = temp_result_5;
  result[5] = temp_result_6;
  result[6] = temp_result_7;
}
template void NumFlux::Q_inv(double[COMPONENT_COUNT_T],double[COMPONENT_COUNT_T],double,double,double);
template void NumFlux::Q_inv(vec,vec,double,double,double);

#pragma region central_and_upwind
void NumFluxCentral::calculate(vec U_L, vec U_R, Point<DIM> quadPoint, Point<DIM> normal, vec& result, ui only_part)
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

void NumFluxUpwind::calculate(vec U_L, vec U_R, Point<DIM> quadPoint, Point<DIM> normal, vec& result, ui only_part)
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
#pragma endregion
#pragma region vijayasundaram
void NumFluxVijayasundaram::calculate(vec U_L, vec U_R, Point<DIM> quadPoint, Point<DIM> normal, vec& result, ui only_part)
{
  d nx = normal(0), ny = normal(1), nz = normal(2);

  d U_L_arr[COMPONENT_COUNT];
  d U_R_arr[COMPONENT_COUNT];
  d w_mean[COMPONENT_COUNT];

  d result_arr[COMPONENT_COUNT];
  d result_temp[COMPONENT_COUNT];

  for (ui i = 0; i < COMPONENT_COUNT; i++)
  {
    w_mean[i] = (U_L[i] + U_R[i]) / 2.;
    U_L_arr[i] = U_L[i];
    U_R_arr[i] = U_R[i];
    result_arr[i] = 0.;
    result_temp[i] = 0.;
  }

  if (only_part == 0 || only_part == 1)
    P_plus(result_temp, w_mean, U_L_arr, nx, ny, nz);
  if (only_part == 0 || only_part == 2)
    P_minus(result_arr, w_mean, U_R_arr, nx, ny, nz);

  for (unsigned int i = 0; i < COMPONENT_COUNT; i++)
    result[i] = result_arr[i] + result_temp[i];
}

d NumFluxVijayasundaram::calculate(vec U_L_prev, vec U_R_prev, dealii::Point<DIM> quadPoint, dealii::Point<DIM> normal, ui comp_i, ui comp_j, ui only_part)
{
  d nx = normal(0), ny = normal(1), nz = normal(2);
  d e[5][5] =
  {
    { 1., 0, 0, 0, 0 },
    { 0, 1., 0, 0, 0 },
    { 0, 0, 1., 0, 0 },
    { 0, 0, 0, 1., 0 },
    { 0, 0, 0, 0, 1. }
  };

  d U_L[COMPONENT_COUNT], U_R[COMPONENT_COUNT];
  d result[COMPONENT_COUNT];
  for (ui i = 0; i < 5; i++)
  {
    U_L[i] = U_L_prev[i];
    U_R[i] = U_R_prev[i];
    result[i] = 0.;
  }


  if (only_part == 1)
    P_plus(result, U_L, e[comp_j], nx, ny, nz);
  if (only_part == 2)
    P_minus(result, U_R, e[comp_j], nx, ny, nz);

  return result[comp_i];
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

void NumFluxVijayasundaram::T(double result[5][5])
{
  this->T_1(result);
  this->T_2(result);
  this->T_3(result);
  this->T_4(result);
  this->T_5(result);
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
  result[0][4] = 1.0;
  result[1][4] = u + a;
  result[2][4] = v;
  result[3][4] = w;
  result[4][4] = (V / 2.0) + (a * a / (KAPPA - 1.0)) + (u * a);
}

void NumFluxVijayasundaram::T_inv(double result[5][5])
{
  this->T_inv_1(result);
  this->T_inv_2(result);
  this->T_inv_3(result);
  this->T_inv_4(result);
  this->T_inv_5(result);
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

#pragma endregion

#define _Pb_  // The magnetic pressure is calculated separately for left and right flux, otherwise the average of B is used
#define SMNUM 1e-8
// HLLD approximate Riemann solver (Mioyshi T., Kusano K., 2005)

void NumFluxHLLD::calculate(vec ul, vec ur, Point<DIM> /*quadPoint*/,
                            Point<DIM> normal, vec& F, ui /*only_part*/)
//void /*class::*/hlld(double *ul,double *ur, double *F)
{
    double srdl,srdr,Fl[COMPONENT_COUNT_T],Fr[COMPONENT_COUNT_T],hl[2],hr[2],Uk,Um,E2,E3,Sl,Sr,pml,pmr,B,B2;
    double Udl[COMPONENT_COUNT_T],Udr[COMPONENT_COUNT_T],Ul[COMPONENT_COUNT_T],Ur[COMPONENT_COUNT_T],cl,cm,cr,ptl,ptr;
    double sp[5],sml,smr,ptst,ptstr,vbstl,vbstr,Bsgnl,Bsgnr,invsumd;

    Q(ul,ul,normal[0],normal[1],normal[2]);
    Q(ur,ur,normal[0],normal[1],normal[2]);
    
    B=0.5*(ul[4]+ur[4]);  // Simple average of mag. field in direction of normal vector
    B2=B*B;

    // Calculate left flux
    hl[0] = 1.0/ul[0];
    Uk=0.5*hl[0]*(ul[1]*ul[1]+ul[2]*ul[2]+ul[3]*ul[3]);
    Um=0.5*(ul[4]*ul[4]+ul[5]*ul[5]+ul[6]*ul[6]);       // ifdef _Pb_ then use B2 instead of ul[4]^2
    hl[1]=(GAMMA-1)*(ul[7]-Uk-Um);
    //E1=hl[0]*(ul[3]*ul[5]-ul[2]*ul[6]);
    E2=hl[0]*(ul[1]*ul[6]-ul[3]*ul[4]);
    E3=hl[0]*(ul[2]*ul[4]-ul[1]*ul[5]);
    
    Fl[0]=ul[1];
    Fl[1]=hl[0]*ul[1]*ul[1]-ul[4]*ul[4]+Um+hl[1];
    Fl[2]=hl[0]*ul[1]*ul[2]-ul[4]*ul[5];
    Fl[3]=hl[0]*ul[1]*ul[3]-ul[4]*ul[6];
    Fl[4]=0.0;
    Fl[5]=-E3;
    Fl[6]= E2;
    Fl[7]=hl[0]*ul[1]*(hl[1]*GAMMA/(GAMMA-1.0)+Uk)+(E2*ul[6]-E3*ul[5]);
   
    // Calculate right flux
    hr[0] = 1.0/ur[0];
    Uk=0.5*hr[0]*(ur[1]*ur[1]+ur[2]*ur[2]+ur[3]*ur[3]);
    Um=0.5*(ur[4]*ur[4]+ur[5]*ur[5]+ur[6]*ur[6]);       // ifdef _Pb_ then use B2 instead of ur[4]^2
    hr[1]=(GAMMA-1)*(ur[7]-Uk-Um);
    //E1=hr[0]*(ur[3]*ur[5]-ur[2]*ur[6]);
    E2=hr[0]*(ur[1]*ur[6]-ur[3]*ur[4]);
    E3=hr[0]*(ur[2]*ur[4]-ur[1]*ur[5]);
    
    Fr[0]=ur[1];
    Fr[1]=hr[0]*ur[1]*ur[1]-ur[4]*ur[4]+Um+hr[1];
    Fr[2]=hr[0]*ur[1]*ur[2]-ur[4]*ur[5];
    Fr[3]=hr[0]*ur[1]*ur[3]-ur[4]*ur[6];
    Fr[4]=0.0;
    Fr[5]=-E3;
    Fr[6]= E2;
    Fr[7]=hr[0]*ur[1]*(hr[1]*GAMMA/(GAMMA-1.0)+Uk)+(E2*ur[6]-E3*ur[5]);
    
#ifdef _Pb_
    pml=0.5*(ul[4]*ul[4]+ul[5]*ul[5]+ul[6]*ul[6]);
    pmr=0.5*(ur[4]*ur[4]+ur[5]*ur[5]+ur[6]*ur[6]);
    // fast magnetoacoustic speed
    cl=GAMMA*hl[1]+2.0*pml;
    cl=sqrt(0.5*hl[0]*(cl+sqrt(cl*cl-4.0*GAMMA*hl[1]*ul[4]*ul[4])));
    cr=GAMMA*hr[1]+2.0*pmr;
    cr=sqrt(0.5*hr[0]*(cr+sqrt(cr*cr-4.0*GAMMA*hr[1]*ur[4]*ur[4])));
    
    ptl=hl[1]+pml;  // total pressure
    ptr=hr[1]+pmr;
#else
    pml=0.5*(B2+ul[5]*ul[5]+ul[6]*ul[6]);
    pmr=0.5*(B2+ur[5]*ur[5]+ur[6]*ur[6]);
      // fast magnetoacoustic speed
    cl=GAMMA*hl[1]+(GAMMA+2.0)*pml;
    cl=sqrt(0.5*hl[0]*(cl+sqrt(cl*cl-4.0*GAMMA*hl[1]*B2)));
    cr=GAMMA*hr[1]+(GAMMA+2.0)*pmr;
    cr=sqrt(0.5*hr[0]*(cr+sqrt(cr*cr-4.0*GAMMA*hr[1]*B2)));
    
    ptl=hl[1]+GAMMA*pml;  // total pressure
    ptr=hr[1]+GAMMA*pmr;
#endif
 
    // maximum of fast magnetoacoustic speeds L/R
    cm=(cl>cr)?cl:cr;
    if (ul[1]*hl[0]<=ur[1]*hr[0]){
      sp[0]=ul[1]*hl[0]-cm;
      sp[4]=ur[1]*hr[0]+cm;
    }else{
      sp[0]=ur[1]*hr[0]-cm;
      sp[4]=ul[1]*hl[0]+cm;
    }
    
  //   cm=sqrt((GAMMA-1.0)/(2.0*GAMMA));
  //   if (sp[0]>ul[1]-cl*cm) sp[0]=(ul[1]-cl*cm);
  //   if (sp[4]<ur[1]+cr*cm) sp[4]=(ur[1]+cr*cm);

    // Upwind flux in the case of supersonic flow
    if (sp[0]>=0.0){  // use F_L
      for(register int j=0;j<COMPONENT_COUNT_T;j++)
        F[j]=Fl[j];
#ifndef _Pb_
      F[1]+=GAMMA*pml-B2;
      F[7]-=GAMMA*ul[1]*hl[0]*pml;
#endif
      Q_inv(F,F,normal[0],normal[1],normal[2]);
      return;
    }
    if (sp[4]<=0.0){  // use F_R
      for(register int j=0;j<COMPONENT_COUNT_T;j++)
        F[j]=Fr[j];
#ifndef _Pb_
      F[1]+=GAMMA*pmr-B2;
      F[7]-=GAMMA*ur[1]*hr[0]*pmr;
#endif
      Q_inv(F,F,normal[0],normal[1],normal[2]);
      return;
    }

    // Determine Alfven and middle speeds
    Sl=sp[0]-ul[1]*hl[0];
    Sr=sp[4]-ur[1]*hr[0];
    sp[2]=( ul[1]*Sl-ur[1]*Sr-ptl+ptr
          )/(ul[0]*Sl-ur[0]*Sr);
    sml=sp[0]-sp[2];
    smr=sp[4]-sp[2];
    
    Ul[0]=ul[0]*Sl/sml;  // Density
    Ur[0]=ur[0]*Sr/smr;
  
#ifdef _Pb_
    Ul[4]=Udl[4]=ul[4];
    Ur[4]=Udr[4]=ur[4];
#else
    Ul[4]=Ur[4]=Udl[4]=Udr[4]=B; // Magnetic field Bx (normal direction)
#endif
    srdl=sqrt(Ul[0]);
    srdr=sqrt(Ur[0]);
    
    sp[1]=sp[2]-fabs(Ul[4])/srdl;  // Sl*
    sp[3]=sp[2]+fabs(Ur[4])/srdr;  // Sr*
  
    ptst = ptl+ ul[0]*Sl*(Sl-sml);
    ptstr = ptr+ ur[0]*Sr*(Sr-smr);
    
    // F*_L
    Ul[1]=Ul[0]*sp[2];
    
    cl=ul[0]*Sl*sml-Ul[4]*Ul[4];
    if (fabs(cl)<SMNUM*ptst){
      Ul[2]=Ul[0]*ul[2]*hl[0];
      Ul[3]=Ul[0]*ul[3]*hl[0];
      
      Ul[5]=ul[5];
      Ul[6]=ul[6];
    }else{
      cl=1.0/cl;
      cm=Ul[4]*(Sl-sml)*cl;
      Ul[2]=Ul[0]*(ul[2]*hl[0]-ul[5]*cm);
      Ul[3]=Ul[0]*(ul[3]*hl[0]-ul[6]*cm);
      cm=(ul[0]*Sl*Sl-Ul[4]*Ul[4])*cl;
      Ul[5]=ul[5]*cm;
      Ul[6]=ul[6]*cm;
    }
    vbstl=(Ul[1]*Ul[4]+Ul[2]*Ul[5]+Ul[3]*Ul[6])/Ul[0];
    
    Ul[7]=(Sl*ul[7]-ptl*ul[1]*hl[0]+ptst*sp[2] + Ul[4]*
            ((ul[1]*ul[4]+ul[2]*ul[5]+ul[3]*ul[6])*hl[0]-vbstl))/sml;
    
    // F*_R
    Ur[1]=Ur[0]*sp[2];
    cl=ur[0]*Sr*smr-Ur[4]*Ur[4];
    if (fabs(cl)<SMNUM*ptstr){
      Ur[2]=Ur[0]*ur[2]*hr[0];
      Ur[3]=Ur[0]*ur[3]*hr[0];
      
      Ur[5]=ur[5];
      Ur[6]=ur[6];
    }else{
      cl=1.0/cl;
      cm=Ur[4]*(Sr-smr)*cl;
      Ur[2]=Ur[0]*(ur[2]*hr[0]-ur[5]*cm);
      Ur[3]=Ur[0]*(ur[3]*hr[0]-ur[6]*cm);
      cm=(ur[0]*Sr*Sr-Ur[4]*Ur[4])*cl;
      Ur[5]=ur[5]*cm;
      Ur[6]=ur[6]*cm;
    }
    vbstr=(Ur[1]*Ur[4]+Ur[2]*Ur[5]+Ur[3]*Ur[6])/Ur[0];
    
    Ur[7]=(Sr*ur[7]-ptr*ur[1]*hr[0]+ ptstr*sp[2]+Ur[4]*
          ((ur[1]*Ur[4]+ur[2]*ur[5]+ur[3]*ur[6])*hr[0]-vbstr))/smr;

    if (sp[1]>=0.0){
      for(register int j=0;j<COMPONENT_COUNT_T;j++)
        F[j]=Fl[j]+sp[0]*(Ul[j]-ul[j]);
#ifndef _Pb_
      F[1]+=GAMMA*pml-B2;
      F[4]=0.0;
      F[7]-=GAMMA*ul[1]*hl[0]*pml;
#endif
      Q_inv(F,F,normal[0],normal[1],normal[2]);
      return;
    }
    if (sp[3]<=0.0 && sp[2]<0.0){
      for(register int j=0;j<COMPONENT_COUNT_T;j++)
        F[j]=Fr[j]+sp[4]*(Ur[j]-ur[j]);
#ifndef _Pb_
      F[1]+=GAMMA*pmr-B2;
      F[4]=0.0;
      F[7]-=GAMMA*ur[1]*hr[0]*pmr;
#endif
      Q_inv(F,F,normal[0],normal[1],normal[2]);
      return;
    }
    
    // F**_L and F**_R
    if (B2<SMNUM*(ptst+ptstr)){
      for(register int j=0;j<COMPONENT_COUNT_T;j++){
        Udl[j]=Ul[j];
        Udr[j]=Ur[j];
      }
    }else{
      invsumd=1.0/(srdl+srdr);
      Bsgnl = (Ul[4]>0.0)?1.0:-1.0;
      Bsgnr = (Ur[4]>0.0)?1.0:-1.0;
      
      Udl[0]=Ul[0];
      Udr[0]=Ur[0];
      
      Udl[1]=Ul[1];
      Udr[1]=Ur[1];
      
      cm=invsumd*(srdl*Ul[2]/Ul[0]+srdr*Ur[2]/Ur[0]);
      cl=invsumd*(Ur[5]-Ul[5]);
      Udl[2]=Ul[0]*(cm+Bsgnl*cl);
      Udr[2]=Ur[0]*(cm+Bsgnr*cl);
      
      cm=invsumd*(srdl*Ul[3]/Ul[0]+srdr*Ur[3]/Ur[0]);
      cl=invsumd*(Ur[6]-Ul[6]);
      Udl[3]=Ul[0]*(cm+Bsgnl*cl);
      Udr[3]=Ur[0]*(cm+Bsgnr*cl);

      cm=invsumd*(srdl*Ur[5]+srdr*Ul[5]);
      cl=invsumd*srdl*srdr*(Ur[2]/Ur[0]-Ul[2]/Ul[0]);
      Udl[5]=cm+Bsgnl*cl;
      Udr[5]=cm+Bsgnr*cl;
      
      cm=invsumd*(srdl*Ur[6]+srdr*Ul[6]);
      cl=invsumd*srdl*srdr*(Ur[3]/Ur[0]-Ul[3]/Ul[0]);
      Udl[6]=cm+Bsgnl*cl;
      Udr[6]=cm+Bsgnr*cl;
      
      Udl[7]=Ul[7]-srdl*Bsgnl*(vbstl-sp[2]*Ul[4]-(Udl[2]*Udl[5]+Udl[3]*Udl[6])/Udl[0]);
      Udr[7]=Ur[7]+srdr*Bsgnr*(vbstr-sp[2]*Ur[4]-(Udr[2]*Udr[5]+Udr[3]*Udr[6])/Udr[0]);
    }

    if (sp[2]>=0.0){
      cm=sp[1]-sp[0];
      for(register int j=0;j<COMPONENT_COUNT_T;j++)
        F[j]=Fl[j]+sp[1]*Udl[j]-sp[0]*ul[j]-cm*Ul[j];
#ifndef _Pb_
      F[1]+=GAMMA*pml-B2;
      F[4]=0.0;
      F[7]-=GAMMA*ul[1]*hl[0]*pml;
#endif
      Q_inv(F,F,normal[0],normal[1],normal[2]);
      return;
    }
    cm=sp[3]-sp[4];
    for(register int j=0;j<COMPONENT_COUNT_T;j++){
      F[j]=Fr[j]+sp[3]*Udr[j]-sp[4]*ur[j]-cm*Ur[j];
    }
#ifndef _Pb_
    F[1]+=GAMMA*pmr-B2;
    F[4]=0.0;
    F[7]-=GAMMA*ur[1]*hr[0]*pmr;
#endif
    Q_inv(F,F,normal[0],normal[1],normal[2]);
    return;
}

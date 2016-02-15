#ifndef NUM_FLUX_H
#define NUM_FLUX_H

#include "definitions.h"

class NumFlux
{
public:
  // only_part: utility parameter, used in Vijayasundaram as:
  // - 1: only take into account the positive eigenvalues
  // - 2: only the negative ones
  virtual void calculate(vec U_L, vec U_R, dealii::Point<DIM> quadPoint, dealii::Point<DIM> normal, vec& result, ui only_part = 0) = 0;
  virtual d calculate(vec U_L_prev, vec U_R_prev, dealii::Point<DIM> quadPoint, dealii::Point<DIM> normal, ui comp_i, ui comp_j, ui only_part = 0) = 0;
  
  /// Rotates the state_vector into the local coordinate system.
  void Q(double result[5], double state_vector[5], double nx, double ny, double nz);

  /// Rotates the state_vector back from the local coordinate system.
  void Q_inv(double result[5], double state_vector[5], double nx, double ny, double nz);

};

class NumFluxCentral : public NumFlux
{
public:
  void calculate(vec U_L, vec U_R, dealii::Point<DIM> quadPoint, dealii::Point<DIM> normal, vec& result, ui only_part = 0);
};

class NumFluxUpwind : public NumFlux
{
public:
  void calculate(vec U_L, vec U_R, dealii::Point<DIM> quadPoint, dealii::Point<DIM> normal, vec& result, ui only_part = 0);
};

class NumFluxVijayasundaram : public NumFlux
{
public:
  void calculate(vec U_L, vec U_R, dealii::Point<DIM> quadPoint, dealii::Point<DIM> normal, vec& result, ui only_part);
  d calculate(vec U_L_prev, vec U_R_prev, dealii::Point<DIM> quadPoint, dealii::Point<DIM> normal, ui comp_i, ui comp_j, ui only_part);

  void P_plus(double* result, double w[5], double param[5], double nx, double ny, d nz);

  void P_minus(double* result, double w[5], double param[5], double nx, double ny, d nz);

  // Also calculates the speed of sound.
  void Lambda_plus(double result[5]);

  // Also calculates the speed of sound.
  void Lambda_minus(double result[5]);

  // Calculates all eigenvalues.
  void Lambda(double result[5]);

  void T(double result[5][5]);
  void T_1(double result[5][5]);
  void T_2(double result[5][5]);
  void T_3(double result[5][5]);
  void T_4(double result[5][5]);
  void T_5(double result[5][5]);

  void T_inv(double result[5][5]);
  void T_inv_1(double result[5][5]);
  void T_inv_2(double result[5][5]);
  void T_inv_3(double result[5][5]);
  void T_inv_4(double result[5][5]);
  void T_inv_5(double result[5][5]);

  double q[5];
private:

  // x-velocity, y-velocity, z-velocity, magnitude.
  double u, v, w, V;

  // Speeds of sound.
  double a;
  double a_L;
  double a_R;
  double a_L_star;
  double a_R_star;
  double a_1;
  double a_3;
  double a_B; // Boundary.
};



#endif

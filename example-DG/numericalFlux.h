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
  
  /// Rotates the state_vector into the local coordinate system.
  template<typename arr>
  void Q(arr result, arr state_vector, double nx, double ny, double nz);

  /// Rotates the state_vector back from the local coordinate system.
  template<typename arr>
  void Q_inv(arr result, arr state_vector, double nx, double ny, double nz);
};

class NumFluxCentral : public NumFlux
{
public:
  void calculate(vec U_L, vec U_R, dealii::Point<DIM> quadPoint, dealii::Point<DIM> normal, vec& result, ui only_part = 0);
};

class NumFluxHLLD : public NumFlux
{
public:
  void calculate(vec U_L, vec U_R, dealii::Point<DIM> quadPoint, dealii::Point<DIM> normal, vec& result, ui only_part = 0);
};

class NumFluxUpwind : public NumFlux
{
public:
  void calculate(vec U_L, vec U_R, dealii::Point<DIM> quadPoint, dealii::Point<DIM> normal, vec& result, ui only_part = 0);
};

#endif

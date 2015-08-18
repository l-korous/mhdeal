#ifndef NUMERICAL_FLUX_H
#define NUMERICAL_FLUX_H

#include "definitions.h"

class NumFlux
{
public:
  virtual void calculate(vec U_L, vec U_R, Point<DIM> normal, vec result) = 0;
};

class NumFluxLaxFriedrichs
{
public:
  void calculate(vec U_L, vec U_R, Point<DIM> normal, vec result);
};

#endif // NUMERICAL_FLUX_H

#ifndef _BOUNDARY_CONDITION_TD_H
#define _BOUNDARY_CONDITION_TD_H

#include "util.h"
#include "parameters.h"
#include "boundaryCondition.h"
#include "equations.h"

template <int dim>
class BoundaryConditionTitovDemoulin : public BoundaryCondition<EquationsTypeMhd, dim>
{
public:
  BoundaryConditionTitovDemoulin(Parameters<dim>&);

  void bc_vector_value(int boundary_no, const Point<dim> &point, InputVector &result, const InputVector &W_plus) const;

private:
  double beta;
  double Lg;
  double invLg;
  double N_t;
  double R_t;
  double d2R_t;
  double L2R_t;
  double q_mag;
  double iSgn;
  double heliFactor;
  double Tc2Tp;
  double t_rho;
  double densGrad;
};

#endif

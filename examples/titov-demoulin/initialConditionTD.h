#ifndef _INITIAL_CONDITION_OT_H
#define _INITIAL_CONDITION_OT_H

#include "util.h"
#include "parameters.h"
#include "initialCondition.h"
#include "equations.h"

template <EquationsType equationsType, int dim>
class InitialConditionTitovDemoulin : public InitialCondition<equationsType, dim>
{
public:
  InitialConditionTitovDemoulin(Parameters<dim>&);
  void vector_value(const std::vector<Point<dim> >&, std::vector<std::array<double, Equations<equationsType, dim>::n_components> >&) const;

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

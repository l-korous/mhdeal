#ifndef _INITIAL_CONDITION_TD_H
#define _INITIAL_CONDITION_TD_H

#include "util.h"
#include "parameters.h"
#include "initialCondition.h"
#include "equations.h"
#include "parametersTD.h"

template <EquationsType equationsType, int dim>
class InitialConditionTitovDemoulin : public InitialCondition<equationsType, dim>
{
public:
  InitialConditionTitovDemoulin(Parameters<dim>&, TitovDemoulinParameters&);
  void vector_value(const std::vector<Point<dim> >&, std::vector<std::array<double, Equations<equationsType, dim>::n_components> >&) const;

private:
  TitovDemoulinParameters& td_parameters;
  double invL_G;
  double iSgn;
  double densGrad;
  double d2R;
  double L2R;
  double R2L;
  double q_mag;
};

#endif

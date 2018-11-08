#ifndef _INITIAL_CONDITION_CS_H
#define _INITIAL_CONDITION_CS_H

#include "util.h"
#include "parameters.h"
#include "initialCondition.h"
#include "equations.h"
#include "parametersCS.h"

template <EquationsType equationsType, int dim>
class InitialConditionCS : public InitialCondition<equationsType, dim>
{
public:
  InitialConditionCS(Parameters<dim>&, CSParameters&);
  void vector_value(const std::vector<Point<dim> >&, std::vector<std::array<double, Equations<equationsType, dim>::n_components> >&) const;

private:
  CSParameters& cs_parameters;
  double invL_G;
  double d2R;
  double L2R;
  double R2L;
  double q_mag;

};

#endif

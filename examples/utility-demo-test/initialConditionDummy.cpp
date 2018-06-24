#include "completeEllipticIntegrals.h"
#include "initialConditionDummy.h"
#include "equationsMhd.h"

template <EquationsType equationsType, int dim>
InitialConditionDummy<equationsType, dim>::InitialConditionDummy(Parameters<dim>& parameters) :
  InitialCondition<equationsType, dim>(parameters)
{
}


template <EquationsType equationsType, int dim>
void InitialConditionDummy<equationsType, dim>::vector_value(const std::vector<Point<dim> > &points,
  std::vector<std::array<double, Equations<equationsType, dim>::n_components> >&result) const
{
  for (unsigned int i = 0; i < points.size(); ++i)
  {
    result[i][0] = 1.;
    result[i][1] = 0.;
    result[i][2] = 0.;
    result[i][3] = 0.;
    result[i][5] = 0.;
    result[i][6] = 0.;
    result[i][7] = 0.;
    result[i][4] = 1.0;
  }
}

template class InitialConditionDummy<EquationsTypeMhd, 3>;

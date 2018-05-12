#include "completeEllipticIntegrals.h"
#include "initialConditionMhdBlast.h"
#include "equationsMhd.h"

template <EquationsType equationsType, int dim>
InitialConditionMhdBlast<equationsType, dim>::InitialConditionMhdBlast(Parameters<dim>& parameters) :
  InitialCondition<equationsType, dim>(parameters)
{
}


template <EquationsType equationsType, int dim>
void InitialConditionMhdBlast<equationsType, dim>::vector_value(const std::vector<Point<dim> > &points,
  std::vector<std::array<double, Equations<equationsType, dim>::n_components> >&result) const
{
  // The result is a two-dimensional array, first dimension is for the integration point, second for the component (density, momentum-x, ...)
  for (unsigned int i = 0; i < points.size(); ++i)
  {
    result[i][0] = 1.;
    result[i][1] = 0.;
    result[i][2] = 0.;
    result[i][3] = 0.;
    result[i][5] = 1.0 / std::sqrt(2.);
    result[i][6] = 1.0 / std::sqrt(2.);
    result[i][7] = 0.;
    if (points[i].norm() < 0.1)
      result[i][4] = 10. / (this->parameters.gas_gamma - 1.0) + 0.5 * (result[i][5] * result[i][5] + result[i][6] * result[i][6] + result[i][7] * result[i][7]);
    else
      result[i][4] = .1 / (this->parameters.gas_gamma - 1.0) + 0.5 * (result[i][5] * result[i][5] + result[i][6] * result[i][6] + result[i][7] * result[i][7]);
  }
}

template class InitialConditionMhdBlast<EquationsTypeMhd, 3>;

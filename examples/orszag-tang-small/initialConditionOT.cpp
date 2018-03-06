#include "complete_elliptic_integrals.h"
#include "initialConditionOT.h"
#include "equationsMhd.h"

template <EquationsType equationsType, int dim>
InitialConditionOT<equationsType, dim>::InitialConditionOT(Parameters<dim>& parameters) :
  InitialCondition<equationsType, dim>(parameters)
{
};


template <EquationsType equationsType, int dim>
void InitialConditionOT<equationsType, dim>::vector_value(const std::vector<Point<dim> > &points,
  std::vector<std::array<double, Equations<equationsType, dim>::n_components> >&result) const
{
  for (unsigned int i = 0; i < points.size(); ++i)
  {
    result[i][0] = 25./(36. * 3.1415926);
    result[i][1] = -result[i][0] * sin(2. * 3.1415926 * points[i][1]);
    result[i][2] = result[i][0] * sin(2. * 3.1415926 * points[i][0]);
    result[i][3] = 0.;
    result[i][4] = (5. / (12. * 3.1415926)) / (this->getParams().gas_gamma - 1.0) + 0.5 * (result[i][5] * result[i][5] + result[i][6] * result[i][6] + result[i][7] * result[i][7]);
    result[i][5] = -std::sqrt(1. / (4. * 3.1415926)) * sin(2. * 3.1415926 * points[i][1]);
    result[i][6] = std::sqrt(1. / (4. * 3.1415926)) * sin(2. * 3.1415926 * points[i][0]);
    result[i][7] = 0.0;
  }
}

template class InitialConditionOT<EquationsTypeMhd, 3>;

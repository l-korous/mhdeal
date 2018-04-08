#include "completeEllipticIntegrals.h"
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
  double pressure = 5. / (12. * M_PI);
  double B_0 = std::sqrt(1. / (4. * M_PI));
  for (unsigned int i = 0; i < points.size(); ++i)
  {
    result[i][0] = 25. / (36. * M_PI);
    result[i][1] = -result[i][0] * sin(2. * M_PI * points[i][1]);
    result[i][2] = result[i][0] * sin(2. * M_PI * points[i][0]);
    result[i][3] = 0.;
    result[i][5] = -B_0 * sin(2. * M_PI * points[i][1]);
    result[i][6] = B_0 * sin(4. * M_PI * points[i][0]);
    result[i][7] = 0.0;
    result[i][4] = (pressure / (this->getParams().gas_gamma - 1.0)) + 0.5 * (result[i][5] * result[i][5] + result[i][6] * result[i][6] + result[i][7] * result[i][7]) + 
      0.5 * (result[i][1] * result[i][1] + result[i][2] * result[i][2] + result[i][3] * result[i][3]) / result[i][0];
  }
}

template class InitialConditionOT<EquationsTypeMhd, 3>;

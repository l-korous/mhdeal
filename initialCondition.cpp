#include "initialCondition.h"
#include "equationsEuler.h"
#include "equationsMhd.h"

template <EquationsType equationsType, int dim>
InitialCondition<equationsType, dim>::InitialCondition(Parameters<dim>& parameters) : Function<dim>(Equations<equationsType, dim>::n_components), parameters(parameters)
{
};

template<>
double InitialCondition<EquationsTypeMhd, 3>::value(const Point<3> &p, const unsigned int component) const
{
  switch (component) {
  case 0:
    return 1.;
    break;
  case 1:
    return 0.;
    break;
  case 2:
    return 0.;
    break;
  case 3:
    return 0.;
    break;
  case 4:
    return 1.0 / std::sqrt(2.);
    break;
  case 5:
    return 1.0 / std::sqrt(2.);
    break;
  case 6:
    return 0.;
    break;
  case 7:
    if (p.norm() < 0.1)
      return 10. / (parameters.gas_gamma - 1.0) + 0.5;
    else
      return 0.1 / (parameters.gas_gamma - 1.0) + 0.5;
    break;
  }
};

template class InitialCondition<EquationsTypeMhd, 3>;
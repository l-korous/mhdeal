#include "initialCondition.h"
#include "equationsEuler.h"
#include "equationsMhd.h"

template <EquationsType equationsType, int dim>
InitialCondition<equationsType, dim>::InitialCondition(Parameters<dim>& parameters) : Function<dim>(Equations<equationsType, dim>::n_components), parameters(parameters)
{
};


template <EquationsType equationsType, int dim>
void InitialCondition<equationsType, dim>::vector_value(const std::vector<Point<dim> > &points, std::vector<Vector<double> > & result) const
{
  // The result is a two-dimensional array, first dimension is for the integration point, second for the component (density, momentum-x, ...)
  for (int i = 0; i < points.size(); i++)
  {
    for (int j = 0; j < Equations<equationsType, dim>::n_components; j++)
    {
      switch (j) {
      case 0:
        result[i][j] = 1.;
        break;
      case 1:
        result[i][j] = 0.;
        break;
      case 2:
        result[i][j] = 0.;
        break;
      case 3:
        result[i][j] = 0.;
        break;
      case 4:
        result[i][j] = 1.0 / std::sqrt(2.);
        break;
      case 5:
        result[i][j] = 1.0 / std::sqrt(2.);
        break;
      case 6:
        result[i][j] = 0.;
        break;
      case 7:
        if (points[i].norm() < 0.1)
          result[i][j] = 10. / (parameters.gas_gamma - 1.0) + 0.5;
        else
          result[i][j] = 0.1 / (parameters.gas_gamma - 1.0) + 0.5;
        break;
      }
    }
  }
}

template class InitialCondition<EquationsTypeMhd, 3>;
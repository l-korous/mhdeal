#include "complete_elliptic_integrals.h"
#include "initialCondition.h"
#include "boundaryConditions.h"
#include "equationsMhd.h"

template <EquationsType equationsType, int dim>
InitialCondition<equationsType, dim>::InitialCondition(Parameters<dim>& parameters) : Function<dim>(Equations<equationsType, dim>::n_components), parameters(parameters)
{
};


template <EquationsType equationsType, int dim>
void InitialCondition<equationsType, dim>::vector_value(const std::vector<Point<dim> > &points, std::vector<Vector<double> > & result) const
{
  for (unsigned int i = 0; i < points.size(); ++i)
  {
    result[i][0] = 1.;
    result[i][1] = 0.;
    result[i][2] = 0.;
    result[i][3] = 0.;
    result[i][4] = 0.;
    result[i][5] = 0.;
    result[i][6] = 0.;
    result[i][7] = compute_energy_from_given_pressure(result[i], 10.);
  }
}

template <EquationsType equationsType, int dim>
double InitialCondition<equationsType, dim>::compute_energy_from_given_pressure(const dealii::Vector<double>& W, double pressure) const
{
  return ((1. / (this->parameters.gas_gamma - 1.0)) * pressure) + Equations<equationsType, dim>::compute_kinetic_energy(W) + Equations<equationsType, dim>::compute_magnetic_energy(W);
}

template <EquationsType equationsType, int dim>
Parameters<dim>& InitialCondition<equationsType, dim>::getParams() const
{
  return parameters;
}

template class InitialCondition<EquationsTypeMhd, 3>;

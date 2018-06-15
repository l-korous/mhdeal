#include "boundaryCondition.h"
#include "equationsMhd.h"

template<>
BoundaryCondition<EquationsTypeMhd, 3>::BoundaryCondition(Parameters<3>& parameters) : parameters(parameters) {};

template <EquationsType equationsType, int dim>
void BoundaryCondition<equationsType, dim>::bc_vector_value(int boundary_no, const Point<dim> &point, InputVector &result, const InputVector &W_plus) const
{
  for (unsigned int di = 0; di < Equations<EquationsTypeMhd, dim>::n_components; ++di)
    result[di] = W_plus[di];
}

template <EquationsType equationsType, int dim>
double BoundaryCondition<equationsType, dim>::compute_energy_from_given_pressure(const InputVector &W, double pressure) const
{
  return ((1. / (this->parameters.gas_gamma - 1.0)) * pressure) + Equations<equationsType, dim>::compute_kinetic_energy(W) + Equations<equationsType, dim>::compute_magnetic_energy(W);
}

template class BoundaryCondition<EquationsTypeMhd, 3>;

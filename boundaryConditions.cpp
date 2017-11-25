#include "boundaryConditions.h"
#include "equationsMhd.h"

template<>
BoundaryConditions<EquationsTypeMhd, 3>::BoundaryConditions(Parameters<3>& parameters) : parameters(parameters) {};

template <EquationsType equationsType, int dim>
void BoundaryConditions<equationsType, dim>::bc_vector_value(int boundary_no, const Point<dim> &point, InputVector &result, const InputVector &W_plus) const
{
  for (unsigned int di = 0; di < Equations<EquationsTypeMhd, 3>::n_components; ++di)
    result[di] = W_plus[di];
}

template <EquationsType equationsType, int dim>
bool BoundaryConditions<equationsType, dim>::should_limit_this_boundary_id(int boundary_no) const
{
  return true;
}

template <EquationsType equationsType, int dim>
double BoundaryConditions<equationsType, dim>::compute_energy_from_given_pressure(const InputVector &W, double pressure) const
{
  return ((1. / (this->parameters.gas_gamma - 1.0)) * pressure) + Equations<equationsType, dim>::compute_kinetic_energy(W) + Equations<equationsType, dim>::compute_magnetic_energy(W);
}

template class BoundaryConditions<EquationsTypeMhd, 3>;

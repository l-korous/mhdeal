#include "boundaryConditions.h"
#include "equationsMhd.h"

template<>
BoundaryConditions<EquationsTypeMhd, 3>::BoundaryConditions(Parameters<3>& parameters) : parameters(parameters) {};

template <EquationsType equationsType, int dim>
template <typename InputVector>
void BoundaryConditions<equationsType, dim>::bc_vector_value(int boundary_no, const Point<dim> &point, InputVector &result, const InputVector &W_plus) const
{
  for (unsigned int di = 0; di < Equations<EquationsTypeMhd, 3>::n_components; ++di)
    result[di] = W_plus[di];
  if(boundary_no == 0)
    result[7] = compute_energy_from_given_pressure(W_plus, 11.);
}

template <EquationsType equationsType, int dim>
template <typename InputVector>
typename InputVector::value_type BoundaryConditions<equationsType, dim>::compute_energy_from_given_pressure(const InputVector &W, double pressure) const
{
  return ((1. / (this->parameters.gas_gamma - 1.0)) * pressure) + Equations<equationsType, dim>::compute_kinetic_energy(W) + Equations<equationsType, dim>::compute_magnetic_energy(W);
}

template class BoundaryConditions<EquationsTypeMhd, 3>;
template void BoundaryConditions<EquationsTypeMhd, 3>::bc_vector_value(int boundary_no, const Point<3> &point, dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1> &result, const dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1> &W) const;
template void BoundaryConditions<EquationsTypeMhd, 3>::bc_vector_value(int boundary_no, const Point<3> &point, dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double>, false, 1> &result, const dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double>, false, 1> &W) const;
template typename dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1>::value_type BoundaryConditions<EquationsTypeMhd, 3>::compute_energy_from_given_pressure(const dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1> &W, double pressure) const;
template typename dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double>, false, 1>::value_type BoundaryConditions<EquationsTypeMhd, 3>::compute_energy_from_given_pressure(const dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double>, false, 1> &W, double pressure) const;
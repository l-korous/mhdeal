#include "boundaryCondition.h"
#include "equationsMhd.h"

template<>
BoundaryCondition<EquationsTypeMhd, 3>::BoundaryCondition(Parameters<3>& parameters) : parameters(parameters) {};

template <EquationsType equationsType, int dim>
void BoundaryCondition<equationsType, dim>::bc_vector_value(int boundary_no, const Point<dim> &point, const Tensor<1, dim> &normal, 
  values_vector &result, const grad_vector &grads, const values_vector &W_plus, double time, typename DoFHandler<dim>::active_cell_iterator&) const
{
  for (unsigned int di = 0; di < Equations<EquationsTypeMhd, dim>::n_components; ++di)
    result[di] = W_plus[di];
}

template class BoundaryCondition<EquationsTypeMhd, 3>;

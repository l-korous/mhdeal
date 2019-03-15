#include "completeEllipticIntegrals.h"
#include "boundaryConditionCS.h"
#include "equationsMhd.h"

template <int dim>
BoundaryConditionCSWithVortices<dim>::BoundaryConditionCSWithVortices(Parameters<dim>& parameters, CSParameters& cs_parameters) :
  BoundaryCondition<EquationsTypeMhd, dim>(parameters), cs_parameters(cs_parameters)
{

}



template <int dim>
void BoundaryConditionCSWithVortices<dim>::bc_vector_value(int boundary_no, const Point<dim> &point, const Tensor<1, dim> &normal, 
  values_vector &result, const grad_vector &grads, const values_vector &values, double time, typename DoFHandler<dim>::active_cell_iterator&) const
{
  // For other than z=0 boundaries, we use do-nothing
  if (point[2] > SMALL)
  {
    for (unsigned int di = 0; di < Equations<EquationsTypeMhd, dim>::n_components; ++di)
      result[di] = values[di];
    return;
  }

  result[0] = std::sin(60*time)*std::sin(60 * time);//density wave
  result[1] = 0;
  result[2] = 0.;
  result[3] = 0.0;

  // energy density
  result[4] = values[4];

  result[5] = 0.;
  result[6] = values[6];
  result[7] = 0.;
}

template <int dim>
BoundaryConditionCSFree<dim>::BoundaryConditionCSFree(Parameters<dim>& parameters, CSParameters& cs_parameters) :
  BoundaryCondition<EquationsTypeMhd, dim>(parameters)
{
}

template <int dim>
void BoundaryConditionCSFree<dim>::bc_vector_value(int boundary_no, const Point<dim> &point, const Tensor<1, dim> &normal, 
  values_vector &result, const grad_vector &grads, const values_vector &values, double time, typename DoFHandler<dim>::active_cell_iterator&) const
{
  for (unsigned int di = 0; di < Equations<EquationsTypeMhd, dim>::n_components; ++di)
    result[di] = values[di];
  return;
}

template <int dim>
BoundaryConditionCSTest<dim>::BoundaryConditionCSTest(Parameters<dim>& parameters, CSParameters& cs_parameters) :
  BoundaryCondition<EquationsTypeMhd, dim>(parameters), cs_parameters(cs_parameters)
{
}

template <int dim>
void BoundaryConditionCSTest<dim>::bc_vector_value(int boundary_no, const Point<dim> &point, const Tensor<1, dim> &normal, 
  values_vector &result, const grad_vector &grads, const values_vector &values, double time, typename DoFHandler<dim>::active_cell_iterator& cell) const
{
  for (unsigned int di = 0; di < Equations<EquationsTypeMhd, dim>::n_components; ++di)
    result[di] = values[di];
  return;
  // Density the same.
  result[0] = values[0];

  // Velocities are zero on the bottom boundary, otherwise the same as inside.
  if (boundary_no != 4)
  {
    result[1] = values[1];
    result[2] = values[2];
    result[3] = values[3];
  }
  else
    result[1] = result[2] = result[3] = 0.;

  // From divergence-free constraint. Here for x-direction:
  // \frac{\partial B^{'}_x}{\partial x} = -\left(\frac{\partial B_y}{\partial y} + \frac{\partial B_z}{\partial z} \right)
  // We want a linear reconstruction B^{'} = B + d * \frac{\partial B^{'}_x}{\partial x}
  // d will be taken as the elementh length in the direction.
  // Assumption: we have cubes.
  double d = std::pow(cell->measure(), 1./3.);
  // In order to have value (and not just the derivative).
  // \left|B^{'}_x}\right| = \left|B_x}\right|

  double derivative;
  if (std::abs(normal[0]) < SMALL)
  {
    // z-direction
    if (std::abs(normal[1]) < SMALL)
    {
      derivative = -grads[5][0] - grads[6][1];
      result[5] = (0.5 * d * grads[5][2] * normal[2]) + values[5];
      result[6] = (0.5 * d * grads[6][2] * normal[2]) + values[6];
      result[7] = (0.5 * d * derivative * normal[2]) + values[7];
    }
    // y-direction
    else
    {
      derivative = -grads[5][0] - grads[7][2];
      result[5] = (0.5 * d * grads[5][1] * normal[1]) + values[5];
      result[6] = (0.5 * d * derivative * normal[1]) + values[6];
      result[7] = (0.5 * d * grads[7][1] * normal[1]) + values[7];
    }
  }
  // x-direction
  else
  {
    derivative = -grads[6][1] - grads[7][2];
    result[5] = (0.5 * d * derivative * normal[0]) + values[5];
    result[6] = (0.5 * d * grads[6][0] * normal[0]) + values[6];
    result[7] = (0.5 * d * grads[7][0] * normal[0]) + values[7];
  }

  result[4] = values[4];

  if (boundary_no = 4) { result[5] = result[7] = 0.;             //wave motion
						result[0]=sin(time*10.)*sin(time*10);
  }
}

template <int dim>
BoundaryConditionCSInitialState<dim>::BoundaryConditionCSInitialState(Parameters<dim>& parameters, CSParameters& cs_parameters) :
  BoundaryCondition<EquationsTypeMhd, dim>(parameters), cs_parameters(cs_parameters),
  ic(parameters, cs_parameters)
{
}

template <int dim>
void BoundaryConditionCSInitialState<dim>::bc_vector_value(int boundary_no, const Point<dim> &p, const Tensor<1, dim> &normal, 
  values_vector &result, const grad_vector &grads, const values_vector &values, double time, typename DoFHandler<dim>::active_cell_iterator& cell) const
{
  std::vector<Point<dim> > points;
  Point<dim> p1;
  for (unsigned int di = 0; di < dim; ++di)
  {
    double length_in_direction = (this->parameters.corner_b[di] - this->parameters.corner_a[di]) / this->parameters.refinements[di];
    p1[di] = p[di] + (0.5 * normal[di] * length_in_direction);
  }
  points.push_back(p1);
  std::array<double, Equations<EquationsTypeMhd, dim>::n_components> vla;
  std::vector<std::array<double, Equations<EquationsTypeMhd, dim>::n_components> > vl;
  vl.push_back(vla);
  ic.vector_value(points, vl);
  for (unsigned int di = 0; di < Equations<EquationsTypeMhd, dim>::n_components; ++di)
    result[di] = vl[0][di];
  //result[0] = values[0];
  result[1] = result[2] = result[3] = 0.;
  //result[4] = Equations<EquationsTypeMhd, dim>::compute_energy_from_pressure(result, this->cs_parameters.beta, this->parameters);
}

template class BoundaryConditionCSWithVortices<3>;
template class BoundaryConditionCSFree<3>;
template class BoundaryConditionCSInitialState<3>;
template class BoundaryConditionCSTest<3>;

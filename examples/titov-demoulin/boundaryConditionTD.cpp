#include "completeEllipticIntegrals.h"
#include "boundaryConditionTD.h"
#include "equationsMhd.h"

template <int dim>
BoundaryConditionTDWithVortices<dim>::BoundaryConditionTDWithVortices(Parameters<dim>& parameters, TitovDemoulinParameters& td_parameters) :
  BoundaryCondition<EquationsTypeMhd, dim>(parameters), td_parameters(td_parameters)
{
  this->eps = sqrt(1. - (td_parameters.d / td_parameters.R) * (td_parameters.d / td_parameters.R));
  y_1 = this->eps * td_parameters.R;
  y_2 = -this->eps * td_parameters.R;
}

template <int dim>
double BoundaryConditionTDWithVortices<dim>::r_1_bar(double x, double y) const
{
  return sqrt(this->eps * this->eps * ((y - this->y_1) * (y - this->y_1)) + (x * x));
}

template <int dim>
double BoundaryConditionTDWithVortices<dim>::r_2_bar(double x, double y) const
{
  return sqrt(this->eps * this->eps * ((y - this->y_2) * (y - this->y_2)) + (x * x));
}

template <int dim>
double BoundaryConditionTDWithVortices<dim>::omega_1(double x, double y) const
{
  return 0.5 * (1. + tanh(1. - r_1_bar(x, y)));
}

template <int dim>
double BoundaryConditionTDWithVortices<dim>::omega_2(double x, double y) const
{
  return 0.5 * (1. + tanh(1. - r_2_bar(x, y)));
}

template <int dim>
double BoundaryConditionTDWithVortices<dim>::omega(double time) const
{
  return (this->td_parameters.omega_0 / 2.) * (1. + tanh((time - this->td_parameters.t_drive) / this->td_parameters.t_ramp));
}

template <int dim>
void BoundaryConditionTDWithVortices<dim>::bc_vector_value(int boundary_no, const Point<dim> &point, const Tensor<1, dim> &normal, 
  values_vector &result, const grad_vector &grads, const values_vector &values, double time, typename DoFHandler<dim>::active_cell_iterator&) const
{
  // For other than z=0 boundaries, we use do-nothing
  if (point[2] > SMALL)
  {
    for (unsigned int di = 0; di < Equations<EquationsTypeMhd, dim>::n_components; ++di)
      result[di] = values[di];
    return;
  }

  double x = point[0], y = point[1];
  result[0] = values[0];
  result[1] = result[0] * (-omega(time) * this->eps * ((y - y_1) * omega_1(x, y) + (y - y_2) * omega_2(x, y)));
  result[2] = result[0] * (omega(time) * (x / this->eps) * (omega_1(x, y) + omega_2(x, y)));
  result[3] = 0.0;

  // energy density
  result[4] = values[4];

  result[5] = values[5];
  result[6] = values[6];
  result[7] = values[7];
}

template <int dim>
BoundaryConditionTDFree<dim>::BoundaryConditionTDFree(Parameters<dim>& parameters, TitovDemoulinParameters& td_parameters) :
  BoundaryCondition<EquationsTypeMhd, dim>(parameters)
{
}

template <int dim>
void BoundaryConditionTDFree<dim>::bc_vector_value(int boundary_no, const Point<dim> &point, const Tensor<1, dim> &normal, 
  values_vector &result, const grad_vector &grads, const values_vector &values, double time, typename DoFHandler<dim>::active_cell_iterator&) const
{
  for (unsigned int di = 0; di < Equations<EquationsTypeMhd, dim>::n_components; ++di)
    result[di] = values[di];
  return;
}

template <int dim>
BoundaryConditionTDTest<dim>::BoundaryConditionTDTest(Parameters<dim>& parameters, TitovDemoulinParameters& td_parameters) :
  BoundaryCondition<EquationsTypeMhd, dim>(parameters)
{
}

template <int dim>
void BoundaryConditionTDTest<dim>::bc_vector_value(int boundary_no, const Point<dim> &point, const Tensor<1, dim> &normal, 
  values_vector &result, const grad_vector &grads, const values_vector &values, double time, typename DoFHandler<dim>::active_cell_iterator& cell) const
{
  // Density the same.
  result[0] = values[0];

  // Velocities are zero on the bottom boundary, otherwise the same as inside.
  if (point[2] > SMALL)
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
      result[5] = values[5];
      result[6] = values[6];
      result[7] = 0.5 * d * derivative * normal[2] + values[7];
    }
    // y-direction
    else
    {
      derivative = -grads[5][0] - grads[7][2];
      result[5] = values[5];
      result[6] = 0.5 * d * derivative * normal[1] + values[6];
      result[7] = values[7];
    }
  }
  // x-direction
  else
  {
    derivative = -grads[6][1] - grads[7][2];
    result[5] = 0.5 * d * derivative * normal[0] + values[5];
    result[6] = values[6];
    result[7] = values[7];
  }

  result[4] = values[4];
}

template <int dim>
BoundaryConditionTDInitialState<dim>::BoundaryConditionTDInitialState(Parameters<dim>& parameters, TitovDemoulinParameters& td_parameters) :
  BoundaryCondition<EquationsTypeMhd, dim>(parameters), td_parameters(td_parameters),
  ic(parameters, td_parameters)
{
}

template <int dim>
void BoundaryConditionTDInitialState<dim>::bc_vector_value(int boundary_no, const Point<dim> &p, const Tensor<1, dim> &normal, 
  values_vector &result, const grad_vector &grads, const values_vector &values, double time, typename DoFHandler<dim>::active_cell_iterator&) const
{
  // For other than z=0 boundaries, we use do-nothing
  if (p[2] > SMALL)
  {
    for (unsigned int di = 0; di < Equations<EquationsTypeMhd, dim>::n_components; ++di)
      result[di] = values[di];
    return;
  }

  std::vector<Point<dim> > points;
  points.push_back(p);
  std::array<double, Equations<EquationsTypeMhd, dim>::n_components> vla;
  std::vector<std::array<double, Equations<EquationsTypeMhd, dim>::n_components> > vl;
  vl.push_back(vla);
  ic.vector_value(points, vl);
  for (unsigned int di = 0; di < Equations<EquationsTypeMhd, dim>::n_components; ++di)
    result[di] = vl[0][di];
}

template class BoundaryConditionTDWithVortices<3>;
template class BoundaryConditionTDFree<3>;
template class BoundaryConditionTDInitialState<3>;
template class BoundaryConditionTDTest<3>;

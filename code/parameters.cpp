#include "parameters.h"

template <EquationsType equationsType, int dim>
InitialCondition<equationsType, dim>::InitialCondition() : Function<dim>(Equations<equationsType, dim>::n_components)
{
};

template <EquationsType equationsType, int dim>
double InitialCondition<equationsType, dim>::value(const Point<dim> &p, const unsigned int  component = 0) const
{
  double x = p(0);
  double y = p(1);
  switch (component) {
  case 0:
    return 0.;
    break;
  case 1:
    return 0.;
    break;
  case 2:
    return 10. * (x < -0.7)*(y > 0.3)*(y < 0.45) + (1 - (x < -0.7)*(y > 0.3)*(y < 0.45));
    break;
  case 3:
    return 2.5 * (1.5 - y);
    break;
  }
};

template <EquationsType equationsType, int dim>
BoundaryConditions<equationsType, dim>::BoundaryConditions()
{
  for (unsigned int di = 0; di < Equations<equationsType, dim>::n_components; ++di)
    kind[0][di] = Equations<equationsType, dim>::outflow_boundary;

  for (unsigned int di = 0; di < dim; ++di)
    kind[1][di] = Equations<equationsType, dim>::no_penetration_boundary;

  for (unsigned int di = dim; di < Equations<equationsType, dim>::n_components; ++di)
    kind[1][di] = Equations<equationsType, dim>::outflow_boundary;

  for (unsigned int di = 0; di < Equations<equationsType, dim>::n_components; ++di)
    kind[2][di] = Equations<equationsType, dim>::outflow_boundary;

  for (unsigned int di = 0; di < dim; ++di)
    kind[3][di] = Equations<equationsType, dim>::no_penetration_boundary;

  for (unsigned int di = dim; di < Equations<equationsType, dim>::n_components; ++di)
    kind[3][di] = Equations<equationsType, dim>::outflow_boundary;

  for (unsigned int di = 0; di < dim; ++di)
    kind[4][di] = Equations<equationsType, dim>::no_penetration_boundary;

  for (unsigned int di = dim; di < Equations<equationsType, dim>::n_components; ++di)
    kind[4][di] = Equations<equationsType, dim>::outflow_boundary;
};

template <EquationsType equationsType, int dim>
void BoundaryConditions<equationsType, dim>::bc_vector_value(int boundary_id, const std::vector<Point<dim> > &points, std::vector<Vector<double> > & result) const
{
  for (int j = 0; j < Equations<equationsType, dim>::n_components; j++)
  {
    for (int i = 0; i < points.size(); i++)
      result[i][j] = 0.;
  }
}

template <int dim>
Parameters<dim>::Parameters()
{
  this->mesh_filename = "slide.inp";
  this->final_time = 10.;
  this->time_step = .01;
  this->theta = 0.5;

  this->output = OutputType::quiet_solver;
  this->solver = SolverType::direct;
  this->linear_residual = 1e-10;
  this->max_iterations = 300;
  this->ilut_fill = 1.5;
  this->ilut_drop = 1e-6;
  this->ilut_atol = 1e-6;
  this->ilut_rtol = 1.0;

  this->polynomial_order = 0;
  this->max_nonlinear_iterations = 30;
  this->nonlinear_residual_norm_threshold = 1e-8;

  this->output_step = 0.01;
  this->schlieren_plot = false;

  this->stabilization_kind = StabilizationKind::constant_stabilization;
  this->stabilization_value = 1.;

  this->is_stationary = false;
}

template class Parameters<2>;
template class Parameters<3>;

template class BoundaryConditions<EquationsTypeEuler, 2>;
template class BoundaryConditions<EquationsTypeEuler, 3>;

template class InitialCondition<EquationsTypeEuler, 2>;
template class InitialCondition<EquationsTypeEuler, 3>;
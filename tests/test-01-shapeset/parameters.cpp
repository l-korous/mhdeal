#include "parameters.h"
#include "equations.h"

template <int dim>
#ifdef HAVE_MPI
Parameters<dim>::Parameters(parallel::distributed::Triangulation<dim> &triangulation)
#else
Parameters<dim>::Parameters(Triangulation<dim> &triangulation)
#endif
{
  this->debug = false;

  // Two corners of the hyper-rectangle
  // - corner A
  this->corner_a = Point<dim>(.0, .0, .0);
  // - and corner B which should be the farthest one from corner A
  this->corner_b = Point<dim>(1., 1., 1.);
  // Refinements in x-, y-, and z- coordinates.
  this->refinements = { 1, 1, 1 };
  // deal.II function that takes the above attributes and returns the triangulation (the first parameter, passed by reference).
  GridGenerator::subdivided_hyper_rectangle(triangulation, this->refinements, this->corner_a, this->corner_b, true);

  this->polynomial_order = 1;
  this->quadrature_order = 2;

  this->output_matrix = true;
  this->output = OutputType::quiet_solver;
  this->output_rhs = true;
  this->output_solution = true;

  this->solver = direct;
  this->linear_residual = 1e-10;
  this->max_iterations = 10000;
  this->ilut_fill = 1.5;
  this->ilut_drop = 1e-6;
  this->ilut_atol = 1e-6;
  this->ilut_rtol = 1.0;
}

template class Parameters<2>;
template class Parameters<3>;

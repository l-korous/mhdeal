#include "parameters.h"
#include "equationsMhd.h"

template <int dim>
#ifdef HAVE_MPI
Parameters<dim>::Parameters(parallel::distributed::Triangulation<dim> &triangulation)
#else
Parameters<dim>::Parameters(Triangulation<dim> &triangulation)
#endif
{
  this->debug = false;
  this->debug_limiter = false;
  this->initCond = 0;

  // Two corners of the hyper-rectangle
  // - corner A
  this->corner_a = Point<dim>(0, 0, 0);
  // - and corner B which should be the farthest one from corner A
  this->corner_b = Point<dim>(.35 , .35, .01);
  // Refinements in x-, y-, and z- coordinates.
  this->refinements = { 60, 60, 1 };
  // deal.II function that takes the above attributes and returns the triangulation (the first parameter, passed by reference).
  GridGenerator::subdivided_hyper_rectangle(triangulation, this->refinements, this->corner_a, this->corner_b, true);

  this->time_step = 1.e-4;
  this->final_time = 10.;

  this->theta = 0.0;

  this->polynomial_order_dg = 1;
  this->polynomial_order_hdiv = 1;

  this->quadrature_order = 5;

  this->output_matrix = false;
  this->output = OutputType::quiet_solver;
  this->output_rhs = false;
  this->output_solution = false;

  this->output_step = -1.;

  this->solver = gmres;
  this->linear_residual = 1e-10;
  this->max_iterations = 10000;
  this->ilut_fill = 1.5;
  this->ilut_drop = 1e-6;
  this->ilut_atol = 1e-6;
  this->ilut_rtol = 1.0;

  this->gas_gamma = 1.4;

  this->newton_damping = 1.;
  this->newton_max_iterations = 30;
  this->newton_residual_norm_threshold = 1e-8;

  this->num_flux_type = hlld;
  this->lax_friedrich_stabilization_value = 1.;

  this->is_stationary = false;

  this->needs_gradients = false;
  this->needs_forcing = false;
}

template class Parameters<2>;
template class Parameters<3>;

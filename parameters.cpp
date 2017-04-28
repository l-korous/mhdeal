#include "parameters.h"
#include "equationsEuler.h"
#include "equationsMhd.h"

template <int dim>
#ifdef HAVE_MPI
Parameters<dim>::Parameters(parallel::distributed::Triangulation<dim> &triangulation)
#else
Parameters<dim>::Parameters(Triangulation<dim> &triangulation)
#endif
{
  this->debug = false;

  this->corner_a = Point<dim>(.0, .0, .0);
  this->corner_b = Point<dim>(.5, .5, .001);
  this->refinements = { 5, 5, 5 };
  GridGenerator::subdivided_hyper_rectangle(triangulation, this->refinements, this->corner_a, this->corner_b);

  this->time_step = 1e-6;
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

  this->initial_step = true;

}

template class Parameters<2>;
template class Parameters<3>;
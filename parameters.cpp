#include "parameters.h"
#include "equationsMhd.h"

template <int dim>
#ifdef HAVE_MPI
Parameters<dim>::Parameters(parallel::distributed::Triangulation<dim> &triangulation)
#else
Parameters<dim>::Parameters(Triangulation<dim> &triangulation)
#endif
{
  this->num_flux_type = hlld;
  this->cfl_constant = .1;
  this->corner_a = Point<dim>(-0.3, -.3, 0.);
  this->corner_b = Point<dim>(.3 , .3, .01);
  this->refinements = { 100, 100, 1 };
  this->quadrature_order = 5;
  this->polynomial_order_dg = 0;
  this->polynomial_order_hdiv = 0;

  this->theta = 0.0;
  this->postprocess_in_newton_loop = true;
  this->patches = 2;
  this->output_step = -1.e-3;

  this->debug = false;
  this->debug_limiter = false;
  this->debug_dofs = false;
  this->initCond = 0;

  this->output_matrix = false;
  this->output = OutputType::quiet_solver;
  this->output_rhs = false;
  this->output_solution = false;

  this->snapshot_step = 1.; 
  
  this->time_step = 1.e-5;
  this->final_time = 10.;
  
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

  this->lax_friedrich_stabilization_value = 0.;

  this->is_stationary = false;

  this->needs_gradients = false;
  this->needs_forcing = false;
  
  GridGenerator::subdivided_hyper_rectangle(triangulation, this->refinements, this->corner_a, this->corner_b, true);
}

template class Parameters<2>;
template class Parameters<3>;

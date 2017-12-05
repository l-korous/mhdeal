#include "parameters.h"
#include "equationsMhd.h"
#include "DealiiExtensions.h"

template <int dim>
#ifdef HAVE_MPI
Parameters<dim>::Parameters(parallel::distributed::Triangulation<dim> &triangulation)
#else
Parameters<dim>::Parameters(Triangulation<dim> &triangulation)
#endif
{
  this->initCond = 0;
  this->num_flux_type = hlld;
  this->cfl_constant = .05;
  this->corner_a = Point<dim>(-0.5, -0.75, 0.);
  this->corner_b = Point<dim>(.5, .75, .01);
  this->refinements = { 100, 100, 1 };
  this->quadrature_order = 5;
  this->initial_quadrature_order = 10;
  this->polynomial_order_dg = 1;
  this->polynomial_order_hdiv = 0;
  this->limit_in_nonlin_loop = false;
  periodic_boundaries = { { 0, 1, 0 },{ 2, 3, 1 } };
  
  this->initial_and_max_newton_damping = .9;
  this->decrease_factor = .8;
  this->increase_factor = 1.1;


  this->patches = 2;
  this->output_step = -1.e-3;

  this->debug = false;
  this->debug_limiter = false;
  this->debug_dofs = false;

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

  this->newton_max_iterations = 100;
  this->newton_residual_norm_threshold = 1e-8;

  this->lax_friedrich_stabilization_value = 0.;

  this->is_stationary = false;

  this->needs_gradients = false;
  this->needs_forcing = false;

  GridGenerator::subdivided_hyper_rectangle(triangulation, this->refinements, this->corner_a, this->corner_b, true);
  std::vector<DealIIExtensions::FacePair<dim> > matched_pairs;
  for (std::vector<std::array<int, 3> >::const_iterator it = periodic_boundaries.begin(); it != periodic_boundaries.end(); it++)
    dealii::GridTools::collect_periodic_faces(triangulation, (*it)[0], (*it)[1], (*it)[2], matched_pairs);
  triangulation.add_periodicity(matched_pairs);
}

template class Parameters<2>;
template class Parameters<3>;
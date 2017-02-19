#ifndef _PARAMETERS_H
#define _PARAMETERS_H

#include "util.h"
#include "equations.h"

template <int dim>
class Parameters
{
public:
#ifdef HAVE_MPI
  Parameters(parallel::distributed::Triangulation<dim> &triangulation);
#else
  Parameters(Triangulation<dim> &triangulation);
#endif

  // Flux
  enum NumFluxType { hlld, lax_friedrich };
  NumFluxType num_flux_type;
  double lax_friedrich_stabilization_value;

  // Output
  double output_step;
  double output_matrix;
  double output_rhs;
  double output_solution;

  double gas_gamma;

  // Nonlinear solver
  double newton_damping;
  int max_nonlinear_iterations;
  double nonlinear_residual_norm_threshold;

  // Linear solver
  enum SolverType { gmres, direct };
  SolverType solver;
  enum  OutputType { quiet_solver, verbose_solver };
  OutputType output;
  double linear_residual;
  int max_iterations;
  double ilut_fill;
  double ilut_atol;
  double ilut_rtol;
  double ilut_drop;

  // Global
  double time_step, final_time, theta;
  bool initial_step;
  double time_step_after_initialization, initialization_time, theta_after_initialization;
  bool is_stationary;
  int polynomial_order_dg;
  int polynomial_order_hdiv;
  bool needs_gradients;

  // Mesh - string input
  std::string mesh_filename;
  int mesh_extrusion_slices;

  // Mesh - generation
  Point<dim> corner_a, corner_b;
  std::vector<unsigned int> refinements;

  // Debugging purposes
  bool debug;
  bool debug_neighbor;
};

#endif
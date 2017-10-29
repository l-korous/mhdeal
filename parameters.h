#ifndef _PARAMETERS_H
#define _PARAMETERS_H

#include "util.h"
#include "equations.h"

template <int dim>
class Parameters
{
public:
  // Parameters constructor takes a triangulation as an attribute (passed by reference), and the constructor is responsible for filling out the triangulation.
#ifdef HAVE_MPI
  Parameters(parallel::distributed::Triangulation<dim> &triangulation);
#else
  Parameters(Triangulation<dim> &triangulation);
#endif

  // Flux enumeration - for potentially adding more fluxes, decision which one to use is then made in Equations<>::numerical_normal_flux.
  enum NumFluxType { hlld, lax_friedrich };
  NumFluxType num_flux_type;
  // A special value for lax_friedrich
  double lax_friedrich_stabilization_value;

  bool postprocess_in_newton_loop;
  
  // Output step - either < 0 (output all steps), or > 0 (time difference between two outputs)
  double output_step, snapshot_step;

  // Output matrix after assemble_system() in Problem::run().
  bool output_matrix;
  // Output rhs after assemble_system() in Problem::run().
  bool output_rhs;
  // Output current_solution after solve() in Problem::run().
  bool output_solution;

  // initial conditions: 0 - MHD Blast, 1 - Titov&Demoulin flux rope
  // Number of patches
  unsigned int initCond, patches; 
  
  // Gas gamma value.
  double gas_gamma;

  // Nonlinear solver
  // Damping factor in Newton
  double newton_damping;
  // Maximum allowed nonlinear iterations count, fail if exceeded
  int newton_max_iterations;
  // Tolerance for nonlinear residual norm, succeed the nonlinear loop if norm < newton_residual_norm_threshold
  double newton_residual_norm_threshold;

  // Linear solver type enumeration
  enum SolverType { gmres, direct };
  // Linear solver type selected
  SolverType solver;
  // Verbosity enumeration
  enum  OutputType { quiet_solver, verbose_solver };
  // Verbosity selected
  OutputType output;

  // Tolerance for linear residual norm, succeed the linear loop if norm < newton_residual_norm_threshold
  double linear_residual;
  // Maximum allowed linear iterations count, succeed the linear loop exceeded
  int max_iterations;
  // Linear solver parameters.
  double ilut_fill;
  double ilut_atol;
  double ilut_rtol;
  double ilut_drop;

  // Global - obvious, theta ~ Theta-scheme for time-discretization
  double time_step, final_time, theta, cfl_constant;
  // Flag whether this is a stationary problem.
  bool is_stationary;
  // Polynomial order for the flow part.
  int polynomial_order_dg;
  // Polynomial order for the mag field part.
  int polynomial_order_hdiv;
  // Quadrature order.
  int quadrature_order, initial_quadrature_order;
  // Do we need calculated gradients? E.g. for resistivity.
  bool needs_gradients;
  // Do we have forcing (source) terms in the equations?
  bool needs_forcing;

  // Mesh - string input
  std::string mesh_filename;
  // If we extrude a 2d mesh into 3d one, how many slices do we want thus created?
  int mesh_extrusion_slices;

  // Mesh - generation - two corners opposite in all dimensions.
  Point<dim> corner_a, corner_b;
  // Mesh - generation - how many refinements (in all dimensions) do we want?
  std::vector<unsigned int> refinements;

  // Debugging purposes
  bool debug, debug_limiter, debug_dofs;
};

#endif

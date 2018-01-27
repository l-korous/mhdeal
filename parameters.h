#ifndef _PARAMETERS_H
#define _PARAMETERS_H

#include "util.h"
#include "equations.h"

template <int dim>
class Parameters
{
public:
  // Parameters constructor takes a triangulation as an attribute (passed by reference), and the constructor is responsible for filling out the triangulation.
  Parameters() {};

  // Flux enumeration - for potentially adding more fluxes, decision which one to use is then made in Equations<>::numerical_normal_flux.
  enum NumFluxType { hlld, lax_friedrich };
  NumFluxType num_flux_type;
  // A special value for lax_friedrich
  double lax_friedrich_stabilization_value;

  // Output step - either < 0 (output all steps), or > 0 (time difference between two outputs)
  double output_step, snapshot_step;

  // Output matrix after assemble_system() in Problem::run().
  bool output_matrix;
  // Output rhs after assemble_system() in Problem::run().
  bool output_rhs;
  // Output limited_solution after solve() in Problem::run().
  bool output_solution;

  // Number of patches
  unsigned int patches; 
  
  // Gas gamma value.
  double gas_gamma;

  // Nonlinear solver
  // Damping factor in Newton
  double initial_and_max_newton_damping;
  double decrease_factor;
  double increase_factor;
  bool limit_in_nonlin_loop;
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

  // Global - obvious
  double time_step, final_time, cfl_constant, initial_and_max_cfl_constant;
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

  // Debugging purposes
  bool debug, debug_limiter, debug_dofs;

  Point<dim> corner_a;
  Point<dim> corner_b;
  std::vector<unsigned int> refinements;
  std::vector<std::array<int, 3> > periodic_boundaries;
};

#endif

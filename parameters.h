#ifndef _PARAMETERS_H
#define _PARAMETERS_H

#include "util.h"
#include "equations.h"

template <int dim>
class Parameters
{
public:
  // Parameters constructor takes a triangulation as an attribute (passed by reference), and the constructor is responsible for filling out the triangulation.
  Parameters() {
    this->limit = true;
    this->output_file_prefix = "";
    this->lax_friedrich_stabilization_value = .5;
    this->snapshot_step = 1.;
    this->time_step = 1.e-5;

    this->debug = false;
    this->output_matrix = false;
    this->output = quiet_solver;
    this->output_rhs = false;
    this->output_solution = false;

    this->solver = gmres;
    this->linear_residual = 1e-10;
    this->max_iterations = 10000;
    this->ilut_fill = 1.5;
    this->ilut_drop = 1e-6;
    this->ilut_atol = 1e-6;
    this->ilut_rtol = 1.0;
  };

  // Use exactly Div-Free space.
  bool use_div_free_space_for_B;

  // Perform iterative improvements.
  bool use_iterative_improvement;

  // Flux enumeration - for potentially adding more fluxes, decision which one to use is then made in Equations<>::numerical_normal_flux.
  enum NumFluxType { hlld, lax_friedrich };
  NumFluxType num_flux_type;
  // A special value for lax_friedrich
  double lax_friedrich_stabilization_value;

  // Output step - either < 0 (output all steps), or > 0 (time difference between two outputs)
  double output_step, snapshot_step;
  // File name
  std::string output_file_prefix;

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
  // Make damping factor variable
  bool automatic_damping;
  // Make CFL coefficient variable
  bool automatic_cfl;
  double decrease_factor;
  double increase_factor;
  bool limit, limit_in_nonlin_loop;
  // Maximum allowed nonlinear iterations count, fail if exceeded
  int newton_max_iterations;
  // Tolerance for nonlinear residual norm, succeed the nonlinear loop if norm < newton_residual_norm_threshold
  double newton_residual_norm_threshold;
  // If relative change between the previous non-linear step and the current one is smaller than this, then we say that the current step was SUCCESSFUL.
  double stagnation_coefficient;
  // If relative change between the previous non-linear step and the current one is larger than this, then we say that the current step was NOT SUCCESSFUL.
  double bad_step_coefficient;

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
  double time_step, final_time, initial_and_max_cfl_coefficient;
  // Polynomial order for the flow part.
  int polynomial_order_dg;
  // Quadrature order.
  int quadrature_order;
  // Debugging purposes
  bool debug;

  Point<dim> corner_a;
  Point<dim> corner_b;
  std::vector<unsigned int> refinements;
  std::vector<std::array<int, 3> > periodic_boundaries;
};

#endif

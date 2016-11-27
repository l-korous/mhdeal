#ifndef _PARAMETERS_H
#define _PARAMETERS_H

#include "util.h"
#include "equations.h"

template <int dim>
class Parameters
{
public:
  Parameters();

  // Flux
  enum StabilizationKind { constant_stabilization, mesh_dependent_stabilization, no_stabilization };
  StabilizationKind stabilization_kind;
  double stabilization_value;

  // Output
  double output_step;

  double gas_gamma;

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
  double time_step, final_time;
  double theta;
  bool is_stationary;
  int polynomial_order;
  int max_nonlinear_iterations;
  double nonlinear_residual_norm_threshold;

  // Mesh
  std::string mesh_filename;
  int MeshSlicesInZDirection;
};

#endif
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
  bool schlieren_plot;
  double output_step;

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

  // Initial condition
  class InitialCondition : public Function<dim>
  {
  public:
    InitialCondition();

    double value(const Point<dim> &p, const unsigned int  component = 0) const;
  };
  InitialCondition initial_condition;

  // Boundary conditions
  static const unsigned int max_n_boundaries = 10;
  class BoundaryConditions
  {
  public:
    typename Equations<dim>::BoundaryKind kind[Equations<dim>::n_components];

    BoundaryConditions();
  };
  BoundaryConditions  boundary_conditions[max_n_boundaries];
  static void bc_vector_value(int boundary_no, const std::vector<Point<dim> > &points, std::vector<Vector<double> >&);
};

#endif
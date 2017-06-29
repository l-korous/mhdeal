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

  // Output matrix after assemble_system() in Problem::run().
  bool output_matrix;
  // Output rhs after assemble_system() in Problem::run().
  bool output_rhs;
  // Output current_solution after solve() in Problem::run().
  bool output_solution;

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

  int polynomial_order;
  // Quadrature order.
  int quadrature_order;

  // Mesh - generation - two corners opposite in all dimensions.
  Point<dim> corner_a, corner_b;
  // Mesh - generation - how many refinements (in all dimensions) do we want?
  std::vector<unsigned int> refinements;

  // Debugging purposes
  bool debug;
};

#endif

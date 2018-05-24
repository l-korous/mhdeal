#ifndef _PROBLEM_H
#define _PROBLEM_H

#include "util.h"
#include "equationsMhd.h"
#include "parameters.h"
#include "initialCondition.h"
#include "boundaryConditions.h"
#include "assemblingUtilities.h"
#include "dealiiExtensions.h"
#include "feDivFree.h"
#include "feTaylor.h"
#include "numericalFlux.h"
#include "slopeLimiter.h"
#include "adaptivity.h"

template <EquationsType equationsType, int dim> class AssemblingUtilities;

// Class that accepts all input from the user, provides interface for output, etc.
// Should not be changed.
template <EquationsType equationsType, int dim>
class Problem
{
public:
  Problem(Parameters<dim>& parameters, Equations<equationsType, dim>& equations,
#ifdef HAVE_MPI
    parallel::distributed::Triangulation<dim>& triangulation,
#else
    Triangulation<dim>& triangulation,
#endif
    InitialCondition<equationsType, dim>& initial_condition, BoundaryConditions<equationsType, dim>& boundary_conditions, Adaptivity<dim>* adaptivity = 0);
  void run();

  // Technical matters done only once after creation.
  void setup_system();

  // Sets up initial values for current & previous solution.
  void setup_initial_solution();

  // Performs a single global assembly.
  void assemble_system(bool assemble_matrix = true);

  void postprocess();

  // Performs a single global assembly.
  void calculate_cfl_condition();

  // Performs a local assembly for all volumetric contributions on the local cell.
  void assemble_cell_term(FullMatrix<double>& cell_matrix, Vector<double>& cell_rhs, bool assemble_matrix, std::vector<double>& JxW);

  // Performs a local assembly for all surface contributions on the local cell.
  // i.e. face terms calculated on all faces - internal and boundary
  void assemble_face_term(const unsigned int face_no, const FEFaceValuesBase<dim> &fe_v, const FEFaceValuesBase<dim> &fe_v_prev, const FEFaceValuesBase<dim> &fe_v_prev_neighbor, const bool external_face,
    const unsigned int boundary_id, Vector<double>& cell_rhs, std::vector<double>& JxW, std::vector<Tensor<1, dim> >& normals);

  void output_base();
  void output_results() const;
  void output_matrix(TrilinosWrappers::SparseMatrix& mat, const char* suffix) const;
  void output_vector(TrilinosWrappers::MPI::Vector& vec, const char* suffix) const;

  // Solves the assembled system
  void solve();

  // Process output, move to next time step.
  void handle_outputs();

  // Triangulation - passed as a constructor parameter
#ifdef HAVE_MPI
  parallel::distributed::Triangulation<dim> triangulation;
  parallel::distributed::Triangulation<dim>& prev_triangulation;
#else
  Triangulation<dim> triangulation;
  Triangulation<dim>& prev_triangulation;
#endif

  // Equations - passed as a constructor parameter
  Equations<equationsType, dim>& equations;

  // Parameters - passed as a constructor parameter
  Parameters<dim>& parameters;

  // Initical conditions - passed as a constructor parameter
  InitialCondition<equationsType, dim>& initial_condition;

  // Boundary conditions - passed as a constructor parameter
  BoundaryConditions<equationsType, dim>& boundary_conditions;

  // Dofs calculated by this MPI process.
  IndexSet locally_owned_dofs, prev_locally_owned_dofs;

  // Dofs calculated by this MPI process + all Dofs on all neighboring cells.
  IndexSet locally_relevant_dofs, prev_locally_relevant_dofs;

  const MappingQ1<dim> mapping;
  const FESystem<dim> fe;
  DoFHandler<dim> *dof_handler, prev_dof_handler;
  const QGauss<dim> quadrature;
  const QGauss<dim - 1> face_quadrature;

  // Sought solution, the previous one, etc.
  TrilinosWrappers::MPI::Vector     current_limited_solution;
  TrilinosWrappers::MPI::Vector     current_unlimited_solution;
  TrilinosWrappers::MPI::Vector     prev_solution;

  // The system being assembled.
  TrilinosWrappers::MPI::Vector system_rhs;
  TrilinosWrappers::SparseMatrix system_matrix;

  ConstraintMatrix constraints, prev_constraints;
  
  double last_output_time, last_snapshot_time, time;
  // Adaptivity step is increased always, but if mesh is refined, time step does not progress. Therefore, there can be more adaptivity steps in one time step.
  int time_step, adaptivity_step;
  // For CFL.
  double cfl_time_step;
  double max_signal_speed;

  DealIIExtensions::PeriodicCellMap<dim> periodic_cell_map;
  FEValuesExtractors::Vector mag;
  void precalculate_global();
  unsigned int dofs_per_cell;
  unsigned short n_quadrature_points_cell, n_quadrature_points_face;

  // NumFlux
  NumFlux<equationsType, dim>* numFlux;

  // Slope limiter
  SlopeLimiter<equationsType, dim>* slopeLimiter;

  // For adaptivity
  enum commonCellRelativeSize
  {
    currentMoreRefined,
    prevMoreRefined,
    equallyRefined
  };

  Adaptivity<dim>* adaptivity;

  // Indication whether mesh has been refined in this step.
  // - optimization.
  bool refined_mesh;

  const UpdateFlags update_flags;
  const UpdateFlags face_update_flags;
  const UpdateFlags neighbor_face_update_flags;
  typename DoFHandler<dim>::cell_iterator cell, prev_cell;
  FEValues<dim> fe_v_cell, fe_v_prev_cell;
  FEFaceValues<dim> fe_v_face, fe_v_prev_face;
  FESubfaceValues<dim> fe_v_subface, fe_v_prev_subface;
  FEFaceValues<dim> fe_v_prev_face_neighbor;
  FESubfaceValues<dim> fe_v_prev_subface_neighbor;
  std::vector<types::global_dof_index> dof_indices;
  std::vector<types::global_dof_index> prev_dof_indices;
  std::vector<types::global_dof_index> prev_dof_indices_neighbor;
  std::vector<std::array<double, Equations<equationsType, dim>::n_components> > normal_fluxes_old;
  std::vector<std::array<double, Equations<equationsType, dim>::n_components> > W_prev, Wplus_old, Wminus_old;
  std::vector<std::array<std::array<double, dim>, Equations<equationsType, dim>::n_components> > fluxes_old;

  std::array <unsigned short, BASIS_FN_COUNT> component_ii;
  std::array <bool, BASIS_FN_COUNT> is_primitive;
  std::array <bool, BASIS_FN_COUNT> basis_fn_is_constant;

  // This is here and not in the loop because of tests - we test by looking at the last res_norm.
  double res_norm;
  // Utility
  AssemblingUtilities<equationsType, dim> assembling_utils;
};
#endif
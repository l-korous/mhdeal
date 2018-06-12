#include "util.h"
#include "equationsMhd.h"
#include "parameters.h"
#include "initialCondition.h"
#include "boundaryConditions.h"
#include "feDivFree.h"
#include "feTaylor.h"
#include "numericalFlux.h"
#include "slopeLimiter.h"
#include "adaptivity.h"

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
    InitialCondition<equationsType, dim>& initial_condition, BoundaryConditions<equationsType, dim>& boundary_conditions);
  void run();

  // Technical matters done only once after creation.
  void setup_system();

  // Performs a single global assembly.
  void assemble_system(bool assemble_matrix = true);

  void postprocess();
  
  void set_adaptivity(Adaptivity<dim>* adaptivity);

  // Performs a single global assembly.
  void calculate_cfl_condition();

  // Performs a local assembly for all volumetric contributions on the local cell.
  void assemble_cell_term(FullMatrix<double>& cell_matrix, Vector<double>& cell_rhs, bool assemble_matrix);
  
  // Performs a local assembly for all surface contributions on the local cell.
  // i.e. face terms calculated on all faces - internal and boundary
  void assemble_face_term(const unsigned int face_no, const FEFaceValuesBase<dim> &fe_v, const FEFaceValuesBase<dim> &fe_v_neighbor, const bool external_face, const unsigned int boundary_id, Vector<double>& cell_rhs);
  
  void output_base();
  void output_results(bool use_prev_solution = false) const;
  void output_matrix(TrilinosWrappers::SparseMatrix& mat, const char* suffix) const;
  void output_vector(TrilinosWrappers::MPI::Vector& vec, const char* suffix) const;

  // Solves the assembled system
  void solve();

  void move_time_step_handle_outputs();

  void perform_reset_after_refinement();
  bool reset_after_refinement;

  // Triangulation - passed as a constructor parameter
#ifdef HAVE_MPI
  parallel::distributed::Triangulation<dim>& triangulation;
#else
  Triangulation<dim>& triangulation;
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
  IndexSet locally_owned_dofs;

  // Dofs calculated by this MPI process + all Dofs on all neighboring cells.
  IndexSet locally_relevant_dofs;

  AztecOO* solver;
  const MappingQ1<dim> mapping;
  const FESystem<dim> fe;
  DoFHandler<dim> dof_handler;
  const QGauss<dim> quadrature;
  const QGauss<dim - 1> face_quadrature;

  // Currently sought solution, the previous one, and the initial solution for newton's loop on the current time level.
  TrilinosWrappers::MPI::Vector     current_limited_solution;
  TrilinosWrappers::MPI::Vector     current_unlimited_solution;
  TrilinosWrappers::MPI::Vector     prev_solution;
  
  // The system being assembled.
  TrilinosWrappers::MPI::Vector system_rhs;
  TrilinosWrappers::SparseMatrix system_matrix;

  // Rest is technical.
  ConstraintMatrix constraints;
  MPI_Comm mpi_communicator;

  double last_output_time, time;
  int time_step;
  double cfl_time_step;
  // For CFL.
  double max_signal_speed;

  FEValuesExtractors::Vector mag;
  void precalculate_global();
  unsigned int dofs_per_cell;
  unsigned short n_quadrature_points_cell, n_quadrature_points_face;

  // NumFlux
  NumFlux<equationsType, dim>* numFlux;

  // Slope limiter
  SlopeLimiter<equationsType, dim>* slopeLimiter;

  // TODO Revise this for adaptivity (subface_flags, ...)
  const UpdateFlags update_flags;
  const UpdateFlags face_update_flags;
  const UpdateFlags neighbor_face_update_flags;
  // DOF indices both on the currently assembled element and the neighbor.
  FEValues<dim> fe_v_cell;
  FEFaceValues<dim> fe_v_face;
  FESubfaceValues<dim> fe_v_subface;
  FEFaceValues<dim> fe_v_face_neighbor;
  FESubfaceValues<dim> fe_v_subface_neighbor;
  std::vector<types::global_dof_index> dof_indices;
  std::vector<types::global_dof_index> dof_indices_neighbor;
  std::vector<std::array<double, Equations<equationsType, dim>::n_components> > Wplus_old, Wminus_old;
  std::vector<std::array<double, Equations<equationsType, dim>::n_components> > normal_fluxes_old;
  std::vector<std::array<double, Equations<equationsType, dim>::n_components> > W_prev;
  std::vector<std::array<std::array<double, dim>, Equations<equationsType, dim>::n_components> > fluxes_old;

  std::array <unsigned short, BASIS_FN_COUNT> component_ii;
  std::array <bool, BASIS_FN_COUNT> is_primitive;
  std::array <bool, BASIS_FN_COUNT> basis_fn_is_constant;

  Adaptivity<dim>* adaptivity;
};

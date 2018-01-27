#include "util.h"
#include "equationsMhd.h"
#include "parameters.h"
#include "initialCondition.h"
#include "boundaryConditions.h"
#include "DealiiExtensions.h"
#include "fe_div_free.h"
#include "fe_taylor.h"

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

private:
  // Technical matters done only once after creation.
  void setup_system();

  // Sets up initial values for current & previous solution.
  void setup_initial_solution();

  // Performs a single global assembly.
  void assemble_system(bool assemble_matrix = true);

  // Performs a single global assembly.
  void postprocess();

  // Performs a single global assembly.
  void calculate_cfl_condition();

  // Performs a local assembly for all volumetric contributions on the local cell.
  void assemble_cell_term(const FEValues<dim> &fe_v, const std::vector<types::global_dof_index>& local_dofs, FullMatrix<double>& cell_matrix, Vector<double>& cell_rhs, bool assemble_matrix);
  
  // Performs a local assembly for all surface contributions on the local cell.
  // i.e. face terms calculated on all faces - internal and boundary
  void assemble_face_term(const unsigned int face_no, const FEFaceValuesBase<dim> &fe_v, const FEFaceValuesBase<dim> &fe_v_neighbor, const std::vector<types::global_dof_index>& local_dofs,
    const std::vector<types::global_dof_index>& local_dofs_neighbor, const bool external_face, const unsigned int boundary_id, const double face_diameter, FullMatrix<double>& cell_matrix,
    Vector<double>& cell_rhs, bool assemble_matrix);
  
  // Output
  void output_results(const char* prefix = "") const;
  void output_matrix(TrilinosWrappers::SparseMatrix& mat, const char* suffix, int time_step, int newton_step = -1) const;
  void output_vector(TrilinosWrappers::MPI::Vector& vec, const char* suffix, int time_step, int newton_step = -1) const;

  // Solves the assembled system
  void solve(TrilinosWrappers::MPI::Vector &newton_update, bool reset_matrix = true);

  void move_time_step_handle_outputs();

  // Triangulation - passed as a constructor parameter
#ifdef HAVE_MPI
  parallel::distributed::Triangulation<dim>& triangulation;
#else
  Triangulation<dim>& triangulation;
#endif

  void save();
  void load();

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

  const FESystem<dim> fe;
  DoFHandler<dim> dof_handler;
  const QGauss<dim> quadrature, initial_quadrature;
  const QGauss<dim - 1> face_quadrature;

  // Currently sought solution, the previous one, and the initial solution for newton's loop on the current time level.
  TrilinosWrappers::MPI::Vector     current_limited_solution;
  TrilinosWrappers::MPI::Vector     current_unlimited_solution;
  TrilinosWrappers::MPI::Vector     prev_solution;
  TrilinosWrappers::MPI::Vector     lin_solution;
  
  // The system being assembled.
  TrilinosWrappers::MPI::Vector system_rhs;
  TrilinosWrappers::SparseMatrix system_matrix;

  // Rest is technical.
  ConditionalOStream verbose_cout;

  ConstraintMatrix constraints;

  const MappingQ1<dim> mapping;

  MPI_Comm mpi_communicator;

  bool initial_step;
  bool assemble_only_rhs;

  double last_output_time, last_snapshot_time, time;
  int time_step;
  double cfl_time_step;
  AztecOO solver;

  DealIIExtensions::PeriodicCellMap<dim> periodic_cell_map;
  FEValuesExtractors::Vector mag;
  void precalculate_global();
  unsigned int dofs_per_cell;

  // TODO Revise this for adaptivity (subface_flags, ...)
  UpdateFlags update_flags;
  UpdateFlags face_update_flags;
  UpdateFlags neighbor_face_update_flags;

  unsigned short n_quadrature_points_cell;

  // May be increased, but for linear functions, it is exactly this.
#define BASIS_FN_COUNT 100
  std::array <unsigned short, BASIS_FN_COUNT> component_ii;
  std::array <bool, BASIS_FN_COUNT> is_primitive;
  std::array <bool, BASIS_FN_COUNT> basis_fn_is_constant;

  // This is here and not in the loop because of tests - we test by looking at the last res_norm.
  public:
  double res_norm;
};

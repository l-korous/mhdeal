#include "util.h"
#include "equationsEuler.h"
#include "equationsMhd.h"
#include "parameters.h"
#include "initialCondition.h"
#include "boundaryConditions.h"

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
  void setup_system();

  void process_initial_condition();

  void assemble_system();
  void assemble_cell_term(const FEValues<dim> &fe_v, const std::vector<types::global_dof_index>& local_dofs, FullMatrix<double>& cell_matrix, Vector<double>& cell_rhs);
  void assemble_face_term(const unsigned int face_no, const FEFaceValuesBase<dim> &fe_v, const FEFaceValuesBase<dim> &fe_v_neighbor, const std::vector<types::global_dof_index>& local_dofs,
    const std::vector<types::global_dof_index>& local_dofs_neighbor, const bool external_face, const unsigned int boundary_id, const double face_diameter, FullMatrix<double>& cell_matrix,
    Vector<double>& cell_rhs, FullMatrix<double>& cell_matrix_neighbor, Vector<double>& cell_rhs_neighbor);

#ifdef HAVE_MPI
  MPI_Comm mpi_communicator;
  parallel::distributed::Triangulation<dim>& triangulation;
#else
  Triangulation<dim>& triangulation;
#endif

  void output_results() const;

  IndexSet locally_owned_dofs;

  IndexSet locally_relevant_dofs;

  ConstraintMatrix constraints;

  const MappingQ1<dim> mapping;

  const FESystem<dim> fe;

  DoFHandler<dim> dof_handler;

  const QGauss<dim> quadrature;

  const QGauss<dim - 1> face_quadrature;

#ifdef HAVE_MPI
  dealii::LinearAlgebraTrilinos::MPI::Vector locally_relevant_solution;

  TrilinosWrappers::MPI::Vector     current_solution;

  TrilinosWrappers::MPI::Vector     old_solution;

  TrilinosWrappers::MPI::Vector newton_initial_guess;

  dealii::LinearAlgebraTrilinos::MPI::Vector system_rhs;

  dealii::LinearAlgebraTrilinos::MPI::SparseMatrix system_matrix;

  void solve(TrilinosWrappers::MPI::Vector &newton_update);
#else
  Vector<double> old_solution;

  Vector<double> current_solution;

  Vector<double> newton_initial_guess;

  Vector<double> system_rhs;

  TrilinosWrappers::SparseMatrix system_matrix;

  void solve(Vector<double> &newton_update);
#endif

  Equations<equationsType, dim>& equations;

  Parameters<dim>& parameters;

  InitialCondition<equationsType, dim>& initial_condition;

  BoundaryConditions<equationsType, dim>& boundary_conditions;

  ConditionalOStream              verbose_cout;
};
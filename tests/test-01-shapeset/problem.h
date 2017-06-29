#include "util.h"
#include "equations.h"
#include "parameters.h"

// Class that accepts all input from the user, provides interface for output, etc.
// Should not be changed.
template <int dim>
class Problem
{
public:
  Problem(Parameters<dim>& parameters, Equations<dim>& equations,
#ifdef HAVE_MPI
    parallel::distributed::Triangulation<dim>& triangulation
#else
    Triangulation<dim>& triangulation
#endif
  );
  void run();

private:
  // Technical matters done only once after creation.
  void setup_system();

  // Performs a single global assembly.
  void assemble_system();

  // Performs a local assembly for all volumetric contributions on the local cell.
  void assemble_cell_term(const FEValues<dim> &fe_v, const std::vector<types::global_dof_index>& local_dofs, FullMatrix<double>& cell_matrix, Vector<double>& cell_rhs);
  
  // Output
  void output_results() const;

  // Solves the assembled system
  void solve(TrilinosWrappers::MPI::Vector &newton_update);

  // Triangulation - passed as a constructor parameter
#ifdef HAVE_MPI
  parallel::distributed::Triangulation<dim>& triangulation;
#else
  Triangulation<dim>& triangulation;
#endif

  // Parameters - passed as a constructor parameter
  Parameters<dim>& parameters;

  // Dofs calculated by this MPI process.
  IndexSet locally_owned_dofs;

  // Dofs calculated by this MPI process + all Dofs on all neighboring cells.
  IndexSet locally_relevant_dofs;

  const FESystem<dim> fe;
  DoFHandler<dim> dof_handler;
  const QGauss<dim> quadrature;

  TrilinosWrappers::MPI::Vector solution;

  // Equations - passed as a constructor parameter
  Equations<dim>& equations;

  // The system being assembled.
  TrilinosWrappers::MPI::Vector system_rhs;
  TrilinosWrappers::SparseMatrix system_matrix;

  // Rest is technical.
  ConditionalOStream              verbose_cout;

  ConstraintMatrix constraints;

  const MappingQ1<dim> mapping;

  MPI_Comm mpi_communicator;
};

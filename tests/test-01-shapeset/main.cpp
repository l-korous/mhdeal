#include "util.h"
#include "problem.h"

// Dimension of the problem - passed as a template parameter to pretty much every class.
#define DIMENSION 3
// Type of must be from the enumeration EquationsType defined in equations.h.
#define EQUATIONS EquationsTypeMhd
// Whether or not we want to delete all previous outputs (vtk, matrix outputs, solution output) on start.
#define DELETE_OUTPUTS_ON_START

int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, dealii::numbers::invalid_unsigned_int);
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);

  try
  {
    // Declaration of triangulation. The triangulation is not initialized here, but rather in the constructor of Parameters class.
#ifdef HAVE_MPI
    parallel::distributed::Triangulation<DIMENSION> triangulation(mpi_communicator, typename Triangulation<DIMENSION>::MeshSmoothing(Triangulation<DIMENSION>::smoothing_on_refinement | Triangulation<DIMENSION>::smoothing_on_coarsening));
#else
    Triangulation<DIMENSION> triangulation;
#endif    
    
    // Initialization of parameters. See parameters.h for description of the individual parameters, set up values in parameters.cpp
    Parameters<DIMENSION> parameters(triangulation);
    // Set up equations - see equations.h, equations.h
    Equations<DIMENSION> equations;
    // Put together the problem.
    Problem<DIMENSION> problem(parameters, equations, triangulation);
    // Run the problem - entire transient problem.
    problem.run();
  }
  catch (std::exception &exc)
  {
    std::cerr << std::endl << std::endl
      << "----------------------------------------------------"
      << std::endl;
    std::cerr << "Exception on processing: " << std::endl
      << exc.what() << std::endl
      << "Aborting!" << std::endl
      << "----------------------------------------------------"
      << std::endl;
    return 1;
  }
  catch (...)
  {
    std::cerr << std::endl << std::endl
      << "----------------------------------------------------"
      << std::endl;
    std::cerr << "Unknown exception!" << std::endl
      << "Aborting!" << std::endl
      << "----------------------------------------------------"
      << std::endl;
    return 1;
  };
  
  return 0;
}

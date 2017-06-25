#include "util.h"
#include "problem.h"

// Dimension of the problem - passed as a template parameter to pretty much every class.
#define DIMENSION 3
// Type of equations, must be from the enumeration EquationsType defined in equations.h.
#define EQUATIONS EquationsTypeMhd
// Whether or not we want to delete all previous outputs (vtk, matrix outputs, solution output) on start.
#define DELETE_OUTPUTS_ON_START

int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, dealii::numbers::invalid_unsigned_int);
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);
  InitialCondition<EQUATIONS, DIMENSION> *initial_condition;

  try
  {
    // The main process will optionally delete outputs.
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
#ifdef DELETE_OUTPUTS_ON_START
#ifdef _MSC_VER
      system("del *.visit");
      system("del *.vtk");
      system("del *.vtu");
      system("del *.pvtu");
      system("del *.newton_update");
      system("del *.current_solution");
      system("del *.matrix");
      system("del *.rhs");
#else
      system("rm *.visit");
      system("rm *.vtk");
      system("rm *.vtu");
      system("rm *.pvtu");
      system("rm *.newton_update");
      system("rm *.current_solution");
      system("rm *.matrix");
      system("rm *.rhs");
#endif
#endif
    }

    // Declaration of triangulation. The triangulation is not initialized here, but rather in the constructor of Parameters class.
#ifdef HAVE_MPI
    parallel::distributed::Triangulation<DIMENSION> triangulation(mpi_communicator, typename Triangulation<DIMENSION>::MeshSmoothing(Triangulation<DIMENSION>::smoothing_on_refinement | Triangulation<DIMENSION>::smoothing_on_coarsening));
#else
    Triangulation<DIMENSION> triangulation;
#endif    
    
    // Initialization of parameters. See parameters.h for description of the individual parameters, set up values in parameters.cpp
    Parameters<DIMENSION> parameters(triangulation);
    // Set up of initial condition. See initialCondition.h for description of methods, set up the specific function in initialCondition.cpp
    switch(parameters.initCond){
      case 0:
        initial_condition = new MHDBlastIC<EQUATIONS, DIMENSION>(parameters);
        break;
      case 1:
        initial_condition = new TitovDemoulinIC<EQUATIONS, DIMENSION>(parameters);
        break;
      case 2:
        initial_condition = new TaylorBasisTestIC<EQUATIONS, DIMENSION>(parameters);
        break;
    }
    // Set up of boundary condition. See boundaryCondition.h for description of methods, set up the specific function in boundaryCondition.cpp
    BoundaryConditions<EQUATIONS, DIMENSION> boundary_conditions;
    // Set up equations - see equations.h, equationsMhd.h
    Equations<EQUATIONS, DIMENSION> equations(parameters);
    // Put together the problem.
    Problem<EQUATIONS, DIMENSION> problem(parameters, equations, triangulation, *initial_condition, boundary_conditions);
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
    delete initial_condition;
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
    delete initial_condition;
    return 1;
  };

  delete initial_condition;
  
  return 0;
}

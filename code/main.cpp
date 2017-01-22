#include "util.h"
#include "problem.h"

#define DIMENSION 3
#define EQUATIONS EquationsTypeMhd
#define DELETE_VTK_ON_START

int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, dealii::numbers::invalid_unsigned_int);
#ifdef HAVE_MPI
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);
#endif

  try
  {
#ifdef HAVE_MPI
  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
#endif
  {
#ifdef DELETE_VTK_ON_START
#ifdef _MSC_VER
    system("del *.vtk");
    system("del *.vtu");
    system("del *.pvtu");
#else
    system("rm *.vtk");
    system("rm *.vtu");
    system("rm *.pvtu");
#endif
#endif
  }

#ifdef HAVE_MPI
    parallel::distributed::Triangulation<DIMENSION> triangulation(mpi_communicator, typename Triangulation<DIMENSION>::MeshSmoothing(Triangulation<DIMENSION>::smoothing_on_refinement | Triangulation<DIMENSION>::smoothing_on_coarsening));
#else
    Triangulation<DIMENSION> triangulation;
#endif    

    Parameters<DIMENSION> parameters(triangulation);

    InitialCondition<EQUATIONS, DIMENSION> initial_condition(parameters);
    BoundaryConditions<EQUATIONS, DIMENSION> boundary_conditions;

    Equations<EQUATIONS, DIMENSION> equations(parameters);
    Problem<EQUATIONS, DIMENSION> problem(parameters, equations, triangulation, initial_condition, boundary_conditions);
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

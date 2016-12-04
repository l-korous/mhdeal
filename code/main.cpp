#include "util.h"
#include "problem.h"

#define DIMENSION 3
#define EQUATIONS EquationsTypeMhd
#define DELETE_VTK_ON_START

int main(int argc, char *argv[])
{
#ifdef DELETE_VTK_ON_START

#ifdef _MSC_VER
    system("del *.vtk");
#else
    system("rm *.vtk");
#endif

#endif

  try
  {
    using namespace dealii;

    if (argc != 2)
    {
      std::cout << "Usage:" << argv[0] << " input_file" << std::endl;
      std::exit(1);
    }

    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, dealii::numbers::invalid_unsigned_int);

    Triangulation<DIMENSION> triangulation;
    
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

#include "util.h"
#include "equations.h"
#include "flux.h"
#include "parameters.h"
#include "problem.h"

int main(int argc, char *argv[])
{
  try
  {
    using namespace dealii;

    if (argc != 2)
    {
      std::cout << "Usage:" << argv[0] << " input_file" << std::endl;
      std::exit(1);
    }

    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, dealii::numbers::invalid_unsigned_int);

    Parameters<2> parameters;
    Problem<2> problem(parameters);
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

#include "definitions.h"
#include "mhdSolver.h"

using namespace dealii;

int main()
{
  if (DELETE_VTK)
  {
#ifdef _MSC_VER
    system("del *.vtk");
#else
    system("rm *.vtk");
#endif
  }

  try
  {
    MHDSolver feProblem;
    feProblem.run();
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

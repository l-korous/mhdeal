#include "util.h"
#include "equations.h"
#include "parameters.h"
#include "problem.h"

#define DIMENSION 2

template<int dim>
void load_mesh(Triangulation<dim>& triangulation, std::string mesh_filename);

template<>
void load_mesh<2>(Triangulation<2>& triangulation, std::string mesh_filename)
{
  GridIn<2> grid_in;
  grid_in.attach_triangulation(triangulation);

  std::ifstream input_file(mesh_filename.c_str());
  Assert(input_file, ExcFileNotOpen(mesh_filename.c_str()));

  grid_in.read_ucd(input_file);
}

template<>
void load_mesh<3>(Triangulation<3>& triangulation, std::string mesh_filename)
{
  Triangulation<2> tempTriangulation;
  GridIn<2> grid_in;

  grid_in.attach_triangulation(tempTriangulation);

  std::ifstream input_file(mesh_filename.c_str());
  Assert(input_file, ExcFileNotOpen(mesh_filename.c_str()));

  grid_in.read_ucd(input_file);

  GridGenerator::extrude_triangulation(tempTriangulation, 3, 2.0, triangulation);
}

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

    Triangulation<DIMENSION> triangulation;

    Parameters<DIMENSION> parameters;
    load_mesh<DIMENSION>(triangulation, parameters.mesh_filename);
    InitialCondition<EquationsTypeEuler, DIMENSION> initial_condition;
    BoundaryConditions<EquationsTypeEuler, DIMENSION> boundary_conditions;
    Problem<EquationsTypeEuler, DIMENSION> problem(parameters, triangulation, initial_condition, boundary_conditions);
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

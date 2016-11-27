#include "util.h"
#include "problem.h"

#define DIMENSION 3
#define EQUATIONS EquationsTypeEuler

template<int dim>
void load_mesh(Triangulation<dim>& triangulation, Parameters<DIMENSION>& parameters);

template<>
void load_mesh<2>(Triangulation<2>& triangulation, Parameters<DIMENSION>& parameters)
{
  GridIn<2> grid_in;
  grid_in.attach_triangulation(triangulation);

  std::ifstream input_file(parameters.mesh_filename.c_str());
  Assert(input_file, ExcFileNotOpen(parameters.mesh_filename.c_str()));

  grid_in.read_ucd(input_file);
}

template<>
void load_mesh<3>(Triangulation<3>& triangulation, Parameters<DIMENSION>& parameters)
{
  Triangulation<2> tempTriangulation;
  GridIn<2> grid_in;

  grid_in.attach_triangulation(tempTriangulation);

  std::ifstream input_file(parameters.mesh_filename.c_str());
  Assert(input_file, ExcFileNotOpen(parameters.mesh_filename.c_str()));

  grid_in.read_ucd(input_file);

  GridGenerator::extrude_triangulation(tempTriangulation, parameters.MeshSlicesInZDirection, 0.5, triangulation);

  GridOut grid_out;
  std::ofstream output_file("extrudedMesh.vtk");
  grid_out.write_vtk(triangulation, output_file);
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
    load_mesh<DIMENSION>(triangulation, parameters);
    Equations<EQUATIONS, DIMENSION> equations(parameters);
    InitialCondition<EQUATIONS, DIMENSION> initial_condition;
    BoundaryConditions<EQUATIONS, DIMENSION> boundary_conditions;
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

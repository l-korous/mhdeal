#include "util.h"
#include "problem.h"
#include "equationsMhd.h"
#include "parameters.h"

// Dimension of the problem - passed as a template parameter to pretty much every class.
#define DIMENSION 3
// Type of equations, must be from the enumeration EquationsType defined in equations.h.
#define EQUATIONS EquationsTypeMhd

#ifdef HAVE_MPI
void set_triangulation(parallel::distributed::Triangulation<DIMENSION>& triangulation, Parameters<DIMENSION>& parameters)
#else
void set_triangulation(Triangulation<DIMENSION>& triangulation, Parameters<DIMENSION>& parameters)
#endif
{
  GridGenerator::subdivided_hyper_rectangle(triangulation, parameters.refinements, parameters.corner_a, parameters.corner_b, true);

  std::vector<DealIIExtensions::FacePair<DIMENSION> > matched_pairs;
  for (std::vector<std::array<int, 3> >::const_iterator it = parameters.periodic_boundaries.begin(); it != parameters.periodic_boundaries.end(); it++)
    dealii::GridTools::collect_periodic_faces(triangulation, (*it)[0], (*it)[1], (*it)[2], matched_pairs);
  triangulation.add_periodicity(matched_pairs);
}

void set_parameters(Parameters<DIMENSION>& parameters)
{
  parameters.corner_a = Point<DIMENSION>(0.0, 0.0, 0.);
  parameters.corner_b = Point<DIMENSION>(0.2, 0.2, 0.001);
  parameters.refinements = { 30, 30, 1 };
  parameters.limit = false;
  parameters.slope_limiter = parameters.barthJespersen;
  parameters.use_div_free_space_for_B = true;
  //parameters.periodic_boundaries = { { 0, 1, 0 },{ 2, 3, 1 } };
  parameters.num_flux_type = Parameters<DIMENSION>::hlld;
  parameters.lax_friedrich_stabilization_value = 0.5;
  parameters.cfl_coefficient = .01;
  parameters.quadrature_order = 5;
  parameters.polynomial_order_dg = 1;

  parameters.use_iterative_improvement = false;
  parameters.patches = 2;
  parameters.output_step = -1.e-3;

  parameters.debug = false;
  parameters.gas_gamma = 1.4;
}

int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, dealii::numbers::invalid_unsigned_int);
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);

  try
  {
    // The main process will optionally delete outputs.
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
#ifdef _MSC_VER
      system("del *.visit *.vtk *.vtu *.pvtu *.newton_update *.current_solution *.matrix *.rhs");
#else
      system("rm *.visit *.vtk *vtu *.newton_update *.current_solution *.matrix *.rhs");
#endif
    }

    // Initialization of parameters. See parameters.h for description of the individual parameters
    Parameters<DIMENSION> parameters;
    set_parameters(parameters);

    // Declaration of triangulation. The triangulation is not initialized here, but rather in the constructor of Parameters class.
#ifdef HAVE_MPI
    parallel::distributed::Triangulation<DIMENSION> triangulation(mpi_communicator, typename Triangulation<DIMENSION>::MeshSmoothing(Triangulation<DIMENSION>::smoothing_on_refinement | Triangulation<DIMENSION>::smoothing_on_coarsening));
#else
    Triangulation<DIMENSION> triangulation;
#endif    
    set_triangulation(triangulation, parameters);

    MHDBlastIC<EQUATIONS, DIMENSION> initial_condition(parameters);
    // Set up of boundary condition. See boundaryCondition.h for description of methods, set up the specific function in boundaryCondition.cpp
    BoundaryConditions<EQUATIONS, DIMENSION> boundary_conditions(parameters);
    // Set up equations - see equations.h, equationsMhd.h
    Equations<EQUATIONS, DIMENSION> equations;
    // Put together the problem.
    Problem<EQUATIONS, DIMENSION> problem(parameters, equations, triangulation, initial_condition, boundary_conditions);
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

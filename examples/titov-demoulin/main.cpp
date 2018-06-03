#include "util.h"
#include "problem.h"
#include "equationsMhd.h"
#include "parameters.h"
#include "initialConditionTD.h"
#include "adaptivityTD.h"

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

  std::vector<dealii::GridTools::PeriodicFacePair< dealii::TriaIterator<dealii::CellAccessor<DIMENSION> > > > matched_pairs;
  for (std::vector<std::array<int, 3> >::const_iterator it = parameters.periodic_boundaries.begin(); it != parameters.periodic_boundaries.end(); it++)
    dealii::GridTools::collect_periodic_faces(triangulation, (*it)[0], (*it)[1], (*it)[2], matched_pairs);
  triangulation.add_periodicity(matched_pairs);
}

// Parameters that are specific for this example.
int max_cells;
int refine_every_nth_time_step;
int perform_n_initial_refinements;
double refine_threshold;
double coarsen_threshold;

void set_parameters(Parameters<DIMENSION>& parameters)  
{
  parameters.slope_limiter = parameters.vertexBased;
  parameters.corner_a = Point<DIMENSION>(-5., -10., 0.);
  parameters.corner_b = Point<DIMENSION>(5., 10., 10.);
  parameters.refinements = { 10, 20, 10 };
  parameters.limit = true;
  parameters.use_div_free_space_for_B = true;
  parameters.num_flux_type = Parameters<DIMENSION>::hlld;
  parameters.lax_friedrich_stabilization_value = 0.5;
  parameters.cfl_coefficient = .001;
  parameters.start_limiting_at = -.05;
  parameters.quadrature_order = 1;
  parameters.polynomial_order_dg = 0;
  parameters.patches = 0;
  parameters.output_step = -1.e-2;
  parameters.final_time = 10.;

  max_cells = 3500;
  refine_every_nth_time_step = 25;
  perform_n_initial_refinements = 25;
  refine_threshold = 0.10;
  coarsen_threshold = 0.15;
}

int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, dealii::numbers::invalid_unsigned_int);
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);

  try
  {
    // Initialization of parameters. See parameters.h for description of the individual parameters
    Parameters<DIMENSION> parameters;
    set_parameters(parameters);
    parameters.delete_old_outputs(mpi_communicator);

    // Declaration of triangulation. The triangulation is not initialized here, but rather in the constructor of Parameters class.
#ifdef HAVE_MPI
    parallel::distributed::Triangulation<DIMENSION> triangulation(mpi_communicator, typename dealii::Triangulation<DIMENSION>::MeshSmoothing(Triangulation<DIMENSION>::none), parallel::distributed::Triangulation<DIMENSION>::no_automatic_repartitioning);
#else
    Triangulation<DIMENSION> triangulation;
#endif    
    set_triangulation(triangulation, parameters);

    InitialConditionTitovDemoulin<EQUATIONS, DIMENSION> initial_condition(parameters);
    // Set up of boundary condition. See boundaryCondition.h for description of methods, set up the specific function in boundaryCondition.cpp
    BoundaryConditions<EQUATIONS, DIMENSION> boundary_conditions(parameters);
    // Set up equations - see equations.h, equationsMhd.h
    Equations<EQUATIONS, DIMENSION> equations;
    // Adaptivity
    AdaptivityTD<DIMENSION> adaptivity(parameters, mpi_communicator, max_cells, refine_every_nth_time_step, perform_n_initial_refinements, refine_threshold, coarsen_threshold);
    // Put together the problem.
    Problem<EQUATIONS, DIMENSION> problem(parameters, equations, triangulation, initial_condition, boundary_conditions);
    // Run the problem - entire transient problem.
    problem.set_adaptivity(&adaptivity);
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

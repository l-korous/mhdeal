#include "util.h"
#include "problem.h"
#include "equationsMhd.h"
#include "parameters.h"
#include "parametersTD.h"
#include "initialConditionTD.h"
#include "boundaryConditionTD.h"
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
}

void set_parameters(Parameters<DIMENSION>& parameters, TitovDemoulinParameters& td_parameters)
{
  parameters.slope_limiter = parameters.vertexBased;
  parameters.corner_a = Point<DIMENSION>(-5., -10., 0.);
  parameters.corner_b = Point<DIMENSION>(5., 10., 10.);
  parameters.refinements = { 50, 100, 50 };
  parameters.limit = false;
  parameters.limitB = false;
  parameters.use_div_free_space_for_B = false;
  parameters.num_flux_type = Parameters<DIMENSION>::hlld;
  parameters.lax_friedrich_stabilization_value = 0.5;
  parameters.cfl_coefficient = .01;
  parameters.start_limiting_at = -1e-6;
  parameters.quadrature_order = 1;
  parameters.polynomial_order_dg = 0;
  parameters.patches = 0;
  parameters.output_step = -1.e-2;
  parameters.final_time = 1.e-8;

  parameters.max_cells = 2500;
  parameters.refine_every_nth_time_step = 25;
  parameters.perform_n_initial_refinements = 25;
  parameters.refine_threshold = 0.5;
  parameters.coarsen_threshold = 0.2;
  parameters.volume_factor = 4;
  parameters.time_interval_max_cells_multiplicator = 0.;

  // plasma beta
  td_parameters.beta = 0.05;

  // coronal height scale
  td_parameters.L_G = 0.;

  // Density
  td_parameters.rho_0 = 1.;

  // Torus winding number
  td_parameters.N_t = 5.;

  // Torus major radius
  td_parameters.R = 3.0;

  // Torus minor radius
  td_parameters.r = 0.5;

  // Magnetic charge separation distance
  td_parameters.L = 2.0;

  // Geometrical factor
  td_parameters.d = 0.5;
  
  // The coronal/prominence temperature ratio
  td_parameters.Tc2Tp = 1.;

  td_parameters.omega_0 = 0.3;

  td_parameters.t_drive = 2.0;

  td_parameters.t_ramp = 1.0;
}

int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, dealii::numbers::invalid_unsigned_int);
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);

  try
  {
    // Initialization of parameters. See parameters.h for description of the individual parameters
    Parameters<DIMENSION> parameters;
    TitovDemoulinParameters td_parameters;
    set_parameters(parameters, td_parameters);
    parameters.delete_old_outputs(mpi_communicator);

    // Declaration of triangulation. The triangulation is not initialized here, but rather in the constructor of Parameters class.
#ifdef HAVE_MPI
    parallel::distributed::Triangulation<DIMENSION> triangulation(mpi_communicator, typename dealii::Triangulation<DIMENSION>::MeshSmoothing(Triangulation<DIMENSION>::limit_level_difference_at_vertices));
#else
    Triangulation<DIMENSION> triangulation(Triangulation<DIMENSION>::limit_level_difference_at_vertices);
#endif    
    set_triangulation(triangulation, parameters);

    InitialConditionTitovDemoulin<EQUATIONS, DIMENSION> initial_condition(parameters, td_parameters);
    // Set up of boundary condition. See boundaryCondition.h for description of methods, set up the specific function in boundaryCondition.cpp
    BoundaryConditionTDInitialState<DIMENSION> boundary_conditions(parameters, td_parameters);
    // Set up equations - see equations.h, equationsMhd.h
    Equations<EQUATIONS, DIMENSION> equations;
    // Adaptivity
    AdaptivityTD<DIMENSION> adaptivity(parameters, mpi_communicator);
    // Put together the problem.
    Problem<EQUATIONS, DIMENSION> problem(parameters, equations, triangulation, initial_condition, boundary_conditions);
    // Set adaptivity
    // problem.set_adaptivity(&adaptivity);
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

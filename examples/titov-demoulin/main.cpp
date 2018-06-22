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

  std::vector<dealii::GridTools::PeriodicFacePair< dealii::TriaIterator<dealii::CellAccessor<DIMENSION> > > > matched_pairs;
  for (std::vector<std::array<int, 3> >::const_iterator it = parameters.periodic_boundaries.begin(); it != parameters.periodic_boundaries.end(); it++)
    dealii::GridTools::collect_periodic_faces(triangulation, (*it)[0], (*it)[1], (*it)[2], matched_pairs);
  triangulation.add_periodicity(matched_pairs);
}

void set_parameters(Parameters<DIMENSION>& parameters, TitovDemoulinParameters& td_parameters)
{
  parameters.slope_limiter = parameters.vertexBased;
  parameters.corner_a = Point<DIMENSION>(-12., -12., 0.);
  parameters.corner_b = Point<DIMENSION>(12., 12., 20.);
  parameters.refinements = { 50, 60, 50 };
  parameters.limit = false;
  parameters.use_div_free_space_for_B = false;
  parameters.num_flux_type = Parameters<DIMENSION>::hlld;
  parameters.lax_friedrich_stabilization_value = 0.5;
  parameters.cfl_coefficient = .01;
  parameters.start_limiting_at = -.05;
  parameters.quadrature_order = 1;
  parameters.polynomial_order_dg = 0;
  parameters.patches = 0;
  parameters.output_step = -1.e-2;
  parameters.final_time = 20.;

  parameters.max_cells = 3500;
  parameters.refine_every_nth_time_step = 25;
  parameters.perform_n_initial_refinements = 25;
  parameters.refine_threshold = 0.10;
  parameters.coarsen_threshold = 0.15;

  // plasma beta
  td_parameters.beta = 0.05;

  // coronal height scale
  td_parameters.L_G = 0.;
  // Shouldn't be 20? According to page 10/456 of https://github.com/l-korous/doctoral-thesis/blob/master/_reference/Modeling%20of%20H%CE%B1%20Eruptive%20Events%20Observed%20at%20the%20Solar.pdf

  // Torus winding number
  td_parameters.N_t = 0.3;

  // Torus major radius
  td_parameters.R= 8.0;

  // Magnetic charge separation distance
  td_parameters.L = 2.0;

  // Geometrical factor
  td_parameters.d = 2.0;
  
  // "Helicity" factor inside tho loop (used later in B_theta_internal calcs)
  td_parameters.H = 2.0 * (td_parameters.N_t * td_parameters.N_t) / (td_parameters.R * td_parameters.R);

  // The coronal/prominence temperature ratio
  td_parameters.Tc2Tp = 1.0;
  // Shouldn't be 10? According to page 10/456 of https://github.com/l-korous/doctoral-thesis/blob/master/_reference/Modeling%20of%20H%CE%B1%20Eruptive%20Events%20Observed%20at%20the%20Solar.pdf

  td_parameters.omega_0 = 0.3;

  td_parameters.t_drive = 2.0;

  td_parameters.t_ramp = 1.0;

  // density jump half-width...
  double t_rho = 0.12;
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
    parallel::distributed::Triangulation<DIMENSION> triangulation(mpi_communicator, typename dealii::Triangulation<DIMENSION>::MeshSmoothing(Triangulation<DIMENSION>::none), parallel::distributed::Triangulation<DIMENSION>::no_automatic_repartitioning);
#else
    Triangulation<DIMENSION> triangulation;
#endif    
    set_triangulation(triangulation, parameters);

    InitialConditionTitovDemoulin<EQUATIONS, DIMENSION> initial_condition(parameters, td_parameters);
    // Set up of boundary condition. See boundaryCondition.h for description of methods, set up the specific function in boundaryCondition.cpp
    BoundaryConditionTDFree<DIMENSION> boundary_conditions(parameters, td_parameters);
    // Set up equations - see equations.h, equationsMhd.h
    Equations<EQUATIONS, DIMENSION> equations;
    // Adaptivity
    AdaptivityTD<DIMENSION> adaptivity(parameters, mpi_communicator);
    // Put together the problem.
    Problem<EQUATIONS, DIMENSION> problem(parameters, equations, triangulation, initial_condition, boundary_conditions);
    // Set adaptivity
    //problem.set_adaptivity(&adaptivity);
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

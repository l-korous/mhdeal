#include "parameters.h"
#include "equationsEuler.h"
#include "equationsMhd.h"

template<int dim>
void load_slide_mesh(Triangulation<dim>& triangulation, Parameters<dim>& parameters);

template<>
void load_slide_mesh<2>(Triangulation<2>& triangulation, Parameters<2>& parameters)
{
  GridIn<2> grid_in;
  grid_in.attach_triangulation(triangulation);

  std::ifstream input_file(parameters.mesh_filename.c_str());
  Assert(input_file, ExcFileNotOpen(parameters.mesh_filename.c_str()));

  grid_in.read_ucd(input_file);
}

template<>
void load_slide_mesh<3>(Triangulation<3>& triangulation, Parameters<3>& parameters)
{
  Triangulation<2> tempTriangulation;
  GridIn<2> grid_in;

  grid_in.attach_triangulation(tempTriangulation);

  std::ifstream input_file(parameters.mesh_filename.c_str());
  Assert(input_file, ExcFileNotOpen(parameters.mesh_filename.c_str()));

  grid_in.read_ucd(input_file);

  GridGenerator::extrude_triangulation(tempTriangulation, parameters.mesh_extrusion_slices, 0.5, triangulation);
}

template<int dim>
void load_cube_mesh(Triangulation<dim>& triangulation, Parameters<dim>& parameters)
{
  GridGenerator::subdivided_hyper_rectangle(triangulation, parameters.refinements, parameters.corner_a, parameters.corner_b);
}

template <int dim>
#ifdef HAVE_MPI
Parameters<dim>::Parameters(parallel::distributed::Triangulation<dim> &triangulation, Triangulation<dim> &sharedTriangulationForInitialCondition)
#else
Parameters<dim>::Parameters(Triangulation<dim> &triangulation)
#endif
{
  this->corner_a = Point<dim>(-.5, -.75, -.01);
  this->corner_b = Point<dim>(.5, .75, .01);
  this->refinements = { 100, 100, 1 };

#ifdef HAVE_MPI
  load_cube_mesh<dim>(sharedTriangulationForInitialCondition, *this);
#endif
  load_cube_mesh<dim>(triangulation, *this);

  this->final_time = 10.;
  this->time_step = 1e-6;
  this->theta = 0.0;
  this->time_step_after_initialization = 1e-4;
  this->theta_after_initialization = 0.;
  this->initialization_time = 0.;

  this->output = OutputType::quiet_solver;
  this->output_matrix = false;

  this->linear_residual = 1e-10;
  this->max_iterations = 10000;
  this->ilut_fill = 1.5;
  this->ilut_drop = 1e-6;
  this->ilut_atol = 1e-6;
  this->ilut_rtol = 1.0;
  this->newton_damping = 1.;

  this->gas_gamma = 1.4;

  this->polynomial_order_dg = 0;
  this->polynomial_order_hdiv = 1;
  this->max_nonlinear_iterations = 30;
  this->nonlinear_residual_norm_threshold = 1e-10;

  this->output_step = 1e-3;

  this->num_flux_type = hlld;
  this->lax_friedrich_stabilization_value = 1.;

  this->is_stationary = false;

  this->needs_gradients = false;
}

template class Parameters<2>;
template class Parameters<3>;
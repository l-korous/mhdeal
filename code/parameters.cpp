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

  GridGenerator::extrude_triangulation(tempTriangulation, parameters.MeshSlicesInZDirection, 0.5, triangulation);
}

template<int dim>
void load_cube_mesh(Triangulation<dim>& triangulation, Parameters<dim>& parameters);

template<>
void load_cube_mesh<2>(Triangulation<2>& triangulation, Parameters<2>& parameters)
{
  GridGenerator::hyper_cube(triangulation, parameters.cube_left, parameters.cube_right);
  triangulation.refine_global(parameters.uniform_refinements);
}

template<>
void load_cube_mesh<3>(Triangulation<3>& triangulation, Parameters<3>& parameters)
{
  GridGenerator::subdivided_hyper_rectangle(triangulation, std::vector<unsigned int>({ parameters.uniform_refinements, parameters.uniform_refinements, 1 }), Point<3>(parameters.cube_left, parameters.cube_left, -.01), Point<3>(parameters.cube_right, parameters.cube_right, 0.01));
}


template <int dim>
Parameters<dim>::Parameters(Triangulation<dim> &triangulation)
{
  this->mesh_filename = "slide.inp";
  this->MeshSlicesInZDirection = 2;
  this->cube_left = -.5;
  this->cube_right = .5;
  this->uniform_refinements = 50;
  load_cube_mesh<dim>(triangulation, *this);

  this->final_time = 10.;
  this->time_step = .0001;
  this->theta = 0.0;
  this->time_step_after_initialization = .0001;
  this->theta_after_initialization = .5;
  this->initialization_time = 0.;

  this->output = OutputType::quiet_solver;
  this->solver = SolverType::direct;
  this->linear_residual = 1e-10;
  this->max_iterations = 300;
  this->ilut_fill = 1.5;
  this->ilut_drop = 1e-6;
  this->ilut_atol = 1e-6;
  this->ilut_rtol = 1.0;

  this->gas_gamma = 1.4;

  this->polynomial_order = 0;
  this->max_nonlinear_iterations = 30;
  this->nonlinear_residual_norm_threshold = 1e-10;

  this->output_step = this->time_step;

  this->num_flux_type = hlld;
  this->lax_friedrich_stabilization_value = 1.;

  this->is_stationary = false;
}

template class Parameters<2>;
template class Parameters<3>;
#include "parameters.h"
#include "equationsEuler.h"
#include "equationsMhd.h"

template <int dim>
Parameters<dim>::Parameters()
{
  this->mesh_filename = "slide.inp";
  this->MeshSlicesInZDirection = 2;
  this->final_time = 10.;
  this->time_step = .01;
  this->theta = 0.5;

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
  this->nonlinear_residual_norm_threshold = 1e-8;

  this->output_step = 0.01;

  this->stabilization_kind = StabilizationKind::constant_stabilization;
  this->stabilization_value = 1.;

  this->is_stationary = false;
}

template class Parameters<2>;
template class Parameters<3>;
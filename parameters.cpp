#include "parameters.h"

template <int dim>
Parameters<dim>::Parameters() {
  // Just a large enough number here.
  this->dof_threshold = 100000000;
  this->start_limiting_at = 0.;
  this->gas_gamma = 5. / 3.;
  this->limit_in_nonlin_loop = false;
  this->automatic_damping = false;
  this->automatic_cfl = false;
  this->initial_and_max_newton_damping = 1.;
  this->decrease_factor = .9;
  this->increase_factor = 1. / this->decrease_factor;
  this->stagnation_coefficient = 1.e-2;
  this->bad_step_coefficient = 2.;
  this->limit = true;
  this->slope_limiter = vertexBased;
  this->output_file_prefix = "";
  this->lax_friedrich_stabilization_value = .5;
  this->snapshot_step = 1.;
  this->time_step = 1.e-6;

  this->debug = false;
  this->output_matrix = false;
  this->output = quiet_solver;
  this->output_rhs = false;
  this->output_solution = false;

  this->solver = gmres;
  this->linear_residual = 1e-10;
  this->max_iterations = 10000;
  this->ilut_fill = 1.5;
  this->ilut_drop = 1e-6;
  this->ilut_atol = 1e-6;
  this->ilut_rtol = 1.0;
  this->newton_max_iterations = 30;
  this->newton_residual_norm_threshold = 1e-8;
}

template <int dim>
void Parameters<dim>::delete_old_outputs(MPI_Comm& mpi_communicator) const
{
  // The main process will optionally delete outputs.
  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
  {
#ifdef _MSC_VER
    std::stringstream ss;
    ss << "del";
    ss << " " << this->output_file_prefix << "*.vtk";
    ss << " " << this->output_file_prefix << "*.current_solution";
    ss << " " << this->output_file_prefix << "*.prev_solution";
    ss << " " << this->output_file_prefix << "*.matrix";
    ss << " " << this->output_file_prefix << "*.rhs";
    ss << " " << this->output_file_prefix << "*.solution";
    system(ss.str().c_str());
#else
    std::stringstream ss;
    ss << "rm";
    ss << " " << this->output_file_prefix << "*.visit";
    ss << " " << this->output_file_prefix << "*.vtu";
    ss << " " << this->output_file_prefix << "*.pvtu";
    ss << " " << this->output_file_prefix << "*.vtk";
    ss << " " << this->output_file_prefix << "*.current_solution";
    ss << " " << this->output_file_prefix << "*.prev_solution";
    ss << " " << this->output_file_prefix << "*.solution";
    ss << " " << this->output_file_prefix << "*.matrix";
    ss << " " << this->output_file_prefix << "*.rhs";
    system(ss.str().c_str());
#endif
  }
}

template class Parameters<3>;
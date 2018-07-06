#include "parameters.h"

template <int dim>
Parameters<dim>::Parameters() {
  this->start_limiting_at = -1.;
  this->gas_gamma = 5. / 3.;
  this->limit = true;
  this->slope_limiter = vertexBased;
  this->output_file_prefix = "";
  this->lax_friedrich_stabilization_value = .5;
  this->current_time_step_length = 1.e-6;

  this->debug = 0;
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

  this->volume_factor = 4;
  this->time_interval_max_cells_multiplicator = 2.;
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

template <int dim>
bool Parameters<dim>::is_periodic_boundary(int boundary_id) const
{
  for (int pb = 0; pb < this->periodic_boundaries.size(); pb++)
    if (this->periodic_boundaries[pb][0] == boundary_id || this->periodic_boundaries[pb][1] == boundary_id)
      return true;
  return false;
}

template class Parameters<3>;
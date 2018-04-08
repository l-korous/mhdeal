#include "parameters.h"

template <int dim>
Parameters<dim>::Parameters() {
  this->gas_gamma = 5. / 3.;
  this->use_iterative_improvement = false;
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
}

template <int dim>
void Parameters<dim>::delete_old_outputs() const
{
  MPI_Comm mpi_communicator(MPI_COMM_WORLD);
  // The main process will optionally delete outputs.
  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
  {
#ifdef _MSC_VER
    std::stringstream ss;
    ss << "del";
    ss << " " << this->output_file_prefix << ".vtk";
    ss << " " << this->output_file_prefix << ".newton_update";
    ss << " " << this->output_file_prefix << ".current_solution";
    ss << " " << this->output_file_prefix << ".matrix";
    ss << " " << this->output_file_prefix << ".rhs";
    system(ss.str().c_str());
#else
    std::stringstream ss;
    ss << "del";
    ss << " " << this->output_file_prefix << ".visit";
    ss << " " << this->output_file_prefix << ".vtu";
    ss << " " << this->output_file_prefix << ".pvtu";
    ss << " " << this->output_file_prefix << ".vtk";
    ss << " " << this->output_file_prefix << ".newton_update";
    ss << " " << this->output_file_prefix << ".current_solution";
    ss << " " << this->output_file_prefix << ".matrix";
    ss << " " << this->output_file_prefix << ".rhs";
    system(ss.str().c_str());
#endif
  }
}

template class Parameters<3>;
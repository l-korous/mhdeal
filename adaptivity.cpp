#include "adaptivity.h"

template <int dim>
Adaptivity<dim>::Adaptivity(Parameters<dim>& parameters, MPI_Comm& mpi_communicator) :
  parameters(parameters), mpi_communicator(mpi_communicator)
{ }

template class Adaptivity<3>;

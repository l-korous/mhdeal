#ifndef _ADAPTIVITY_H
#define _ADAPTIVITY_H

#include "util.h"
#include "parameters.h"

template <int dim>
class Adaptivity
{
public:
  Adaptivity(Parameters<dim>& parameters, MPI_Comm& mpi_communicator);
  virtual bool refine_mesh(int time_step, double time, TrilinosWrappers::MPI::Vector& solution, const DoFHandler<dim>& dof_handler,
#ifdef HAVE_MPI
    parallel::distributed::Triangulation<dim>& triangulation
#else
    Triangulation<dim>& triangulation
#endif
    , const Mapping<dim>& mapping) = 0;

  protected:
  Parameters<dim>& parameters;
  MPI_Comm& mpi_communicator;
};
#endif
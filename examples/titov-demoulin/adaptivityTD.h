#ifndef _ADAPTIVITY_TD_H
#define _ADAPTIVITY_TD_H

#include "adaptivity.h"
#include "parametersTD.h"

// Initial condition
template <int dim>
class AdaptivityTD : public Adaptivity<dim>
{
public:
  AdaptivityTD(Parameters<dim>&, MPI_Comm& mpi_communicator);
  bool refine_mesh(int time_step, double time, TrilinosWrappers::MPI::Vector& solution, const DoFHandler<dim>& dof_handler,
#ifdef HAVE_MPI
    parallel::distributed::Triangulation<dim>& triangulation
#else
    Triangulation<dim>& triangulation
#endif
    , const Mapping<dim>& mapping);
  
  void calculate_jumps(TrilinosWrappers::MPI::Vector& solution, const DoFHandler<dim>& dof_handler, const Mapping<dim>& mapping, Vector<double>& gradient_indicator);
  int last_time_step;
  int adaptivity_step;

  int max_cells;
  int refine_every_nth_time_step;
  int perform_n_initial_refinements;
  double refine_threshold;
  double coarsen_threshold;
};

#endif
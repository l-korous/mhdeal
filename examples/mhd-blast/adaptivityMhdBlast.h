#ifndef _ADAPTIVITY_MHD_BLAST_H
#define _ADAPTIVITY_MHD_BLAST_H

#include "adaptivity.h"

// Initial condition
template <int dim>
class AdaptivityMhdBlast : public Adaptivity<dim>
{
public:
  AdaptivityMhdBlast(Parameters<dim>&, MPI_Comm& mpi_communicator, int max_cells, int refine_every_nth_time_step, int perform_n_initial_refinements, double refine_threshold, double coarsen_threshold);
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
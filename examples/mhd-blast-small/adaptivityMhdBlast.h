#ifndef _ADAPTIVITY_MHD_BLAST_H
#define _ADAPTIVITY_MHD_BLAST_H

#include "adaptivity.h"

// Initial condition
template <int dim>
class AdaptivityMhdBlast : public Adaptivity<dim>
{
public:
  AdaptivityMhdBlast(Parameters<dim>&,
#ifdef HAVE_MPI
    parallel::distributed::Triangulation<dim>& triangulation
#else
    Triangulation<dim>& triangulation
#endif
  );
  bool refine_mesh(int time_step, double time, TrilinosWrappers::MPI::Vector& solution, const DoFHandler<dim>& dof_handler, const Mapping<dim>& mapping);
  bool process_element(const typename Triangulation<dim>::active_cell_iterator& cell, int ith_cell, int time_step) const;
  void refine_prev_mesh(const DoFHandler<dim>& prev_dof_handler,
#ifdef HAVE_MPI
    parallel::distributed::Triangulation<dim>& triangulation
#else
    Triangulation<dim>& triangulation
#endif
  ) const;
private:
  Vector<double> prev_gradient_indicator[2];
  bool prev_adapted[2];
  int prev_max_cells[2];
  void calculate_jumps(TrilinosWrappers::MPI::Vector& solution, const DoFHandler<dim>& dof_handler, const Mapping<dim>& mapping, Vector<double>& gradient_indicator);
  int last_time_step;
  int adaptivity_step;
};

#endif
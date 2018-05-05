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
  bool refine_mesh(int time_step, TrilinosWrappers::MPI::Vector& solution, DoFHandler<dim>& dof_handler);
  bool process_element(const typename Triangulation<dim>::active_cell_iterator& cell, int ith_cell, int time_step) const;
};

#endif
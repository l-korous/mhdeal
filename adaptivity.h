#ifndef _ADAPTIVITY_H
#define _ADAPTIVITY_H

#include "util.h"
#include "parameters.h"

template <int dim>
class Adaptivity
{
public:
  Adaptivity(Parameters<dim>& parameters,
#ifdef HAVE_MPI
    parallel::distributed::Triangulation<dim>& triangulation
#else
    Triangulation<dim>& triangulation
#endif
  );
  virtual bool refine_mesh(int time_step, double time, TrilinosWrappers::MPI::Vector& solution, const DoFHandler<dim>& dof_handler, const Mapping<dim>& mapping);
  virtual bool process_element(const typename Triangulation<dim>::active_cell_iterator& cell, int ith_cell, int time_step) const = 0;

  protected:
#ifdef HAVE_MPI
  parallel::distributed::Triangulation<dim>& triangulation;
#else
  Triangulation<dim>& triangulation;
#endif
  Parameters<dim>& parameters;
};
#endif
#include "adaptivity.h"

template <int dim>
Adaptivity<dim>::Adaptivity(Parameters<dim>& parameters,
#ifdef HAVE_MPI
  parallel::distributed::Triangulation<dim>& triangulation
#else
  Triangulation<dim>& triangulation
#endif
  ) :
  parameters(parameters),
  triangulation(triangulation)
{
}

template <int dim>
bool Adaptivity<dim>::refine_mesh(int time_step, double time, TrilinosWrappers::MPI::Vector& solution, const DoFHandler<dim>& dof_handler, const Mapping<dim>& mapping)
{
  bool toReturn = false;
  int ith_cell = 0;
  for (typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(); cell != triangulation.end(); ++cell)
  {
    if (!cell->is_locally_owned())
      continue;
    if(this->process_element(cell, ith_cell++, time_step))
      toReturn = true;
  }
  return toReturn;
}

template class Adaptivity<3>;

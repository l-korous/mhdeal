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
bool Adaptivity<dim>::refine_mesh(int time_step, TrilinosWrappers::MPI::Vector& solution, DoFHandler<dim>& dof_handler)
{
  bool toReturn = false;
  int ith_cell = 0;
  for (Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(); cell != triangulation.end(); ++cell)
  {
    if (!cell->is_locally_owned())
      continue;
    if(this->process_element(cell, ith_cell++, time_step))
      toReturn = true;
  }
  return toReturn;
}

template class Adaptivity<3>;

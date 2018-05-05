#include "adaptivityMhdBlast.h"

template <int dim>
AdaptivityMhdBlast<dim>::AdaptivityMhdBlast(Parameters<dim>& parameters,
#ifdef HAVE_MPI
  parallel::distributed::Triangulation<dim>& triangulation
#else
  Triangulation<dim>& triangulation
#endif
) :
  Adaptivity<dim>(parameters, triangulation)
{
}

template <int dim>
bool AdaptivityMhdBlast<dim>::refine_mesh(int time_step, TrilinosWrappers::MPI::Vector& solution, DoFHandler<dim>& dof_handler)
{
  if (dof_handler.n_dofs() > this->parameters.dof_threshold)
    return false;

  Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
  KellyErrorEstimator<dim>::estimate(dof_handler, QGauss<dim - 1>(this->parameters.quadrature_order), typename FunctionMap<dim>::type(), solution, estimated_error_per_cell);
  if(estimated_error_per_cell.l1_norm() < SMALL)
    return false;
  if (estimated_error_per_cell.all_zero())
    return false;
  
#ifdef HAVE_MPI
  parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(triangulation, estimated_error_per_cell, 0.3, 0.03);
#else
  GridRefinement::refine_and_coarsen_fixed_fraction(triangulation, estimated_error_per_cell, 0.3, 0.1);
#endif

  return true;
}

template <int dim>
bool AdaptivityMhdBlast<dim>::process_element(const typename Triangulation<dim>::active_cell_iterator& cell, int ith_cell, int time_step) const
{
  bool toReturn = false;
  bool refine = true;// ((time_step % 10 == 0) && time_step < 70);
  if (refine)
  {
    for (unsigned int vertex_i = 0; vertex_i < GeometryInfo<dim>::vertices_per_cell; ++vertex_i)
    {
      if (std::abs(cell->vertex(vertex_i).norm() - 0.1) < 0.03)
      {
        cell->set_refine_flag(RefinementPossibilities<dim>::cut_xy);
        toReturn = true;
      }
    }
  }
  return toReturn;
}

template class AdaptivityMhdBlast<3>;

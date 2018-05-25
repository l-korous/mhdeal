#include "adaptivityMhdBlast.h"
#include "dealiiExtensions.h"

template <int dim>
AdaptivityMhdBlast<dim>::AdaptivityMhdBlast(Parameters<dim>& parameters, int max_cells, int refine_every_nth_time_step, int perform_n_initial_refinements, double refine_threshold, double coarsen_threshold
) :
  Adaptivity<dim>(parameters),
  last_time_step(0),
  adaptivity_step(0),
  max_cells(max_cells), refine_every_nth_time_step(refine_every_nth_time_step), perform_n_initial_refinements(perform_n_initial_refinements), refine_threshold(refine_threshold), coarsen_threshold(coarsen_threshold)
{
  prev_adapted[1] = prev_adapted[0] = false;
}

template <int dim>
void AdaptivityMhdBlast<dim>::calculate_jumps(TrilinosWrappers::MPI::Vector& solution, const DoFHandler<dim>& dof_handler, const Mapping<dim>& mapping, Vector<double>& gradient_indicator)
{
  FEValuesExtractors::Scalar scalars[2];
  scalars[0].component = 0;
  scalars[1].component = 4;
  const QGauss<dim - 1> face_quadrature(1);
  UpdateFlags face_update_flags = UpdateFlags(update_values | update_JxW_values);
  FEFaceValues<dim> fe_v_face(mapping, dof_handler.get_fe(), face_quadrature, face_update_flags);
  FESubfaceValues<dim> fe_v_subface(mapping, dof_handler.get_fe(), face_quadrature, face_update_flags);
  FEFaceValues<dim> fe_v_face_neighbor(mapping, dof_handler.get_fe(), face_quadrature, update_values);
  for (typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
  {
    if (!cell->is_locally_owned())
      continue;
    Point<dim> jump;
    Point<dim> area;
    for (unsigned int face_no = 0; face_no < 4; ++face_no)
    {
      typename DoFHandler<dim>::face_iterator face = cell->face(face_no);
      if (!face->at_boundary())
      {
        Assert(cell->neighbor(face_no).state() == IteratorState::valid, ExcInternalError());
        typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face_no);
        std::vector<double> u(fe_v_face.n_quadrature_points);
        std::vector<double> u_neighbor(fe_v_face.n_quadrature_points);
        if (face->has_children())
        {
          unsigned int neighbor2 = cell->neighbor_face_no(face_no);
          for (unsigned int subface_no = 0; subface_no < face->number_of_children(); ++subface_no)
          {
            typename DoFHandler<dim>::cell_iterator neighbor_child = cell->neighbor_child_on_subface(face_no, subface_no);
            Assert(!neighbor_child->has_children(), ExcInternalError());
            fe_v_subface.reinit(cell, face_no, subface_no);
            fe_v_face_neighbor.reinit(neighbor_child, neighbor2);
            const std::vector<double> &JxW = fe_v_subface.get_JxW_values();
            for (unsigned int x = 0; x < fe_v_subface.n_quadrature_points; ++x)
              area[face_no / 2] += JxW[x];
            for (int scalar_i = 0; scalar_i < 2; scalar_i++)
            {
              fe_v_subface[scalars[scalar_i]].get_function_values(solution, u);
              fe_v_face_neighbor[scalars[scalar_i]].get_function_values(solution, u_neighbor);

              for (unsigned int x = 0; x < fe_v_subface.n_quadrature_points; ++x)
                jump[face_no / 2] += std::fabs(u[x] - u_neighbor[x])*JxW[x];
            }
          }
        }
        else
        {
          if (!cell->neighbor_is_coarser(face_no))
          {
            unsigned int neighbor2 = cell->neighbor_of_neighbor(face_no);
            fe_v_face.reinit(cell, face_no);
            fe_v_face_neighbor.reinit(neighbor, neighbor2);
            const std::vector<double> &JxW = fe_v_face.get_JxW_values();
            for (unsigned int x = 0; x < fe_v_face.n_quadrature_points; ++x)
              area[face_no / 2] += JxW[x];
            for (int scalar_i = 0; scalar_i < 2; scalar_i++)
            {
              fe_v_face[scalars[scalar_i]].get_function_values(solution, u);
              fe_v_face_neighbor[scalars[scalar_i]].get_function_values(solution, u_neighbor);

              for (unsigned int x = 0; x < fe_v_face.n_quadrature_points; ++x)
                jump[face_no / 2] += std::fabs(u[x] - u_neighbor[x])*JxW[x];
            }
          }
          else //i.e. neighbor is coarser than cell
          {
            std::pair<unsigned int, unsigned int> neighbor_face_subface
              = cell->neighbor_of_coarser_neighbor(face_no);
            Assert(neighbor_face_subface.first < GeometryInfo<dim>::faces_per_cell, ExcInternalError());
            Assert(neighbor_face_subface.second < neighbor->face(neighbor_face_subface.first)->number_of_children(),
              ExcInternalError());
            Assert(neighbor->neighbor_child_on_subface(neighbor_face_subface.first, neighbor_face_subface.second)
              == cell, ExcInternalError());
            fe_v_face.reinit(cell, face_no);
            fe_v_subface.reinit(neighbor, neighbor_face_subface.first,
              neighbor_face_subface.second);
            const std::vector<double> &JxW = fe_v_face.get_JxW_values();
            for (unsigned int x = 0; x < fe_v_face.n_quadrature_points; ++x)
              area[face_no / 2] += JxW[x];
            for (int scalar_i = 0; scalar_i < 2; scalar_i++)
            {
              fe_v_face[scalars[scalar_i]].get_function_values(solution, u);
              fe_v_subface[scalars[scalar_i]].get_function_values(solution, u_neighbor);

              for (unsigned int x = 0; x < fe_v_face.n_quadrature_points; ++x)
                jump[face_no / 2] += std::fabs(u[x] - u_neighbor[x])*JxW[x];
            }
          }
        }
      }
    }
    double average_jumps[dim];
    double sum_of_average_jumps = 0.;
    for (unsigned int i = 0; i < 2; ++i)
    {
      average_jumps[i] = jump(i) / area(i);
      sum_of_average_jumps += average_jumps[i];
    }
    gradient_indicator(cell->active_cell_index()) = sum_of_average_jumps * cell->diameter() * cell->diameter() * cell->diameter();
  }
  for (int i = 0; i < gradient_indicator.size(); i++)
    if (gradient_indicator[i] < SMALL)
      gradient_indicator[i] = 0.;
}

template <int dim>
bool AdaptivityMhdBlast<dim>::refine_prev_mesh(const DoFHandler<dim>& prev_dof_handler,
#ifdef HAVE_MPI
  parallel::distributed::Triangulation<dim>& prev_triangulation
#else
  Triangulation<dim>& prev_triangulation
#endif
) const
{
  if (this->parameters.debug & this->parameters.Adaptivity)
    LOGL(2, "prev adapted: " << (prev_adapted[1] ? "YES" : "NO"));
  if (!prev_adapted[1])
    return false;

  GridRefinement::refine_and_coarsen_fixed_fraction(prev_triangulation, prev_gradient_indicator[1], this->refine_threshold, this->coarsen_threshold, prev_max_cells[1]);

  // If possible, use aniso (only for non-distributed triangulation).
#ifndef HAVE_MPI
  for (typename DoFHandler<dim>::active_cell_iterator cell = prev_dof_handler.begin_active(); cell != prev_dof_handler.end(); ++cell)
    if (cell->is_locally_owned())
      if (cell->refine_flag_set())
        cell->set_refine_flag(RefinementPossibilities<dim>::cut_xy);
#endif

  // Fix for small errors.
  for (typename DoFHandler<dim>::active_cell_iterator cell = prev_dof_handler.begin_active(); cell != prev_dof_handler.end(); ++cell)
    if (cell->is_locally_owned())
      if (cell->refine_flag_set() && (prev_gradient_indicator[1](cell->active_cell_index()) < NEGLIGIBLE))
        cell->clear_refine_flag();

  prev_triangulation.prepare_coarsening_and_refinement();

  // Fix for periodic boundaries.
  if (this->parameters.periodic_boundaries.size() > 0)
  {
    DealIIExtensions::PeriodicCellMap<dim> periodic_cell_map;
    for (std::vector<std::array<int, 3> >::const_iterator it = this->parameters.periodic_boundaries.begin(); it != this->parameters.periodic_boundaries.end(); it++)
      DealIIExtensions::make_periodicity_map_dg(prev_dof_handler, (*it)[0], (*it)[1], (*it)[2], periodic_cell_map);
    for (typename DoFHandler<dim>::active_cell_iterator cell = prev_dof_handler.begin_active(); cell != prev_dof_handler.end(); ++cell)
    {
      if (!cell->is_locally_owned())
        continue;

      for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
      {
        if (this->parameters.is_periodic_boundary(cell->face(face_no)->boundary_id()))
        {
          if (cell->refine_flag_set())
          {
            if ((this->parameters.debug & this->parameters.Adaptivity) && (this->parameters.debug & this->parameters.PeriodicBoundaries))
              LOGL(2, "prev_cell refined: " << cell->active_cell_index());
            const DealIIExtensions::FacePair<dim>&  face_pair = periodic_cell_map.find(std::make_pair(cell, face_no))->second;
            typename DoFHandler<dim>::active_cell_iterator neighbor(cell);
            auto this_cell_index = cell->active_cell_index();
            auto zeroth_found_cell_index = (*(face_pair.cell[0])).active_cell_index();
            neighbor = ((zeroth_found_cell_index == this_cell_index && face_no == face_pair.face_idx[0]) ? face_pair.cell[1] : face_pair.cell[0]);
            if (cell->refine_flag_set() || neighbor->refine_flag_set())
            {
              if ((this->parameters.debug & this->parameters.Adaptivity) && (this->parameters.debug & this->parameters.PeriodicBoundaries))
                LOGL(2, "prev_neighbor refined: " << neighbor->active_cell_index());
              cell->set_refine_flag(RefinementPossibilities<dim>::cut_xyz);
              neighbor->clear_coarsen_flag();
              neighbor->set_refine_flag(RefinementPossibilities<dim>::cut_xyz);
            }
            else
              neighbor->clear_coarsen_flag();
          }
          cell->clear_coarsen_flag();
        }
        if (cell->level() > 0)
        {
          if (this->parameters.is_periodic_boundary(cell->parent()->face(face_no)->boundary_id()))
          {
            for (int i = 0; i < cell->parent()->n_children(); i++)
              cell->parent()->child(i)->clear_coarsen_flag();
          }
        }
      }
    }
    prev_triangulation.prepare_coarsening_and_refinement();
  }

  return true;
}

template <int dim>
bool AdaptivityMhdBlast<dim>::refine_mesh(int time_step, double time, TrilinosWrappers::MPI::Vector& solution, const DoFHandler<dim>& dof_handler,
#ifdef HAVE_MPI
  parallel::distributed::Triangulation<dim>& triangulation
#else
  Triangulation<dim>& triangulation
#endif
  , const Mapping<dim>& mapping)
{
  if (time_step % this->refine_every_nth_time_step)
  {
    prev_adapted[1] = prev_adapted[0];
    prev_gradient_indicator[1] = prev_gradient_indicator[0];
    prev_max_cells[1] = prev_max_cells[0];
    prev_adapted[0] = false;
    return false;
  }

  if (++adaptivity_step > (time_step == 0 ? this->perform_n_initial_refinements : 1))
  {
    adaptivity_step = 0;
    prev_adapted[1] = prev_adapted[0];
    prev_gradient_indicator[1] = prev_gradient_indicator[0];
    prev_max_cells[1] = prev_max_cells[0];
    prev_adapted[0] = false;
    return false;
  }

  Vector<double> gradient_indicator(triangulation.n_active_cells());
  calculate_jumps(solution, dof_handler, mapping, gradient_indicator);
  prev_gradient_indicator[1] = prev_gradient_indicator[0];
  prev_gradient_indicator[0] = gradient_indicator;
  prev_adapted[1] = prev_adapted[0];
  prev_adapted[0] = true;
  prev_max_cells[1] = prev_max_cells[0];
  prev_max_cells[0] = max_cells + (int)std::floor((time / 0.5) * 10000.);

  GridRefinement::refine_and_coarsen_fixed_fraction(triangulation, gradient_indicator, this->refine_threshold, this->coarsen_threshold, prev_max_cells[0]);
  
  triangulation.prepare_coarsening_and_refinement();

  // If possible, use aniso (only for non-distributed triangulation).
#ifndef HAVE_MPI
  for (typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
    if (cell->is_locally_owned())
      if (cell->refine_flag_set())
        cell->set_refine_flag(RefinementPossibilities<dim>::cut_xy);
#endif

  // Fix for small errors.
  for (typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
    if (cell->is_locally_owned())
      if (cell->refine_flag_set() && (gradient_indicator(cell->active_cell_index()) < NEGLIGIBLE))
        cell->clear_refine_flag();

  // Fix for periodic boundaries.
  if (this->parameters.periodic_boundaries.size() > 0)
  {
    DealIIExtensions::PeriodicCellMap<dim> periodic_cell_map;
    for (std::vector<std::array<int, 3> >::const_iterator it = this->parameters.periodic_boundaries.begin(); it != this->parameters.periodic_boundaries.end(); it++)
      DealIIExtensions::make_periodicity_map_dg(dof_handler, (*it)[0], (*it)[1], (*it)[2], periodic_cell_map);
    for (typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
    {
      if (!cell->is_locally_owned())
        continue;

      for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
      {
        if (this->parameters.is_periodic_boundary(cell->face(face_no)->boundary_id()))
        {
          if (cell->refine_flag_set())
          {
            if ((this->parameters.debug & this->parameters.Adaptivity) && (this->parameters.debug & this->parameters.PeriodicBoundaries))
              LOGL(2, "cell refined: " << cell->active_cell_index());
            const DealIIExtensions::FacePair<dim>&  face_pair = periodic_cell_map.find(std::make_pair(cell, face_no))->second;
            typename DoFHandler<dim>::active_cell_iterator neighbor(cell);
            auto this_cell_index = cell->active_cell_index();
            auto zeroth_found_cell_index = (*(face_pair.cell[0])).active_cell_index();
            neighbor = ((zeroth_found_cell_index == this_cell_index && face_no == face_pair.face_idx[0]) ? face_pair.cell[1] : face_pair.cell[0]);
            if (cell->refine_flag_set() || neighbor->refine_flag_set())
            {
              if ((this->parameters.debug & this->parameters.Adaptivity) && (this->parameters.debug & this->parameters.PeriodicBoundaries))
                LOGL(2, "neighbor refined: " << neighbor->active_cell_index());
              cell->set_refine_flag(RefinementPossibilities<dim>::cut_xyz);
              neighbor->clear_coarsen_flag();
              neighbor->set_refine_flag(RefinementPossibilities<dim>::cut_xyz);
            }
            else
              neighbor->clear_coarsen_flag();
          }
          cell->clear_coarsen_flag();
        }
        if (cell->level() > 0)
        {
          if (this->parameters.is_periodic_boundary(cell->parent()->face(face_no)->boundary_id()))
          {
            for (int i = 0; i < cell->parent()->n_children(); i++)
              cell->parent()->child(i)->clear_coarsen_flag();
          }
        }
      }
    }
    triangulation.prepare_coarsening_and_refinement();
  }

  return true;
}

template <int dim>
bool AdaptivityMhdBlast<dim>::process_element(const typename Triangulation<dim>::active_cell_iterator& cell, int ith_cell, int time_step) const
{
  bool toReturn = false;
  bool refine = true;
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

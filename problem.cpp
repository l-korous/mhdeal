#include "problem.h"

template <EquationsType equationsType, int dim>
Problem<equationsType, dim>::Problem(Parameters<dim>& parameters, Equations<equationsType, dim>& equations,
#ifdef HAVE_MPI
  parallel::distributed::Triangulation<dim>& triangulation_,
#else
  Triangulation<dim>& triangulation_,
#endif
  InitialCondition<equationsType, dim>& initial_condition, BoundaryConditions<equationsType, dim>& boundary_conditions, Adaptivity<dim>* adaptivity) :
  mpi_communicator(MPI_COMM_WORLD),
  parameters(parameters),
  equations(equations),
  prev_triangulation(triangulation_),
#ifdef HAVE_MPI
  triangulation(MPI_COMM_WORLD),
#endif
  initial_condition(initial_condition),
  boundary_conditions(boundary_conditions),
  mapping(),
  fe(parameters.use_div_free_space_for_B ? FESystem<dim>(FE_DG_Taylor<dim>(parameters.polynomial_order_dg), 5, FE_DG_DivFree<dim>(), 1) : FESystem<dim>(FE_DG_Taylor<dim>(parameters.polynomial_order_dg), 8)),
  prev_dof_handler(prev_triangulation),
  quadrature(parameters.quadrature_order),
  face_quadrature(parameters.quadrature_order),
  verbose_cout(std::cout, false),
  last_output_time(0.), last_snapshot_time(0.), time(0.),
  time_step(0),
  adaptivity_step(0),
  mag(dim + 2),
  update_flags(update_values | update_JxW_values | update_gradients),
  face_update_flags(update_values | update_JxW_values | update_normal_vectors | update_q_points),
  neighbor_face_update_flags(update_values | update_JxW_values | update_q_points | update_normal_vectors),
  fe_v_cell(mapping, fe, quadrature, update_flags),
  fe_v_prev_cell(mapping, fe, quadrature, update_flags),
  fe_v_face(mapping, fe, face_quadrature, face_update_flags),
  fe_v_prev_face(mapping, fe, face_quadrature, face_update_flags),
  fe_v_prev_face_neighbor(mapping, fe, face_quadrature, neighbor_face_update_flags),
  fe_v_subface(mapping, fe, face_quadrature, face_update_flags),
  fe_v_prev_subface(mapping, fe, face_quadrature, face_update_flags),
  fe_v_prev_subface_neighbor(mapping, fe, face_quadrature, neighbor_face_update_flags),
  adaptivity(adaptivity),
  refined_mesh(false)
{
  assembling_utils.set_problem(this);
  triangulation.copy_triangulation(prev_triangulation);
  dof_handler = new DoFHandler<dim>(triangulation);

  n_quadrature_points_cell = quadrature.get_points().size();
  fluxes_old.resize(n_quadrature_points_cell);
  W_prev.resize(n_quadrature_points_cell);
  n_quadrature_points_face = face_quadrature.get_points().size();
  Wplus_old.resize(n_quadrature_points_face);
  Wminus_old.resize(n_quadrature_points_face);
  normal_fluxes_old.resize(n_quadrature_points_face);

  if (parameters.num_flux_type == parameters.hlld)
    this->numFlux = new NumFluxHLLD<equationsType, dim>(this->parameters);
  else if (parameters.num_flux_type == parameters.lax_friedrich)
    this->numFlux = new NumFluxLaxFriedrich<equationsType, dim>(this->parameters);

  if (parameters.slope_limiter == parameters.vertexBased)
    this->slopeLimiter = new VertexBasedSlopeLimiter<equationsType, dim>(parameters, mapping, fe, *dof_handler, periodic_cell_map, dofs_per_cell, triangulation, dof_indices, component_ii, is_primitive);
  else if (parameters.slope_limiter == parameters.barthJespersen)
    this->slopeLimiter = new BarthJespersenSlopeLimiter<equationsType, dim>(parameters, mapping, fe, *dof_handler, periodic_cell_map, dofs_per_cell, triangulation, dof_indices, component_ii, is_primitive);
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::setup_system()
{
  dof_handler->clear();
  dof_handler->distribute_dofs(fe);
  prev_dof_handler.clear();
  prev_dof_handler.distribute_dofs(fe);
  DoFTools::extract_locally_relevant_dofs(*dof_handler, locally_relevant_dofs);
  DoFTools::extract_locally_relevant_dofs(prev_dof_handler, prev_locally_relevant_dofs);
  prev_locally_owned_dofs = prev_dof_handler.locally_owned_dofs();
  locally_owned_dofs = dof_handler->locally_owned_dofs();
  
  constraints.clear();
  constraints.reinit(locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(*dof_handler, constraints);
  DynamicSparsityPattern dsp(locally_relevant_dofs);
  DoFTools::make_sparsity_pattern(*dof_handler, dsp, constraints, false);
  
  prev_constraints.clear();
  prev_constraints.reinit(prev_locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(prev_dof_handler, prev_constraints);
  
  if (this->parameters.periodic_boundaries.size() > 0)
  {
    DynamicSparsityPattern prev_dsp(prev_locally_relevant_dofs);
    DoFTools::make_flux_sparsity_pattern(prev_dof_handler, prev_dsp, prev_constraints, false);

    std::vector<GridTools::PeriodicFacePair<typename DoFHandler<dim>::cell_iterator> > matched_pairs;
    for (std::vector<std::array<int, 3> >::const_iterator it = parameters.periodic_boundaries.begin(); it != parameters.periodic_boundaries.end(); it++)
      GridTools::collect_periodic_faces(prev_dof_handler, (*it)[0], (*it)[1], (*it)[2], matched_pairs);

    periodic_cell_map.clear();
    DealIIExtensions::make_periodicity_map_dg<DoFHandler<dim> >(matched_pairs, periodic_cell_map);
    DealIIExtensions::make_sparser_flux_sparsity_pattern(prev_dof_handler, prev_dsp, prev_constraints, parameters, periodic_cell_map);
  }

  prev_constraints.close();

  constraints.close();

#ifdef HAVE_MPI 
  SparsityTools::distribute_sparsity_pattern(dsp, dof_handler->n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);
#endif

  system_rhs.reinit(locally_owned_dofs, mpi_communicator);
  system_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);
  precalculate_global();
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::postprocess()
{
  current_limited_solution = current_unlimited_solution;
  constraints.distribute(current_limited_solution);

  this->slopeLimiter->postprocess(current_limited_solution, current_unlimited_solution);
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::calculate_cfl_condition()
{
  cfl_time_step = parameters.cfl_coefficient * GridTools::minimal_cell_diameter(this->triangulation) / this->max_signal_speed;
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::precalculate_global()
{
  dofs_per_cell = dof_handler->get_fe().dofs_per_cell;
  this->dof_indices.resize(dofs_per_cell);
  this->prev_dof_indices.resize(dofs_per_cell);
  this->prev_dof_indices_neighbor.resize(dofs_per_cell);

  for (unsigned int i = 0; i < dofs_per_cell; ++i)
  {
    const unsigned int component_i = this->fe.system_to_base_index(i).first.first;
    is_primitive[i] = this->fe.is_primitive(i);
    if (is_primitive[i])
    {
      component_ii[i] = this->fe.system_to_component_index(i).first;
      basis_fn_is_constant[i] = dynamic_cast<const FiniteElementIsConstantInterface<dim>*>(&(this->fe.base_element(component_i)))->is_constant(this->fe.system_to_component_index(i).second);
    }
    else
    {
      component_ii[i] = 999;
      basis_fn_is_constant[i] = false;
    }
  }
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::assemble_system(bool assemble_matrix)
{
  this->max_signal_speed = 0.;

  // Local (cell) matrices and rhs - for the currently assembled element and the neighbor
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double> cell_rhs(dofs_per_cell);

  // Loop through all cells.
  int ith_cell = 0;

  const std::list<std::pair<typename DoFHandler<dim>::cell_iterator, typename DoFHandler<dim>::cell_iterator> > cell_list = 
    GridTools::get_finest_common_cells(*dof_handler, prev_dof_handler);

  typename std::list<std::pair<typename DoFHandler<dim>::cell_iterator, typename DoFHandler<dim>::cell_iterator> >::const_iterator cell_iter;
  // The refinement case and integer representing the child's index are only used for the case when
  // - the current element is finer than the previous one
  // - the previous element has a neighbor across a face that is more refined (there are more neighbors)
  // => we of course need to take into account only some of the neighbors of the previous element
  typename std::list<std::tuple<typename DoFHandler<dim>::cell_iterator, typename DoFHandler<dim>::cell_iterator, commonCellRelativeSize, RefinementCase<dim>, int > >::const_iterator subcell_iter;

  for (cell_iter = cell_list.begin(); cell_iter != cell_list.end(); ++cell_iter)
  {
    cell = cell_iter->first;
    prev_cell = cell_iter->second;

    std::list<std::tuple<typename DoFHandler<dim>::cell_iterator, typename DoFHandler<dim>::cell_iterator, commonCellRelativeSize, RefinementCase<dim>, int > > subcells;

    if (cell->active())
    {
      // Here we can take a shortcut, since the cell is active, we can make this call.
      if (!cell->is_locally_owned())
        continue;
      else if (prev_cell->active())
        subcells.push_back(std::make_tuple(cell, prev_cell, equallyRefined, RefinementPossibilities<dim>::no_refinement, 0));
      else
      {
        for (unsigned int child = 0; child < prev_cell->n_children(); ++child)
          subcells.push_back(std::make_tuple(cell, prev_cell->child(child), prevMoreRefined, RefinementPossibilities<dim>::no_refinement, 0));
      }
    }
    else
    {
      for (unsigned int child = 0; child < cell->n_children(); ++child)
        if (cell->child(child)->is_locally_owned())
          subcells.push_back(std::make_tuple(cell->child(child), prev_cell, currentMoreRefined, cell->refinement_case(), child));
    }

    for (subcell_iter = subcells.begin(); subcell_iter != subcells.end(); ++subcell_iter)
    {
      cell = std::get<0>(*subcell_iter);
      // Now we know this one is active for sure.
      if (!cell->is_locally_owned())
        continue;

      prev_cell = std::get<1>(*subcell_iter);

      std::vector<double> JxW;
      std::vector<Tensor<1, dim> > normals;
      // Cell is active, we want its quadrature points.
      if (std::get<2>(*subcell_iter) == currentMoreRefined)
      {
        fe_v_cell.reinit(cell);
        fe_v_prev_cell.reinit(prev_cell, fe_v_cell.get_quadrature_points());
        JxW = fe_v_cell.get_JxW_values();
      }
      else
      {
        fe_v_prev_cell.reinit(prev_cell);
        fe_v_cell.reinit(cell, fe_v_prev_cell.get_quadrature_points());
        JxW = fe_v_prev_cell.get_JxW_values();
      }

      if (assemble_matrix)
        cell_matrix = 0;
      cell_rhs = 0;

      cell->get_dof_indices(dof_indices);
      prev_cell->get_dof_indices(prev_dof_indices);

      if (parameters.debug)
      {
        std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": Current cell: " << ith_cell << std::endl;
        std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": - level: " << cell->level() << std::endl;
        for (unsigned int v_i = 0; v_i < GeometryInfo<dim>::vertices_per_cell; ++v_i)
          std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": \t vertex " << v_i << ": " << cell->vertex(v_i) << std::endl;
        std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": Previous cell: " << ith_cell << std::endl;
        std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": - level: " << prev_cell->level() << std::endl;
        for (unsigned int v_i = 0; v_i < GeometryInfo<dim>::vertices_per_cell; ++v_i)
          std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": \t vertex " << v_i << ": " << prev_cell->vertex(v_i) << std::endl;
      }

      ith_cell++;

      // Assemble the volumetric integrals.
      assemble_cell_term(cell_matrix, cell_rhs, assemble_matrix, JxW);

      // Assemble the face integrals.
      for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
      {
        if (parameters.debug)
          std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": Face: " << face_no << std::endl;

        // Face that is internal wrt. prev_cell
        if (std::get<2>(*subcell_iter) == currentMoreRefined && assembling_utils.is_refinement_within_current_cell(face_no, std::get<3>(*subcell_iter), std::get<4>(*subcell_iter)))
        {
          if (parameters.debug)
            std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": \t\tFace within cell" << std::endl;
          fe_v_face.reinit(cell, face_no);
          fe_v_prev_face.reinit(prev_cell, face_no, fe_v_face.get_quadrature_points());
          JxW = fe_v_face.get_JxW_values();
          normals = fe_v_face.get_all_normal_vectors();
          prev_cell->get_dof_indices(prev_dof_indices_neighbor);
          assemble_face_term(face_no, fe_v_face, fe_v_prev_face, fe_v_prev_face, false, numbers::invalid_unsigned_int, cell_rhs, JxW, normals);
        }
        // Boundary face - here we pass the boundary id
        else if (prev_cell->at_boundary(face_no))
        {
          if (parameters.debug)
            std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": \tboundary" << std::endl;

          if (parameters.is_periodic_boundary(prev_cell->face(face_no)->boundary_id()))
          {
            if (parameters.debug)
              std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": \t\tperiodic" << std::endl;
            
            const DealIIExtensions::FacePair<dim>&  face_pair = periodic_cell_map.find(std::make_pair(prev_cell, face_no))->second;
            typename DoFHandler<dim>::active_cell_iterator neighbor(prev_cell);
            auto this_cell_index = prev_cell->active_cell_index();
            auto zeroth_found_cell_index = (*(face_pair.cell[0])).active_cell_index();
            neighbor = ((zeroth_found_cell_index == this_cell_index && face_no == face_pair.face_idx[0]) ? face_pair.cell[1] : face_pair.cell[0]);
            const unsigned int neighbor_face = ((zeroth_found_cell_index == this_cell_index && face_no == face_pair.face_idx[0]) ? face_pair.face_idx[1] : face_pair.face_idx[0]);

            neighbor->get_dof_indices(prev_dof_indices_neighbor);
            if (parameters.debug)
            {
              std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": \tNeighbor: " << std::endl;
              std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": - level: " << neighbor->level() << std::endl;
              for (unsigned int v_i = 0; v_i < GeometryInfo<dim>::vertices_per_cell; ++v_i)
                std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": \t vertex " << v_i << ": " << neighbor->vertex(v_i) << std::endl;
            }

            if (std::get<2>(*subcell_iter) == currentMoreRefined)
            {
              fe_v_face.reinit(cell, face_no);
              fe_v_prev_face.reinit(prev_cell, face_no, fe_v_face.get_quadrature_points());
              fe_v_prev_face_neighbor.reinit(neighbor, neighbor_face, fe_v_face.get_quadrature_points());
              JxW = fe_v_face.get_JxW_values();
              normals = fe_v_face.get_all_normal_vectors();
            }
            else
            {
              fe_v_prev_face.reinit(prev_cell, face_no);
              fe_v_face.reinit(cell, face_no, fe_v_prev_face.get_quadrature_points());
              fe_v_prev_face_neighbor.reinit(neighbor, neighbor_face, fe_v_prev_face.get_quadrature_points());
              JxW = fe_v_prev_face.get_JxW_values();
              normals = fe_v_prev_face.get_all_normal_vectors();
            }

            assemble_face_term(face_no, fe_v_face, fe_v_prev_face, fe_v_prev_face_neighbor, false, cell->face(face_no)->boundary_id(), cell_rhs, JxW, normals);
          }
          else
          {
            if (std::get<2>(*subcell_iter) == currentMoreRefined)
            {
              fe_v_face.reinit(cell, face_no);
              fe_v_prev_face.reinit(prev_cell, face_no, fe_v_face.get_quadrature_points());
              JxW = fe_v_face.get_JxW_values();
              normals = fe_v_face.get_all_normal_vectors();
            }
            else
            {
              fe_v_prev_face.reinit(prev_cell, face_no);
              fe_v_face.reinit(cell, face_no, fe_v_prev_face.get_quadrature_points());
              JxW = fe_v_prev_face.get_JxW_values();
              normals = fe_v_prev_face.get_all_normal_vectors();
            }
            prev_cell->get_dof_indices(prev_dof_indices_neighbor);
            assemble_face_term(face_no, fe_v_face, fe_v_prev_face, fe_v_prev_face, true, cell->face(face_no)->boundary_id(), cell_rhs, JxW, normals);
          }
        }
        // Internal face, which corresponds to a real face in prev_triangulation
        else
        {
          if (prev_cell->face(face_no)->has_children())
          {
            unsigned int neighbor_face_no = prev_cell->neighbor_face_no(face_no);

            for (unsigned int subface_no = 0; subface_no < prev_cell->face(face_no)->number_of_children(); ++subface_no)
            {
              const typename DoFHandler<dim>::active_cell_iterator neighbor_child = prev_cell->neighbor_child_on_subface(face_no, subface_no);
              neighbor_child->get_dof_indices(prev_dof_indices_neighbor);

              if (parameters.debug)
              {
                std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": \t\tFace NOT within cell, neighbor: " << std::endl;
                for (unsigned int v_i = 0; v_i < GeometryInfo<dim>::vertices_per_cell; ++v_i)
                  std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": \t\t vertex " << v_i << ": " << neighbor_child->vertex(v_i) << std::endl;
              }

              Assert(!neighbor_child->has_children(), ExcInternalError());

              if ((std::get<2>(*subcell_iter) == equallyRefined) || (std::get<2>(*subcell_iter) == prevMoreRefined))
              {
                if (parameters.debug)
                  std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": \t\t - Current equally or less refined" << std::endl;
                fe_v_prev_face_neighbor.reinit(neighbor_child, neighbor_face_no);
                fe_v_prev_subface.reinit(prev_cell, face_no, subface_no, fe_v_prev_face_neighbor.get_quadrature_points());
                fe_v_face.reinit(cell, face_no, fe_v_prev_face_neighbor.get_quadrature_points());
                JxW = fe_v_prev_face_neighbor.get_JxW_values();
                normals = fe_v_prev_face_neighbor.get_all_normal_vectors();
                // Now we need to invert the normals since they are taken from the neighbor's perspective.
                for (int n_i = 0; n_i < normals.size(); n_i++)
                  for (int n_i_j = 0; n_i_j < dim; n_i_j++)
                    normals[n_i][n_i_j] *= -1.;
                assemble_face_term(face_no, fe_v_face, fe_v_prev_subface, fe_v_prev_face_neighbor, false, numbers::invalid_unsigned_int, cell_rhs, JxW, normals);
              }
              else
              {
                if (parameters.debug)
                  std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": \t\t - Current more refined" << std::endl;
                int asdf = prev_cell->face(face_no)->number_of_children();
                if (assembling_utils.is_refinement_compatible_with_subface(face_no, std::get<3>(*subcell_iter), std::get<4>(*subcell_iter), subface_no, prev_cell->face(face_no)->refinement_case()))
                {
                  if (parameters.debug)
                    std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": \t\t - Neighbor compatible with subface" << std::endl;
                  fe_v_face.reinit(cell, face_no);
                  fe_v_prev_face.reinit(prev_cell, face_no, fe_v_face.get_quadrature_points());
                  fe_v_prev_face_neighbor.reinit(neighbor_child, neighbor_face_no, fe_v_face.get_quadrature_points());
                  JxW = fe_v_face.get_JxW_values();
                  normals = fe_v_face.get_all_normal_vectors();
                  assemble_face_term(face_no, fe_v_face, fe_v_prev_face, fe_v_prev_face_neighbor, false, numbers::invalid_unsigned_int, cell_rhs, JxW, normals);
                }
                else if (parameters.debug)
                  std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": \t\t - Neighbor NOT compatible with subface" << std::endl;

              }
            }
          }
          else if (prev_cell->neighbor_is_coarser(face_no))
          {
            if (parameters.debug)
            {
              std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": \tLess neighbors, neighbor: " << std::endl;
              std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": - level: " << prev_cell->neighbor(face_no)->level() << std::endl;
              for (unsigned int v_i = 0; v_i < GeometryInfo<dim>::vertices_per_cell; ++v_i)
                std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": \t vertex " << v_i << ": " << prev_cell->neighbor(face_no)->vertex(v_i) << std::endl;
            }

            const typename DoFHandler<dim>::cell_iterator neighbor = prev_cell->neighbor(face_no);
            neighbor->get_dof_indices(prev_dof_indices_neighbor);

            const std::pair<unsigned int, unsigned int> faceno_subfaceno = prev_cell->neighbor_of_coarser_neighbor(face_no);
            const unsigned int neighbor_face_no = faceno_subfaceno.first;

            if (std::get<2>(*subcell_iter) == equallyRefined)
            {
              fe_v_prev_face.reinit(prev_cell, face_no);
              fe_v_face.reinit(cell, face_no, fe_v_prev_face.get_quadrature_points());
              fe_v_prev_face_neighbor.reinit(neighbor, neighbor_face_no, fe_v_prev_face.get_quadrature_points());
              JxW = fe_v_prev_face.get_JxW_values();
              normals = fe_v_prev_face.get_all_normal_vectors();
            }
            else if (std::get<2>(*subcell_iter) == currentMoreRefined)
            {
              fe_v_face.reinit(cell, face_no);
              fe_v_prev_face.reinit(prev_cell, face_no, fe_v_face.get_quadrature_points());
              fe_v_prev_face_neighbor.reinit(neighbor, neighbor_face_no, fe_v_face.get_quadrature_points());
              JxW = fe_v_face.get_JxW_values();
              normals = fe_v_face.get_all_normal_vectors();
            }
            else if (std::get<2>(*subcell_iter) == prevMoreRefined)
            {
              fe_v_prev_face.reinit(prev_cell, face_no);
              fe_v_face.reinit(cell, face_no, fe_v_prev_face.get_quadrature_points());
              fe_v_prev_face_neighbor.reinit(neighbor, neighbor_face_no, fe_v_prev_face.get_quadrature_points());
              JxW = fe_v_prev_face.get_JxW_values();
              normals = fe_v_prev_face.get_all_normal_vectors();
            }

            assemble_face_term(face_no, fe_v_face, fe_v_prev_face, fe_v_prev_face_neighbor, false, numbers::invalid_unsigned_int, cell_rhs, JxW, normals);
          }
          else
          {
            const typename DoFHandler<dim>::cell_iterator neighbor = prev_cell->neighbor(face_no);
            unsigned int neighbor_face_no = prev_cell->neighbor_face_no(face_no);
            neighbor->get_dof_indices(prev_dof_indices_neighbor);

            if (parameters.debug)
            {
              std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": \tSame neighbor, neighbor: " << std::endl;
              for (unsigned int v_i = 0; v_i < GeometryInfo<dim>::vertices_per_cell; ++v_i)
                std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": \t vertex " << v_i << ": " << neighbor->vertex(v_i) << std::endl;
            }
            if (std::get<2>(*subcell_iter) == equallyRefined)
            {
              fe_v_prev_face.reinit(prev_cell, face_no);
              fe_v_face.reinit(cell, face_no, fe_v_prev_face.get_quadrature_points());
              fe_v_prev_face_neighbor.reinit(neighbor, neighbor_face_no, fe_v_prev_face.get_quadrature_points());
              JxW = fe_v_prev_face.get_JxW_values();
              normals = fe_v_prev_face.get_all_normal_vectors();
            }
            else if (std::get<2>(*subcell_iter) == currentMoreRefined)
            {
              fe_v_face.reinit(cell, face_no);
              fe_v_prev_face.reinit(prev_cell, face_no, fe_v_face.get_quadrature_points());
              fe_v_prev_face_neighbor.reinit(neighbor, neighbor_face_no, fe_v_face.get_quadrature_points());
              JxW = fe_v_face.get_JxW_values();
              normals = fe_v_face.get_all_normal_vectors();
            }
            else if (std::get<2>(*subcell_iter) == prevMoreRefined)
            {
              fe_v_prev_face.reinit(prev_cell, face_no);
              fe_v_face.reinit(cell, face_no, fe_v_prev_face.get_quadrature_points());
              fe_v_prev_face_neighbor.reinit(neighbor, neighbor_face_no, fe_v_prev_face.get_quadrature_points());
              JxW = fe_v_prev_face.get_JxW_values();
              normals = fe_v_prev_face.get_all_normal_vectors();
            }

            assemble_face_term(face_no, fe_v_face, fe_v_prev_face, fe_v_prev_face_neighbor, false, numbers::invalid_unsigned_int, cell_rhs, JxW, normals);
          }
        }
      }

      if (assemble_matrix)
        constraints.distribute_local_to_global(cell_matrix, cell_rhs, dof_indices, system_matrix, system_rhs);
      else
        constraints.distribute_local_to_global(cell_rhs, dof_indices, system_rhs);
    }
  }
  if (assemble_matrix)
    system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
}

template <EquationsType equationsType, int dim>
void
Problem<equationsType, dim>::assemble_cell_term(FullMatrix<double>& cell_matrix, Vector<double>& cell_rhs, bool assemble_matrix, std::vector<double>& JxW)
{
  if (assemble_matrix)
  {
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      for (unsigned int j = 0; j < dofs_per_cell; ++j)
      {
        double val = 0.;
        if (is_primitive[i])
        {
          if (is_primitive[j] && (component_ii[i] == component_ii[j]))
          {
            for (unsigned int q = 0; q < n_quadrature_points_cell; ++q)
              val += JxW[q] * fe_v_cell.shape_value(i, q) * fe_v_cell.shape_value(j, q);
          }
          else if (!is_primitive[j] && (component_ii[i] >= 5))
          {
            for (unsigned int q = 0; q < n_quadrature_points_cell; ++q)
            {
              Tensor<1, dim> fe_v_value = fe_v_cell[mag].value(j, q);
              val += JxW[q] * fe_v_cell.shape_value(i, q) * fe_v_value[component_ii[i] - 5];
            }
          }
        }
        else
        {
          if (is_primitive[j] && (component_ii[j] >= 5))
          {
            for (unsigned int q = 0; q < n_quadrature_points_cell; ++q)
            {
              Tensor<1, dim> fe_v_value = fe_v_cell[mag].value(i, q);
              val += JxW[q] * fe_v_cell.shape_value(j, q) * fe_v_value[component_ii[j] - 5];
            }
          }
          else if (!is_primitive[j])
          {
            for (unsigned int q = 0; q < n_quadrature_points_cell; ++q)
            {
              Tensor<1, dim> fe_v_value_i = fe_v_cell[mag].value(i, q);
              Tensor<1, dim> fe_v_value_j = fe_v_cell[mag].value(j, q);
              val += JxW[q] * fe_v_value_i * fe_v_value_j;
            }
          }
        }

        cell_matrix(i, j) += val;
      }
    }
  }

  if (time_step == 0)
    initial_condition.vector_value(fe_v_cell.get_quadrature_points(), W_prev);
  else
  {
    // Now we calculate the previous values. For this we need to employ the previous FEValues
    for (unsigned int q = 0; q < n_quadrature_points_cell; ++q)
    {
      for (unsigned int c = 0; c < Equations<equationsType, dim>::n_components; ++c)
        W_prev[q][c] = 0.;

      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        if (!is_primitive[i])
        {
          Tensor<1, dim> fe_v_value = fe_v_prev_cell[mag].value(i, q);
          for (unsigned int d = 0; d < dim; d++)
            W_prev[q][5 + d] += prev_solution(prev_dof_indices[i]) * fe_v_value[d];
        }
        else
          W_prev[q][component_ii[i]] += prev_solution(prev_dof_indices[i]) * fe_v_prev_cell.shape_value(i, q);
      }
    }
  }

  for (unsigned int q = 0; q < n_quadrature_points_cell; ++q)
    equations.compute_flux_matrix(W_prev[q], fluxes_old[q], this->parameters);

  for (unsigned int i = 0; i < dofs_per_cell; ++i)
  {
    double val = 0.;
    for (unsigned int q = 0; q < n_quadrature_points_cell; ++q)
    {
      if (is_primitive[i])
      {
        val += JxW[q] * W_prev[q][component_ii[i]] * fe_v_cell.shape_value(i, q);
        if (!basis_fn_is_constant[i])
          for (int d = 0; d < dim; d++)
            val += JxW[q] * parameters.time_step * fluxes_old[q][component_ii[i]][d] * fe_v_cell.shape_grad(i, q)[d];
      }
      else
      {
        Tensor<1, dim> fe_v_value = fe_v_cell[mag].value(i, q);
        Tensor<2, dim> fe_v_grad = fe_v_cell[mag].gradient(i, q);
        for (unsigned int d = 0; d < dim; d++)
        {
          val += JxW[q] * W_prev[q][5 + d] * fe_v_value[d];
          for (int e = 0; e < dim; e++)
            val += JxW[q] * parameters.time_step * fluxes_old[q][5 + d][e] * fe_v_grad[d][e];
        }
      }
    }
    if (std::isnan(val))
    {
      std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": isnan: " << val << std::endl;
      std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": i: " << i << ", ci: " << (!is_primitive[i] ? 1 : fe_v_cell.get_fe().system_to_component_index(i).first) << std::endl;
      std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": point: " << fe_v_cell.quadrature_point(0)[0] << ", " << fe_v_cell.quadrature_point(0)[1] << ", " << fe_v_cell.quadrature_point(0)[2] << std::endl;
      for (int j = 0; j < Equations<equationsType, dim>::n_components; j++)
        std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": W [" << j << "]: " << (double)W_prev[0][j] << std::endl;
      exit(1);
    }
    cell_rhs(i) += val;
  }
}

template <EquationsType equationsType, int dim>
void
Problem<equationsType, dim>::assemble_face_term(const unsigned int face_no, const FEFaceValuesBase<dim> &fe_v, const FEFaceValuesBase<dim> &fe_v_prev, const FEFaceValuesBase<dim> &fe_v_prev_neighbor,
  const bool external_face, const unsigned int boundary_id, Vector<double>& cell_rhs, std::vector<double>& JxW, std::vector<Tensor<1, dim> >& normals)
{
  // This loop is preparation - calculate all states (Wplus on the current element side of the currently assembled face, Wminus on the other side).
  if (time_step == 0)
  {
    initial_condition.vector_value(fe_v.get_quadrature_points(), Wplus_old);
    if (!external_face)
      for (unsigned int q = 0; q < n_quadrature_points_face; ++q)
        for (unsigned int c = 0; c < Equations<equationsType, dim>::n_components; ++c)
          Wminus_old[q][c] = Wplus_old[q][c];
  }
  else
  {
    for (unsigned int q = 0; q < n_quadrature_points_face; ++q)
    {
      for (unsigned int c = 0; c < Equations<equationsType, dim>::n_components; ++c)
        Wplus_old[q][c] = Wminus_old[q][c] = 0.;

      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        if (fe_v.get_fe().has_support_on_face(i, face_no))
        {
          if (!is_primitive[i])
          {
            Tensor<1, dim> fe_v_value = fe_v_cell[mag].value(i, q);
            for (int d = 0; d < dim; d++)
              Wplus_old[q][5 + d] += prev_solution(prev_dof_indices[i]) * fe_v_value[d];
          }
          else
            Wplus_old[q][component_ii[i]] += prev_solution(prev_dof_indices[i]) * fe_v.shape_value(i, q);
          if (!external_face)
          {
            if (!is_primitive[i])
            {
              Tensor<1, dim> fe_v_value_neighbor = fe_v_prev_neighbor[mag].value(i, q);
              for (int d = 0; d < dim; d++)
                Wminus_old[q][5 + d] += prev_solution(prev_dof_indices_neighbor[i]) * fe_v_value_neighbor[d];
            }
            else
              Wminus_old[q][component_ii[i]] += prev_solution(prev_dof_indices_neighbor[i]) * fe_v_prev_neighbor.shape_value(i, q);
          }
        }
      }
    }
  }

  for (unsigned int q = 0; q < n_quadrature_points_face; ++q)
  {
    if (external_face)
      boundary_conditions.bc_vector_value(boundary_id, fe_v.quadrature_point(q), Wminus_old[q], Wplus_old[q]);

    // Once we have the states on both sides of the face, we need to calculate the numerical flux.
    this->numFlux->numerical_normal_flux(normals[q], Wplus_old[q], Wminus_old[q], normal_fluxes_old[q], max_signal_speed);

    // Some debugging outputs.
    if (parameters.debug)
    {
      std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": \t\tpoint_i: " << q << std::endl;
      std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": \t\tq: " << fe_v.quadrature_point(q) << ", n: " << normals[q][0] << ", " << normals[q][1] << ", " << normals[q][2] << std::endl;
      std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": \t\tWplus: ";
      for (unsigned int i = 0; i < Equations<equationsType, dim>::n_components; i++)
        std::cout << Wplus_old[q][i] << (i < 7 ? ", " : "");
      std::cout << std::endl;

      std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": \t\tWminus: ";
      for (unsigned int i = 0; i < Equations<equationsType, dim>::n_components; i++)
        std::cout << Wminus_old[q][i] << (i < 7 ? ", " : "");
      std::cout << std::endl;

      std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": \t\tNum F: ";
      for (unsigned int i = 0; i < Equations<equationsType, dim>::n_components; i++)
        std::cout << normal_fluxes_old[q][i] << (i < 7 ? ", " : "");
      std::cout << std::endl;
    }
  }

  for (unsigned int i = 0; i < dofs_per_cell; ++i)
  {
    if (fe_v.get_fe().has_support_on_face(i, face_no))
    {
      double val = 0.;
      for (unsigned int q = 0; q < n_quadrature_points_face; ++q)
      {
        if (!is_primitive[i])
        {
          Tensor<1, dim> fe_v_value = fe_v_cell[mag].value(i, q);
          val += this->parameters.time_step * (normal_fluxes_old[q][5] * fe_v_value[0] + normal_fluxes_old[q][6] * fe_v_value[1] + normal_fluxes_old[q][7] * fe_v_value[2]) * JxW[q];
        }
        else
          val += this->parameters.time_step * normal_fluxes_old[q][component_ii[i]] * fe_v.shape_value(i, q) * JxW[q];

        if (std::isnan(val))
        {
          numFlux->numerical_normal_flux(normals[q], Wplus_old[q], Wminus_old[q], normal_fluxes_old[q], max_signal_speed);
          std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": isnan: " << val << std::endl;
          std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": i: " << i << ", ci: " << (!is_primitive[i] ? 1 : fe_v.get_fe().system_to_component_index(i).first) << std::endl;
          std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": point: " << fe_v.quadrature_point(q)[0] << ", " << fe_v.quadrature_point(q)[1] << ", " << fe_v.quadrature_point(q)[2] << std::endl;
          std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": normal: " << normals[q][0] << ", " << normals[q][1] << ", " << normals[q][2] << std::endl;
          for (int j = 0; j < Equations<equationsType, dim>::n_components; j++)
            std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": W+ [" << j << "]: " << (double)Wplus_old[q][j] << ", W- [" << j << "]: " << (double)Wminus_old[q][j] << ", F [" << j << "]: " << (double)normal_fluxes_old[q][j] << std::endl;
          exit(1);
        }
      }

      cell_rhs(i) -= val;
    }
  }
}

template <EquationsType equationsType, int dim>
void
Problem<equationsType, dim>::solve()
{
  // Direct solver is only usable without MPI, as it is not distributed.
#ifndef HAVE_MPI
  if (parameters.solver == parameters.direct)
  {
    SolverControl solver_control(1, 0);
    TrilinosWrappers::SolverDirect::AdditionalData data(parameters.output == Parameters<dim>::verbose_solver);
    TrilinosWrappers::SolverDirect direct(solver_control, data);
    direct.solve(system_matrix, current_unlimited_solution, system_rhs);
    return;
  }
  else
#endif
  {
    AztecOO solver;
    LinearAlgebraTrilinos::MPI::Vector completely_distributed_solution(locally_owned_dofs, mpi_communicator);

    Epetra_Vector x(View, system_matrix.trilinos_matrix().DomainMap(), completely_distributed_solution.begin());
    Epetra_Vector b(View, system_matrix.trilinos_matrix().RangeMap(), system_rhs.begin());

    solver.SetAztecOption(AZ_output, (parameters.output == Parameters<dim>::quiet_solver ? AZ_none : AZ_all));
    solver.SetAztecOption(AZ_solver, AZ_gmres);
    solver.SetRHS(&b);
    solver.SetLHS(&x);

    solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
    solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
    solver.SetAztecOption(AZ_overlap, 0);
    solver.SetAztecOption(AZ_reorder, 0);
    solver.SetAztecParam(AZ_drop, parameters.ilut_drop);
    solver.SetAztecParam(AZ_ilut_fill, parameters.ilut_fill);
    solver.SetAztecParam(AZ_athresh, parameters.ilut_atol);
    solver.SetAztecParam(AZ_rthresh, parameters.ilut_rtol);

    solver.SetUserMatrix(const_cast<Epetra_CrsMatrix *> (&system_matrix.trilinos_matrix()));

    solver.Iterate(parameters.max_iterations, parameters.linear_residual);

    constraints.distribute(completely_distributed_solution);
    current_unlimited_solution = completely_distributed_solution;
  }
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::output_base()
{
  for (int i = 0; i < prev_solution.size(); i++)
  {
    typename Equations<equationsType, dim>::Postprocessor postprocessor(parameters);
    DataOut<dim> data_out;
    data_out.attach_dof_handler(*dof_handler);

    for (int j = 0; j < prev_solution.size(); j++)
      prev_solution(j) = 0.;
    prev_solution(i) = 1.;

    // Solution components.
    data_out.add_data_vector(prev_solution, equations.component_names(), DataOut<dim>::type_dof_data, equations.component_interpretation());

    data_out.build_patches(7);

    static unsigned int output_file_number_base = 0;

    std::string filename = "solution-" + Utilities::int_to_string(output_file_number_base, 3) + ".vtk";
    std::ofstream output(filename.c_str());
    data_out.write_vtk(output);

    ++output_file_number_base;
  }
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::output_results() const
{
  typename Equations<equationsType, dim>::Postprocessor postprocessor(parameters);

  DataOut<dim> data_out;
  data_out.attach_dof_handler(*dof_handler);

  // Solution components.
  data_out.add_data_vector(current_limited_solution, equations.component_names(), DataOut<dim>::type_dof_data, equations.component_interpretation());

  // Derived quantities.
  data_out.add_data_vector(current_limited_solution, postprocessor);

#ifdef HAVE_MPI
  // Subdomains.
  Vector<float> subdomain(triangulation.n_active_cells());
  for (unsigned int i = 0; i < subdomain.size(); ++i)
    subdomain(i) = triangulation.locally_owned_subdomain();
  data_out.add_data_vector(subdomain, "subdomain");
#endif
  
  data_out.build_patches(this->parameters.patches);

  static unsigned int output_file_number = 0;

#ifdef HAVE_MPI
  const std::string filename_base = (parameters.output_file_prefix.length() > 0 ? parameters.output_file_prefix : "solution") + "-" + Utilities::int_to_string(output_file_number, 3);

  const std::string filename = (filename_base + "-" + Utilities::int_to_string(triangulation.locally_owned_subdomain(), 4));

  std::ofstream output_vtu((filename + ".vtu").c_str());
  data_out.write_vtu(output_vtu);

  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
  {
    std::vector<std::string> filenames;
    for (unsigned int i = 0; i < Utilities::MPI::n_mpi_processes(mpi_communicator); ++i)
      filenames.push_back(filename_base + "-" + Utilities::int_to_string(i, 4) + ".vtu");

    std::ofstream pvtu_master_output((filename_base + ".pvtu").c_str());
    data_out.write_pvtu_record(pvtu_master_output, filenames);

    std::ofstream visit_master_output((filename_base + ".visit").c_str());
    data_out.write_pvtu_record(visit_master_output, filenames);
}
#else
  std::string filename = (parameters.output_file_prefix.length() > 0 ? parameters.output_file_prefix : "solution") + "-" + Utilities::int_to_string(output_file_number, 3) + ".vtk";
  std::ofstream output(filename.c_str());
  data_out.write_vtk(output);
#endif

  ++output_file_number;
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::setup_initial_solution()
{
  prev_solution = 0;
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::output_matrix(TrilinosWrappers::SparseMatrix& mat, const char* suffix) const
{
  std::ofstream m;
  std::stringstream ssm;
  ssm << time_step << "-" << Utilities::MPI::this_mpi_process(mpi_communicator) << "." << suffix;
  m.open(ssm.str());
  mat.print(m);
  m.close();
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::output_vector(TrilinosWrappers::MPI::Vector& vec, const char* suffix) const
{
  std::ofstream n;
  std::stringstream ssn;
  ssn << time_step << "-" << Utilities::MPI::this_mpi_process(mpi_communicator) << "." << suffix;
  n.open(ssn.str());
  vec.print(n, 10, false, false);
  n.close();
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::run()
{
  // Preparations.
  setup_system();
  current_limited_solution.reinit(locally_relevant_dofs, mpi_communicator);
  current_unlimited_solution.reinit(locally_relevant_dofs, mpi_communicator);
  setup_initial_solution();

#ifdef OUTPUT_BASE
  output_base();
  exit(1);
#endif

  while (time < parameters.final_time)
  {
    // Some output.
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": Step: " << time_step << ", T: " << time << std::endl;
      std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ": - adaptivity step: " << adaptivity_step << std::endl;
      std::cout << Utilities::MPI::this_mpi_process(mpi_communicator) << ":    Number of active cells:       " << triangulation.n_global_active_cells() << std::endl << "   Number of degrees of freedom: " << dof_handler->n_dofs() << std::endl << std::endl;
    }

    // Assemble
    system_rhs = 0;
    system_matrix = 0;

    assemble_system();

    // Output matrix & rhs if required.
    if (parameters.output_matrix)
      output_matrix(system_matrix, "matrix");
    if (parameters.output_rhs)
      output_vector(system_rhs, "rhs");

    // Solve
    solve();

    // Output solution if required
    if (parameters.output_solution)
      output_vector(current_unlimited_solution, "current_unlimited_solution");

    // Postprocess if required
    if ((this->time > this->parameters.start_limiting_at) && parameters.limit && parameters.polynomial_order_dg > 0)
      postprocess();
    else
      current_limited_solution = current_unlimited_solution;

    // Output VTK, calculate CFL, etc.
    handle_outputs();

    this->refined_mesh = this->adaptivity->refine_mesh(this->time_step, this->time, this->current_limited_solution, *dof_handler, triangulation, this->mapping);
    if (this->refined_mesh)
    {
      // Use the old dof_handler to prepare for solution transfer.
#ifdef HAVE_MPI
      parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector> soltrans(prev_dof_handler);
#else
      SolutionTransfer<dim, TrilinosWrappers::MPI::Vector> soltrans(prev_dof_handler);
#endif
      // Replace the old triangulation by the current one before the current one gets refined.
      bool prev_mesh_refined = this->adaptivity->refine_prev_mesh(prev_dof_handler, prev_triangulation);

      // If this is the initial step (where we are looking for a good mesh to capture the initial condition, 
      // then we do not need to transfer anything.
      if (prev_mesh_refined)
      {
        if (time_step > 0)
          soltrans.prepare_for_coarsening_and_refinement(prev_solution);

        // Only now update the prev mesh (after potential call to SolutionTransfer).
        prev_triangulation.execute_coarsening_and_refinement();
      }

      // Refine the current triangulation.
      triangulation.execute_coarsening_and_refinement();

      // Process the refinement (DOFs, matrix, rhs, vectors)
      this->setup_system();

      current_limited_solution.reinit(locally_relevant_dofs, mpi_communicator);
      current_unlimited_solution.reinit(locally_relevant_dofs, mpi_communicator);

      // The rest is not meaningful for the initial step.
      // What we do is that we do not progress the solution, we instead keep the previous one, and interpolate it on the adapted mesh.
      // If that mesh is not adapted, we do not do it of course.
      if ((time_step > 0) && prev_mesh_refined)
      {
        // Now interpolate the solution
        TrilinosWrappers::MPI::Vector interpolated_solution;
        interpolated_solution.reinit(prev_locally_owned_dofs, mpi_communicator);
#ifdef HAVE_MPI
        soltrans.interpolate(interpolated_solution);
#else
        soltrans.interpolate(prev_solution, interpolated_solution);
#endif

        // Put the interpolated solution back to the previous one.
        prev_solution.reinit(prev_locally_relevant_dofs, mpi_communicator);
        this->prev_solution = interpolated_solution;
        prev_constraints.distribute(prev_solution);
      }
    }
    else
    {
      if (this->adaptivity->refine_prev_mesh(prev_dof_handler, prev_triangulation))
        prev_triangulation.execute_coarsening_and_refinement();
      this->setup_system();
      
      LinearAlgebraTrilinos::MPI::Vector completely_distributed_solution(prev_locally_owned_dofs, mpi_communicator);
      completely_distributed_solution = this->current_limited_solution;
      prev_constraints.distribute(completely_distributed_solution);
      prev_solution.reinit(prev_locally_relevant_dofs, mpi_communicator);
      this->prev_solution = completely_distributed_solution;

      current_limited_solution.reinit(locally_relevant_dofs, mpi_communicator);
      current_unlimited_solution.reinit(locally_relevant_dofs, mpi_communicator);

      // We do not move the solution if we are adapting. It may be a waste, but we do not move forward with a 'bad' solution.
      ++time_step;
      time += parameters.time_step;
    }
    adaptivity_step++;
  }
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::handle_outputs()
{
  if (parameters.output_solution)
    output_vector(current_limited_solution, "solution");

  if ((parameters.output_step < 0) || (time - last_output_time >= parameters.output_step))
  {
    output_results();
    last_output_time = time;
  }

  calculate_cfl_condition();
  double global_cfl_time_step = Utilities::MPI::min(this->cfl_time_step, this->mpi_communicator);
  parameters.time_step = global_cfl_time_step;
}

template class Problem<EquationsTypeMhd, 3>;

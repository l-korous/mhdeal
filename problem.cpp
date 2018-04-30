#include "problem.h"

template <EquationsType equationsType, int dim>
Problem<equationsType, dim>::Problem(Parameters<dim>& parameters, Equations<equationsType, dim>& equations,
#ifdef HAVE_MPI
  parallel::distributed::Triangulation<dim>& triangulation,
#else
  Triangulation<dim>& triangulation,
#endif
  InitialCondition<equationsType, dim>& initial_condition, BoundaryConditions<equationsType, dim>& boundary_conditions) :
  mpi_communicator(MPI_COMM_WORLD),
  parameters(parameters),
  equations(equations),
  triangulation(triangulation),
  prev_triangulation(triangulation),
  initial_condition(initial_condition),
  boundary_conditions(boundary_conditions),
  mapping(),
  fe(parameters.use_div_free_space_for_B ? FESystem<dim>(FE_DG_Taylor<dim>(parameters.polynomial_order_dg), 5, FE_DG_DivFree<dim>(), 1) : FESystem<dim>(FE_DG_Taylor<dim>(parameters.polynomial_order_dg), 8)),
  dof_handler(triangulation),
  prev_dof_handler(prev_triangulation),
  quadrature(parameters.quadrature_order),
  face_quadrature(parameters.quadrature_order),
  verbose_cout(std::cout, false),
  initial_step(true),
  last_output_time(0.), last_snapshot_time(0.), time(0.),
  time_step(0),
  mag(dim + 2),
  update_flags(update_values | update_JxW_values | update_gradients),
  face_update_flags(update_values | update_JxW_values | update_normal_vectors | update_q_points),
  neighbor_face_update_flags(update_values | update_q_points),
  fe_v_cell(mapping, fe, quadrature, update_flags),
  fe_v_prev_cell(mapping, fe, quadrature, update_flags),
  fe_v_face(mapping, fe, face_quadrature, face_update_flags),
  fe_v_face_neighbor(mapping, fe, face_quadrature, neighbor_face_update_flags),
  fe_v_subface(mapping, fe, face_quadrature, face_update_flags),
  fe_v_subface_neighbor(mapping, fe, face_quadrature, neighbor_face_update_flags)
{
  n_quadrature_points_cell = quadrature.get_points().size();
  fluxes_old.resize(n_quadrature_points_cell);
  W_prev.resize(n_quadrature_points_cell);
  n_quadrature_points_face = face_quadrature.get_points().size();
  normal_fluxes_old.resize(n_quadrature_points_face);

  if (parameters.num_flux_type == parameters.hlld)
    this->numFlux = new NumFluxHLLD<equationsType, dim>(this->parameters);
  else if (parameters.num_flux_type == parameters.lax_friedrich)
    this->numFlux = new NumFluxLaxFriedrich<equationsType, dim>(this->parameters);

  if (parameters.slope_limiter == parameters.vertexBased)
    this->slopeLimiter = new VertexBasedSlopeLimiter<equationsType, dim>(parameters, mapping, fe, dof_handler, periodic_cell_map, dofs_per_cell, triangulation, dof_indices, component_ii, is_primitive);
  else if (parameters.slope_limiter == parameters.barthJespersen)
    this->slopeLimiter = new BarthJespersenSlopeLimiter<equationsType, dim>(parameters, mapping, fe, dof_handler, periodic_cell_map, dofs_per_cell, triangulation, dof_indices, component_ii, is_primitive);
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::setup_system()
{
  // This function body should not be changed - it does the usual deal.II setup
  dof_handler.clear();
  dof_handler.distribute_dofs(fe);
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);
  locally_owned_dofs = dof_handler.locally_owned_dofs();
  constraints.clear();
  constraints.reinit(locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  DynamicSparsityPattern dsp(locally_relevant_dofs);
  DoFTools::make_flux_sparsity_pattern(dof_handler, dsp, constraints, false);
  if (this->parameters.periodic_boundaries.size() > 0)
  {
    for (std::vector<std::array<int, 3> >::const_iterator it = this->parameters.periodic_boundaries.begin(); it != parameters.periodic_boundaries.end(); it++)
      DealIIExtensions::make_periodicity_map_dg(dof_handler, (*it)[0], (*it)[1], (*it)[2], periodic_cell_map);
    DealIIExtensions::make_sparser_flux_sparsity_pattern(dof_handler, dsp, constraints, parameters.periodic_boundaries, periodic_cell_map);
  }
  constraints.close();

#ifdef HAVE_MPI 
  SparsityTools::distribute_sparsity_pattern(dsp, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);
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
  dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
  this->dof_indices.resize(dofs_per_cell);
  this->dof_indices_neighbor.resize(dofs_per_cell);

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
    GridTools::get_finest_common_cells(dof_handler, prev_dof_handler);

  typename std::list<std::pair<typename DoFHandler<dim>::cell_iterator, typename DoFHandler<dim>::cell_iterator> >::const_iterator cell_iter;
  typename std::list<std::tuple<typename DoFHandler<dim>::cell_iterator, typename DoFHandler<dim>::cell_iterator, FullMatrix<double> > >::const_iterator subcell_iter;

  // Used for assembling on different meshes, as starting point for prolongation matrices.
  FullMatrix<double> unit_matrix(fe.dofs_per_cell);
  for (unsigned int i = 0; i < unit_matrix.m(); ++i)
    unit_matrix(i, i) = 1;

  for (cell_iter = cell_list.begin(); cell_iter != cell_list.end(); ++cell_iter)
  {
    cell = cell_iter->first;
    prev_cell = cell_iter->second;

    std::list<std::tuple<typename DoFHandler<dim>::cell_iterator, typename DoFHandler<dim>::cell_iterator, FullMatrix<double> > > subcells;

    if (cell->active())
    {
      // Here we can take a shortcut, since the cell is active, we can make this call.
      if (!cell->is_locally_owned())
        continue;
      else if (prev_cell->active())
        subcells.push_back(std::make_tuple(cell, prev_cell, unit_matrix));
      else
      {
        for (unsigned int child = 0; child < GeometryInfo<dim>::max_children_per_cell; ++child)
        {
          FullMatrix<double> new_matrix(fe.dofs_per_cell, fe.dofs_per_cell);
          fe.get_prolongation_matrix(child).mmult(new_matrix, unit_matrix);
          subcells.push_back(std::make_tuple(cell, prev_cell->child(child), new_matrix));
        }
      }
    }
    else
    {
      for (unsigned int child = 0; child < GeometryInfo<dim>::max_children_per_cell; ++child)
      {
        FullMatrix<double> new_matrix(fe.dofs_per_cell, fe.dofs_per_cell);
        fe.get_prolongation_matrix(child).mmult(new_matrix, unit_matrix);
        subcells.push_back(std::make_tuple(cell->child(child), prev_cell, new_matrix));
      }
    }

    for (subcell_iter = subcells.begin(); subcell_iter != subcells.end(); ++subcell_iter)
    {
      cell = std::get<0>(*subcell_iter);
      // Now we know this one is active for sure.
      if (!cell->is_locally_owned())
        continue;

      prev_cell = std::get<1>(*subcell_iter);

      fe_v_cell.reinit(cell);
      fe_v_prev_cell.reinit(prev_cell);

      if (assemble_matrix)
        cell_matrix = 0;
      cell_rhs = 0;

      cell_iter->first->get_dof_indices(dof_indices);

      if (parameters.debug)
        std::cout << "NEW CELL: " << ith_cell << std::endl;
      ith_cell++;

      // Assemble the volumetric integrals.
      assemble_cell_term(cell_matrix, cell_rhs, assemble_matrix, std::get<2>(*subcell_iter));

      // Assemble the face integrals - ONLY if this is not the initial step (where we do the projection of the initial condition).
      if (!initial_step)
      {
        for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
        {
          // Boundary face - here we pass the boundary id
          if (cell->at_boundary(face_no))
          {
            if (is_periodic_boundary(cell->face(face_no)->boundary_id(), this->parameters))
            {
              const DealIIExtensions::FacePair<dim>&  face_pair = periodic_cell_map.find(std::make_pair(cell, face_no))->second;
              typename DoFHandler<dim>::active_cell_iterator neighbor(cell);
              auto this_cell_index = cell->active_cell_index();
              auto zeroth_found_cell_index = (*(face_pair.cell[0])).active_cell_index();
              neighbor = ((zeroth_found_cell_index == this_cell_index && face_no == face_pair.face_idx[0]) ? face_pair.cell[1] : face_pair.cell[0]);
              const unsigned int neighbor_face = ((zeroth_found_cell_index == this_cell_index && face_no == face_pair.face_idx[0]) ? face_pair.face_idx[1] : face_pair.face_idx[0]);

              neighbor->get_dof_indices(dof_indices_neighbor);

              fe_v_face.reinit(cell, face_no);
              fe_v_face_neighbor.reinit(neighbor, neighbor_face);

              assemble_face_term(face_no, fe_v_face, fe_v_face_neighbor, false, cell->face(face_no)->boundary_id(), cell_rhs, std::get<2>(*subcell_iter));
            }
            else
            {
              fe_v_face.reinit(cell, face_no);
              assemble_face_term(face_no, fe_v_face, fe_v_face, true, cell->face(face_no)->boundary_id(), cell_rhs, std::get<2>(*subcell_iter));
            }
          }
          else
          {
            // Here the neighbor face is more split than the current one (has children with respect to the current face of the current element), we need to assemble sub-face by sub-face
            // Not performed if there is no adaptivity involved.
            if (cell->neighbor(face_no)->has_children())
            {
              const unsigned int neighbor2 = cell->neighbor_of_neighbor(face_no);

              for (unsigned int subface_no = 0; subface_no < cell->face(face_no)->n_children(); ++subface_no)
              {
                const typename DoFHandler<dim>::active_cell_iterator neighbor_child = cell->neighbor_child_on_subface(face_no, subface_no);

                Assert(neighbor_child->face(neighbor2) == cell->face(face_no)->child(subface_no), ExcInternalError());
                Assert(neighbor_child->has_children() == false, ExcInternalError());

                fe_v_subface.reinit(cell, face_no, subface_no);
                fe_v_face_neighbor.reinit(neighbor_child, neighbor2);
                neighbor_child->get_dof_indices(dof_indices_neighbor);

                assemble_face_term(face_no, fe_v_subface, fe_v_face_neighbor, false, numbers::invalid_unsigned_int, cell_rhs, std::get<2>(*subcell_iter));
              }
            }
            // Here the neighbor face is less split than the current one, there is some transformation needed.
            // Not performed if there is no adaptivity involved.
            else if (cell->neighbor(face_no)->level() != cell->level())
            {
              const typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face_no);
              Assert(neighbor->level() == cell->level() - 1, ExcInternalError());
              neighbor->get_dof_indices(dof_indices_neighbor);

              const std::pair<unsigned int, unsigned int> faceno_subfaceno = cell->neighbor_of_coarser_neighbor(face_no);
              const unsigned int neighbor_face_no = faceno_subfaceno.first, neighbor_subface_no = faceno_subfaceno.second;

              Assert(neighbor->neighbor_child_on_subface(neighbor_face_no, neighbor_subface_no) == cell, ExcInternalError());

              fe_v_face.reinit(cell, face_no);
              fe_v_subface_neighbor.reinit(neighbor, neighbor_face_no, neighbor_subface_no);

              assemble_face_term(face_no, fe_v_face, fe_v_subface_neighbor, false, numbers::invalid_unsigned_int, cell_rhs, std::get<2>(*subcell_iter));
            }
            // Here the neighbor face fits exactly the current face of the current element, this is the 'easy' part.
            // This is the only face assembly case performed without adaptivity.
            else
            {
              const typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face_no);
              neighbor->get_dof_indices(dof_indices_neighbor);

              fe_v_face.reinit(cell, face_no);
              fe_v_face_neighbor.reinit(neighbor, cell->neighbor_of_neighbor(face_no));
              assemble_face_term(face_no, fe_v_face, fe_v_face_neighbor, false, numbers::invalid_unsigned_int, cell_rhs, std::get<2>(*subcell_iter));
            }
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
Problem<equationsType, dim>::assemble_cell_term(FullMatrix<double>& cell_matrix, Vector<double>& cell_rhs, bool assemble_matrix, const FullMatrix<double> prolongation_matrix)
{
  // In initial step we assume that the 'previous' mesh is the same as the current one.
  if (initial_step)
  {
    initial_condition.vector_value(fe_v_cell.get_quadrature_points(), W_prev);

    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      double rhs_val = 0.;
      for (unsigned int q = 0; q < n_quadrature_points_cell; ++q)
      {
        if (is_primitive[i])
          rhs_val += fe_v_cell.JxW(q) * W_prev[q][component_ii[i]] * fe_v_cell.shape_value(i, q);
        else
        {
          dealii::Tensor<1, dim> fe_v_value = fe_v_cell[mag].value(i, q);
          for (unsigned int d = 0; d < dim; d++)
            rhs_val += fe_v_cell.JxW(q) * W_prev[q][5 + d] * fe_v_value[d];
        }
      }
      cell_rhs(i) += rhs_val;

      if (assemble_matrix)
      {
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
        {
          double matrix_val = 0.;
          for (unsigned int q = 0; q < n_quadrature_points_cell; ++q)
          {
            if (is_primitive[i])
            {
              if (is_primitive[j] && (component_ii[i] == component_ii[j]))
                matrix_val += fe_v_cell.JxW(q) * fe_v_cell.shape_value(i, q) * fe_v_cell.shape_value(j, q);
              else if (!is_primitive[j] && (component_ii[i] >= 5))
              {
                dealii::Tensor<1, dim> fe_v_value = fe_v_cell[mag].value(j, q);
                matrix_val += fe_v_cell.JxW(q) * fe_v_cell.shape_value(i, q) * fe_v_value[component_ii[i] - 5];
              }
            }
            else
            {
              if (is_primitive[j] && (component_ii[j] >= 5))
              {
                dealii::Tensor<1, dim> fe_v_value = fe_v_cell[mag].value(i, q);
                matrix_val += fe_v_cell.JxW(q) * fe_v_cell.shape_value(j, q) * fe_v_value[component_ii[j] - 5];
              }
              else if (!is_primitive[j])
              {
                dealii::Tensor<1, dim> fe_v_value_i = fe_v_cell[mag].value(i, q);
                dealii::Tensor<1, dim> fe_v_value_j = fe_v_cell[mag].value(j, q);
                matrix_val += fe_v_cell.JxW(q) * fe_v_value_i * fe_v_value_j;
              }
            }
          }
          cell_matrix(i, j) += matrix_val;
        }
      }
    }
  }
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
          dealii::Tensor<1, dim> fe_v_value = fe_v_cell[mag].value(i, q);
          for (unsigned int d = 0; d < dim; d++)
            W_prev[q][5 + d] += prev_solution(dof_indices[i]) * fe_v_value[d];
        }
        else
          W_prev[q][component_ii[i]] += prev_solution(dof_indices[i]) * fe_v_cell.shape_value(i, q);
      }

      equations.compute_flux_matrix(W_prev[q], fluxes_old[q], this->parameters);
    }

    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      double val = 0.;
      for (unsigned int q = 0; q < n_quadrature_points_cell; ++q)
      {
        if (is_primitive[i])
        {
          val += fe_v_cell.JxW(q) * W_prev[q][component_ii[i]] * fe_v_cell.shape_value(i, q);
          if (!basis_fn_is_constant[i])
            for (int d = 0; d < dim; d++)
              val += fe_v_cell.JxW(q) * parameters.time_step * fluxes_old[q][component_ii[i]][d] * fe_v_cell.shape_grad(i, q)[d];
        }
        else
        {
          dealii::Tensor<1, dim> fe_v_value = fe_v_cell[mag].value(i, q);
          dealii::Tensor<2, dim> fe_v_grad = fe_v_cell[mag].gradient(i, q);
          for (unsigned int d = 0; d < dim; d++)
          {
            val += fe_v_cell.JxW(q) * W_prev[q][5 + d] * fe_v_value[d];
            for (int e = 0; e < dim; e++)
              val += fe_v_cell.JxW(q) * parameters.time_step * fluxes_old[q][5 + d][e] * fe_v_grad[d][e];
          }
        }
      }
      cell_rhs(i) += val;
    }
  }
}

template <EquationsType equationsType, int dim>
void
Problem<equationsType, dim>::assemble_face_term(const unsigned int face_no, const FEFaceValuesBase<dim> &fe_v, const FEFaceValuesBase<dim> &fe_v_neighbor,
  const bool external_face, const unsigned int boundary_id, Vector<double>& cell_rhs, const FullMatrix<double> prolongation_matrix)
{
  if (parameters.debug)
    std::cout << "\tEdqe: " << face_no << std::endl;

  // This loop is preparation - calculate all states (Wplus on the current element side of the currently assembled face, Wminus on the other side).
  for (unsigned int q = 0; q < n_quadrature_points_face; ++q)
  {
    for (unsigned int c = 0; c < Equations<equationsType, dim>::n_components; ++c)
      Wplus_old[c] = Wminus_old[c] = 0.;

    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      if (fe_v.get_fe().has_support_on_face(i, face_no))
      {
        if (!is_primitive[i])
        {
          dealii::Tensor<1, dim> fe_v_value = fe_v_cell[mag].value(i, q);
          for (int d = 0; d < dim; d++)
            Wplus_old[5 + d] += prev_solution(dof_indices[i]) * fe_v_value[d];
        }
        else
          Wplus_old[component_ii[i]] += prev_solution(dof_indices[i]) * fe_v.shape_value(i, q);
        if (!external_face)
        {
          if (!is_primitive[i])
          {
            dealii::Tensor<1, dim> fe_v_value_neighbor = fe_v_neighbor[mag].value(i, q);
            for (int d = 0; d < dim; d++)
              Wminus_old[5 + d] += prev_solution(dof_indices_neighbor[i]) * fe_v_value_neighbor[d];
          }
          else
            Wminus_old[component_ii[i]] += prev_solution(dof_indices_neighbor[i]) * fe_v_neighbor.shape_value(i, q);
        }
      }
    }

    if (external_face)
      boundary_conditions.bc_vector_value(boundary_id, fe_v.quadrature_point(q), Wminus_old, Wplus_old);

    // Once we have the states on both sides of the face, we need to calculate the numerical flux.
    this->numFlux->numerical_normal_flux(fe_v.normal_vector(q), Wplus_old, Wminus_old, normal_fluxes_old[q], max_signal_speed);

    // Some debugging outputs.
    if (parameters.debug)
    {
      std::cout << "\t\tpoint_i: " << q << std::endl;
      std::cout << "\t\tq: " << fe_v.quadrature_point(q) << ", n: " << fe_v.normal_vector(q)[0] << ", " << fe_v.normal_vector(q)[1] << ", " << fe_v.normal_vector(q)[2] << std::endl;
      std::cout << "\t\tWplus: ";
      for (unsigned int i = 0; i < Equations<equationsType, dim>::n_components; i++)
        std::cout << Wplus_old[i] << (i < 7 ? ", " : "");
      std::cout << std::endl;

      std::cout << "\t\tWminus: ";
      for (unsigned int i = 0; i < Equations<equationsType, dim>::n_components; i++)
        std::cout << Wminus_old[i] << (i < 7 ? ", " : "");
      std::cout << std::endl;

      std::cout << "\t\tNum F: ";
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
          dealii::Tensor<1, dim> fe_v_value = fe_v_cell[mag].value(i, q);
          val += this->parameters.time_step * (normal_fluxes_old[q][5] * fe_v_value[0] + normal_fluxes_old[q][6] * fe_v_value[1] + normal_fluxes_old[q][7] * fe_v_value[2]) * fe_v.JxW(q);
        }
        else
          val += this->parameters.time_step * normal_fluxes_old[q][component_ii[i]] * fe_v.shape_value(i, q) * fe_v.JxW(q);

        if (std::isnan(val))
        {
          numFlux->numerical_normal_flux(fe_v.normal_vector(q), Wplus_old, Wminus_old, normal_fluxes_old[q], max_signal_speed);
          std::cout << "isnan: " << val << std::endl;
          std::cout << "i: " << i << ", ci: " << (!is_primitive[i] ? 1 : fe_v.get_fe().system_to_component_index(i).first) << std::endl;
          std::cout << "point: " << fe_v.quadrature_point(q)[0] << ", " << fe_v.quadrature_point(q)[1] << ", " << fe_v.quadrature_point(q)[2] << std::endl;
          std::cout << "normal: " << fe_v.normal_vector(q)[0] << ", " << fe_v.normal_vector(q)[1] << ", " << fe_v.normal_vector(q)[2] << std::endl;
          for (int j = 0; j < Equations<equationsType, dim>::n_components; j++)
            std::cout << "W+ [" << j << "]: " << (double)Wplus_old[j] << ", W- [" << j << "]: " << (double)Wminus_old[j] << ", F [" << j << "]: " << (double)normal_fluxes_old[q][j] << std::endl;
          exit(1);
        }
      }

      cell_rhs(i) -= val;
    }
  }
}

template <EquationsType equationsType, int dim>
void
Problem<equationsType, dim>::solve(TrilinosWrappers::MPI::Vector &newton_update, bool reset_matrix)
{
  // Direct solver is only usable without MPI, as it is not distributed.
#ifndef HAVE_MPI
  if (parameters.solver == parameters.direct)
  {
    SolverControl solver_control(1, 0);
    TrilinosWrappers::SolverDirect::AdditionalData data(parameters.output == Parameters<dim>::verbose_solver);
    TrilinosWrappers::SolverDirect direct(solver_control, data);
    direct.solve(system_matrix, newton_update, system_rhs);
    return;
  }
  else
#endif
  {
    dealii::LinearAlgebraTrilinos::MPI::Vector completely_distributed_solution(locally_owned_dofs, mpi_communicator);

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

    if (reset_matrix)
      solver.SetUserMatrix(const_cast<Epetra_CrsMatrix *> (&system_matrix.trilinos_matrix()));

    solver.Iterate(parameters.max_iterations, parameters.linear_residual);

    constraints.distribute(completely_distributed_solution);
    newton_update = completely_distributed_solution;
  }
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::output_base()
{
  for (int i = 0; i < prev_solution.size(); i++)
  {
    typename Equations<equationsType, dim>::Postprocessor postprocessor(parameters);
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);

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
  data_out.attach_dof_handler(dof_handler);

  // Solution components.
  data_out.add_data_vector(prev_solution, equations.component_names(), DataOut<dim>::type_dof_data, equations.component_interpretation());

  // Derived quantities.
  data_out.add_data_vector(prev_solution, postprocessor);

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
  prev_solution.reinit(locally_relevant_dofs, mpi_communicator);
  current_limited_solution.reinit(locally_owned_dofs, mpi_communicator);
  current_unlimited_solution.reinit(locally_relevant_dofs, mpi_communicator);

#ifdef HAVE_MPI
  bool should_load_from_file = false;
  double _time, _time_step;
  std::ifstream history("history");
  if (history.is_open())
  {
    std::string line;
    double corner_a_test[3];
    double corner_b_test[3];
    int ref_test[3];
    while (getline(history, line))
    {
      std::istringstream ss(line);
      ss >> _time >> _time_step;
      if (!should_load_from_file)
      {
        bool dimensions_match = true;
        for (int i = 0; i < dim; i++)
          ss >> corner_a_test[i] >> corner_b_test[i] >> ref_test[i];
        for (int i = 0; i < dim; i++)
        {
          if (this->parameters.corner_a[i] != corner_a_test[i])
            dimensions_match = false;
          if (this->parameters.corner_b[i] != corner_b_test[i])
            dimensions_match = false;
          if (this->parameters.refinements[i] != ref_test[i])
            dimensions_match = false;
        }
        if (dimensions_match)
          should_load_from_file = true;
      }
    }
    history.close();
  }

  if (should_load_from_file)
  {
    this->time = _time;
    this->last_output_time = _time;
    this->last_snapshot_time = _time;
    this->time_step = _time_step;
    load();
  }
  else {
    prev_solution = 0;
    remove("triangulation");
    remove("triangulation.info");
    remove("history");
  }
#else
  prev_solution = 0;
#endif

  current_unlimited_solution = prev_solution;
  current_limited_solution = prev_solution;
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::output_matrix(TrilinosWrappers::SparseMatrix& mat, const char* suffix, int time_step, int newton_step) const
{
  std::ofstream m;
  std::stringstream ssm;
  if (newton_step >= 0)
    ssm << time_step << "-" << newton_step << "-" << Utilities::MPI::this_mpi_process(mpi_communicator) << "." << suffix;
  else
    ssm << time_step << "-" << Utilities::MPI::this_mpi_process(mpi_communicator) << "." << suffix;
  m.open(ssm.str());
  mat.print(m);
  m.close();
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::output_vector(TrilinosWrappers::MPI::Vector& vec, const char* suffix, int time_step, int newton_step) const
{
  std::ofstream n;
  std::stringstream ssn;
  if (newton_step >= 0)
    ssn << time_step << "-" << newton_step << "-" << Utilities::MPI::this_mpi_process(mpi_communicator) << "." << suffix;
  else
    ssn << time_step << "-" << Utilities::MPI::this_mpi_process(mpi_communicator) << "." << suffix;
  n.open(ssn.str());
  vec.print(n, 10, false, false);
  n.close();
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::save()
{
#ifdef HAVE_MPI
  std::ofstream history("history");
  history << this->time << " " << this->time_step;
  for (int i = 0; i < dim; i++)
    history << " " << this->parameters.corner_a[i] << " " << this->parameters.corner_b[i] << " " << this->parameters.refinements[i];
  history << std::endl;
  history.close();

  parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector> sol_trans(dof_handler);
  sol_trans.prepare_serialization(this->prev_solution);
  this->triangulation.save("triangulation");
#endif
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::load()
{
#ifdef HAVE_MPI
  this->triangulation.load("triangulation");
  parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector> sol_trans(dof_handler);
  sol_trans.deserialize(this->prev_solution);
#endif
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::run()
{
  // Preparations.
  setup_system();
  setup_initial_solution();
  TrilinosWrappers::MPI::Vector newton_update(locally_relevant_dofs, mpi_communicator);

  // Time loop.
  double newton_damping = this->parameters.initial_and_max_newton_damping;

#ifdef OUTPUT_BASE
  output_base();
  exit(1);
#endif

  while (time < parameters.final_time)
  {
    // Some output.
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      std::cout << "Step: " << time_step << ", T: " << time << std::endl;
      if (initial_step)
        std::cout << "   Number of active cells:       " << triangulation.n_active_cells() << std::endl << "   Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl << std::endl;
    }

    // Assemble
    system_rhs = 0;
    if (initial_step)
      system_matrix = 0;

    assemble_system(initial_step);

    // Output matrix & rhs if required.
    if (parameters.output_matrix)
      output_matrix(system_matrix, "matrix", time_step);
    if (parameters.output_rhs)
      output_vector(system_rhs, "rhs", time_step);

    // Solve
    solve(current_unlimited_solution, initial_step);

    // Output solution if required
    if (parameters.output_solution)
      output_vector(current_unlimited_solution, "current_unlimited_solution", time_step);

    // Postprocess if required
    if ((this->time > this->parameters.start_limiting_at) && !initial_step && parameters.limit && parameters.polynomial_order_dg > 0)
      postprocess();
    else
      current_limited_solution = current_unlimited_solution;

    prev_solution = current_limited_solution;

    move_time_step_handle_outputs();
}
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::move_time_step_handle_outputs()
{
  if (parameters.output_solution)
    output_vector(prev_solution, "prev_solution", time_step);

  if ((parameters.output_step < 0) || (time - last_output_time >= parameters.output_step))
  {
    output_results();
    last_output_time = time;
  }

  if ((parameters.snapshot_step < 0) || (time - last_snapshot_time >= parameters.snapshot_step))
  {
    save();
    last_snapshot_time = time;
  }

  if (!initial_step)
  {
    calculate_cfl_condition();
    double global_cfl_time_step = dealii::Utilities::MPI::min(this->cfl_time_step, this->mpi_communicator);
    parameters.time_step = global_cfl_time_step;
  }

  ++time_step;
  time += parameters.time_step;
  initial_step = false;
}

template <EquationsType equationsType, int dim>
bool Problem<equationsType, dim>::is_periodic_boundary(int boundary_id, const Parameters<dim>& parameters)
{
  for (int pb = 0; pb < parameters.periodic_boundaries.size(); pb++)
    if (parameters.periodic_boundaries[pb][0] == boundary_id || parameters.periodic_boundaries[pb][1] == boundary_id)
      return true;
  return false;
}

template class Problem<EquationsTypeMhd, 3>;

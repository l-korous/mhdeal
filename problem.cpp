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
  initial_condition(initial_condition),
  boundary_conditions(boundary_conditions),
  mapping(),
  fe(parameters.use_div_free_space_for_B ? FESystem<dim>(FE_DG_Taylor<dim>(parameters.polynomial_order_dg), 5, FE_DG_DivFree<dim>(), 1) : FESystem<dim>(FE_DG_Taylor<dim>(parameters.polynomial_order_dg), 8)), dof_handler(triangulation),
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
  fe_v_face(mapping, fe, face_quadrature, face_update_flags),
  fe_v_face_neighbor(mapping, fe, face_quadrature, neighbor_face_update_flags)
{
  fe_v_cell = new FEValues<dim>(mapping, fe, quadrature, update_flags);
  n_quadrature_points_cell = quadrature.get_points().size();
  fluxes_old.resize(n_quadrature_points_cell);
  W_prev.resize(n_quadrature_points_cell);
  n_quadrature_points_face = face_quadrature.get_points().size();
  normal_fluxes_old.resize(n_quadrature_points_face);
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
  // This is what we return.
  current_limited_solution = current_unlimited_solution;
  constraints.distribute(current_limited_solution);

  int cell_count = 0;
  // Loop through all cells.
  for (typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
  {
    if (!cell->is_locally_owned())
      continue;

    bool u_c_set[Equations<equationsType, dim>::n_components];
    for (int i = 0; i < Equations<equationsType, dim>::n_components; i++)
      u_c_set[i] = false;

    double u_c[Equations<equationsType, dim>::n_components];
    cell->get_dof_indices(dof_indices);

    PostprocessData* data = 0;
    auto it = this->postprocessData.find(cell->active_cell_index());
    if (it != this->postprocessData.end())
      data = &(it->second);
    else
    {
      data = &(((postprocessData.insert(std::make_pair(cell->active_cell_index(), PostprocessData()))).first)->second);
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        if (!is_primitive[i])
          data->lambda_indices_to_multiply_all_B_components.push_back(dof_indices[i]);
        else
        {
          if (!u_c_set[component_ii[i]])
            u_c_set[component_ii[i]] = true;
          else
            data->lambda_indices_to_multiply[component_ii[i]].push_back(dof_indices[i]);
        }
      }
      for (int i = 0; i < Equations<equationsType, dim>::n_components; i++)
        u_c_set[i] = false;

      data->center = cell->center();
      for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i)
      {
        data->vertexIndex[i] = cell->vertex_index(i);
        data->vertexPoint[i] = data->center + (1. - NEGLIGIBLE) * (cell->vertex(i) - data->center);

        unsigned short neighbor_i = 0;
        for (auto neighbor_element : GridTools::find_cells_adjacent_to_vertex(triangulation, data->vertexIndex[i]))
        {
          typename DoFHandler<dim>::active_cell_iterator neighbor(&triangulation, neighbor_element->level(), neighbor_element->index(), &dof_handler);
          if (neighbor->active_cell_index() != cell->active_cell_index())
          {
            data->neighbor_dof_indices[i][neighbor_i].resize(dofs_per_cell);
            neighbor->get_dof_indices(data->neighbor_dof_indices[i][neighbor_i++]);
          }
        }
      }
    }

    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      if (is_primitive[i])
      {
        if (!u_c_set[component_ii[i]])
        {
          u_c[component_ii[i]] = current_unlimited_solution(dof_indices[i]);
          u_c_set[component_ii[i]] = true;
        }
      }
    }

    if (parameters.debug)
      std::cout << "cell: " << ++cell_count << " - center: " << data->center << ", values: " << u_c[0] << ", " << u_c[1] << ", " << u_c[2] << ", " << u_c[3] << ", " << u_c[4] << std::endl;

    double alpha_e[Equations<equationsType, dim>::n_components];
    for (int i = 0; i < Equations<equationsType, dim>::n_components; i++)
      alpha_e[i] = 1.;
    for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i)
    {
      std::set<unsigned int> visited_faces;

      // (!!!) Find out u_i
      Vector<double> u_i(Equations<equationsType, dim>::n_components);
      VectorTools::point_value(dof_handler, current_unlimited_solution, data->vertexPoint[i], u_i);

      if (parameters.debug)
      {
        std::cout << "\tv_i: " << cell->vertex(i) << ", values: ";
        for (int i = 0; i < Equations<equationsType, dim>::n_components; i++)
          std::cout << u_i[i] << (i == Equations<equationsType, dim>::n_components - 1 ? "" : ", ");
        std::cout << std::endl;
      }

      // Init u_i_min, u_i_max
      double u_i_min[Equations<equationsType, dim>::n_components];
      double u_i_max[Equations<equationsType, dim>::n_components];
      for (int k = 0; k < Equations<equationsType, dim>::n_components; k++)
      {
        u_i_min[k] = u_c[k];
        u_i_max[k] = u_c[k];
      }

      // For all vertices -> v_i
      for (auto dof_indices_neighbor : data->neighbor_dof_indices[i])
      {
        if (dof_indices_neighbor.size() == 0)
          continue;
        bool u_i_extrema_set[Equations<equationsType, dim>::n_components];
          for (int i = 0; i < Equations<equationsType, dim>::n_components; i++)
            u_i_extrema_set[i] = false;

        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          if (is_primitive[i])
          {
            if (!u_i_extrema_set[component_ii[i]])
            {
              double val = current_unlimited_solution(dof_indices_neighbor[i]);
              if (parameters.debug)
              {
                if (val < u_i_min[component_ii[i]])
                  std::cout << "\tdecreasing u_i_min to: " << val << std::endl;
                if (val > u_i_max[component_ii[i]])
                  std::cout << "\tincreasing u_i_max to: " << val << std::endl;
              }
              u_i_min[component_ii[i]] = std::min(u_i_min[component_ii[i]], val);
              u_i_max[component_ii[i]] = std::max(u_i_max[component_ii[i]], val);
              u_i_extrema_set[component_ii[i]] = true;
            }
          }
        }
      }

      // Based on u_i_min, u_i_max, u_i, get alpha_e
      for (int k = 0; k < Equations<equationsType, dim>::n_components; k++)
        if (std::abs((u_c[k] - u_i[k]) / u_c[k]) > NEGLIGIBLE)
        {
          alpha_e[k] = std::min(alpha_e[k], ((u_i[k] - u_c[k]) > 0.) ? std::min(1.0, (u_i_max[k] - u_c[k]) / (u_i[k] - u_c[k])) : std::min(1.0, (u_i_min[k] - u_c[k]) / (u_i[k] - u_c[k])));
          if (parameters.debug)
            std::cout << "\talpha_e[" << k << "]: " << alpha_e[k] << std::endl;
        }
    }

    for (int k = 0; k < Equations<equationsType, dim>::n_components; k++)
      for (int i = 0; i < data->lambda_indices_to_multiply[k].size(); i++)
        current_limited_solution(data->lambda_indices_to_multiply[k][i]) *= alpha_e[k];

    double alpha_e_B = std::min(std::min(alpha_e[5], alpha_e[6]), alpha_e[7]);
    for (int i = 0; i < data->lambda_indices_to_multiply_all_B_components.size(); i++)
      current_limited_solution(data->lambda_indices_to_multiply_all_B_components[i]) *= alpha_e_B;
  }
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::calculate_cfl_condition()
{
  cfl_time_step = parameters.cfl_constant * GridTools::minimal_cell_diameter(this->triangulation) / this->equations.max_signal_speed;
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
  this->cfl_time_step = 1.e6;

  // Local (cell) matrices and rhs - for the currently assembled element and the neighbor
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double> cell_rhs(dofs_per_cell);

  // Loop through all cells.
  int ith_cell = 0;
  for (typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
  {
    // Only assemble what belongs to this process.
    if (!cell->is_locally_owned())
      continue;

    if (this->initial_step)
      fe_v_cell->reinit(cell);

    if (assemble_matrix)
      cell_matrix = 0;
    cell_rhs = 0;

    cell->get_dof_indices(dof_indices);

    if (parameters.debug)
      std::cout << "NEW CELL: " << ith_cell++ << std::endl;

    // Assemble the volumetric integrals.
    assemble_cell_term(cell_matrix, cell_rhs, assemble_matrix);

    // Assemble the face integrals - ONLY if this is not the initial step (where we do the projection of the initial condition).
    if (!initial_step)
    {
      for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
      {
        // Boundary face - here we pass the boundary id
        if (cell->at_boundary(face_no))
        {
          if (is_periodic_boundary(cell->face(face_no)->boundary_id()))
          {
            const DealIIExtensions::FacePair<dim>&  face_pair = periodic_cell_map.find(std::make_pair(cell, face_no))->second;
            typename DoFHandler<dim>::active_cell_iterator neighbor(cell);
            neighbor = ((*(face_pair.cell[0])).active_cell_index() == (*cell).active_cell_index()) ? face_pair.cell[1] : face_pair.cell[0];
            const unsigned int neighbor_face = ((*(face_pair.cell[0])).active_cell_index() == (*cell).active_cell_index()) ? face_pair.face_idx[1] : face_pair.face_idx[0];

            neighbor->get_dof_indices(dof_indices_neighbor);

            fe_v_face.reinit(cell, face_no);
            fe_v_face_neighbor.reinit(neighbor, neighbor_face);

            assemble_face_term(face_no, fe_v_face, fe_v_face_neighbor, false, cell->face(face_no)->boundary_id(), cell_rhs);
          }
          else
          {
            fe_v_face.reinit(cell, face_no);
            assemble_face_term(face_no, fe_v_face, fe_v_face, true, cell->face(face_no)->boundary_id(), cell_rhs);
          }
        }
        else
        {
          const typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face_no);
          neighbor->get_dof_indices(dof_indices_neighbor);

          fe_v_face.reinit(cell, face_no);
          fe_v_face_neighbor.reinit(neighbor, cell->neighbor_of_neighbor(face_no));

          assemble_face_term(face_no, fe_v_face, fe_v_face_neighbor, false, numbers::invalid_unsigned_int, cell_rhs);
        }
      }
    }

    if (assemble_matrix)
      constraints.distribute_local_to_global(cell_matrix, cell_rhs, dof_indices, system_matrix, system_rhs);
    else
      constraints.distribute_local_to_global(cell_rhs, dof_indices, system_rhs);
  }

  if (assemble_matrix)
    system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
}

template <EquationsType equationsType, int dim>
void
Problem<equationsType, dim>::assemble_cell_term(FullMatrix<double>& cell_matrix, Vector<double>& cell_rhs, bool assemble_matrix)
{
  if (initial_step)
  {
    initial_condition.vector_value(fe_v_cell->get_quadrature_points(), W_prev);

    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      double val = 0.;
      for (unsigned int q = 0; q < n_quadrature_points_cell; ++q)
      {
        if (is_primitive[i])
          val += fe_v_cell->JxW(q) * W_prev[q][component_ii[i]] * fe_v_cell->shape_value(i, q);
        else
        {
          dealii::Tensor<1, dim> fe_v_value = (*fe_v_cell)[mag].value(i, q);
          for (unsigned int d = 0; d < dim; d++)
            val += fe_v_cell->JxW(q) * W_prev[q][5 + d] * fe_v_value[d];
        }

        if (assemble_matrix)
        {
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
          {
            if (is_primitive[i])
            {
              if (is_primitive[j] && (component_ii[i] == component_ii[j]))
                cell_matrix(i, j) += fe_v_cell->JxW(q) * fe_v_cell->shape_value(i, q) * fe_v_cell->shape_value(j, q);
              else if (!is_primitive[j] && (component_ii[i] >= 5))
              {
                dealii::Tensor<1, dim> fe_v_value = (*fe_v_cell)[mag].value(j, q);
                cell_matrix(i, j) += fe_v_cell->JxW(q) * fe_v_cell->shape_value(i, q) * fe_v_value[component_ii[i] - 5];
              }
            }
            else
            {
              if (is_primitive[j] && (component_ii[j] >= 5))
              {
                dealii::Tensor<1, dim> fe_v_value = (*fe_v_cell)[mag].value(i, q);
                cell_matrix(i, j) += fe_v_cell->JxW(q) * fe_v_cell->shape_value(j, q) * fe_v_value[component_ii[j] - 5];
              }
              else if (!is_primitive[j])
                cell_matrix(i, j) += fe_v_cell->JxW(q) * fe_v_cell->shape_value(i, q) * fe_v_cell->shape_value(j, q);
            }
          }
        }
      }
      cell_rhs(i) += val;
    }
  }
  else
  {
    for (unsigned int q = 0; q < n_quadrature_points_cell; ++q)
    {
      for (unsigned int c = 0; c < Equations<equationsType, dim>::n_components; ++c)
        W_prev[q][c] = W_lin[c] = 0.;

      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        if (!is_primitive[i])
        {
          dealii::Tensor<1, dim> fe_v_value = (*fe_v_cell)[mag].value(i, q);
          for (unsigned int d = 0; d < dim; d++)
          {
            W_prev[q][5 + d] += prev_solution(dof_indices[i]) * fe_v_value[d];
            W_lin[5 + d] += lin_solution(dof_indices[i]) * fe_v_value[d];
          }
        }
        else
        {
          W_prev[q][component_ii[i]] += prev_solution(dof_indices[i]) * fe_v_cell->shape_value(i, q);
          W_lin[component_ii[i]] += lin_solution(dof_indices[i]) * fe_v_cell->shape_value(i, q);
        }
      }
      equations.compute_flux_matrix(W_lin, fluxes_old[q]);
    }

    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      double val = 0.;
      for (unsigned int q = 0; q < n_quadrature_points_cell; ++q)
      {
        if (is_primitive[i])
        {
          val += fe_v_cell->JxW(q) * W_prev[q][component_ii[i]] * fe_v_cell->shape_value(i, q);
          if (!basis_fn_is_constant[i])
            for (int d = 0; d < dim; d++)
              val += fe_v_cell->JxW(q) * parameters.time_step * fluxes_old[q][component_ii[i]][d] * fe_v_cell->shape_grad(i, q)[d];
        }
        else
        {
          dealii::Tensor<1, dim> fe_v_value = (*fe_v_cell)[mag].value(i, q);
          dealii::Tensor<2, dim> fe_v_grad = (*fe_v_cell)[mag].gradient(i, q);
          for (unsigned int d = 0; d < dim; d++)
          {
            val += fe_v_cell->JxW(q) * W_prev[q][5 + d] * fe_v_value[d];
            for (int e = 0; e < dim; e++)
              val += fe_v_cell->JxW(q) * parameters.time_step * fluxes_old[q][5 + d][e] * fe_v_grad[d][e];
          }
        }
      }
      cell_rhs(i) += val;
    }
  }
}

template <EquationsType equationsType, int dim>
void
Problem<equationsType, dim>::assemble_face_term(const unsigned int face_no, const FEFaceValuesBase<dim> &fe_v, const FEFaceValuesBase<dim> &fe_v_neighbor, const bool external_face, const unsigned int boundary_id, Vector<double>& cell_rhs)
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
      if (fe_v.get_fe().has_support_on_face(i, face_no) == true)
      {
        if (!is_primitive[i])
        {
          dealii::Tensor<1, dim> fe_v_value = (*fe_v_cell)[mag].value(i, q);
          for (int d = 0; d < dim; d++)
            Wplus_old[5 + d] += lin_solution(dof_indices[i]) * fe_v_value[d];
        }
        else
          Wplus_old[component_ii[i]] += lin_solution(dof_indices[i]) * fe_v.shape_value(i, q);

        if (!external_face)
        {
          if (!is_primitive[i])
          {
            dealii::Tensor<1, dim> fe_v_value_neighbor = fe_v_neighbor[mag].value(i, q);
            for (int d = 0; d < dim; d++)
              Wminus_old[5 + d] += lin_solution(dof_indices_neighbor[i]) * fe_v_value_neighbor[d];
          }
          else
            Wminus_old[component_ii[i]] += lin_solution(dof_indices_neighbor[i]) * fe_v_neighbor.shape_value(i, q);
        }
      }
    }

    if (external_face)
      boundary_conditions.bc_vector_value(boundary_id, fe_v.quadrature_point(q), Wminus_old, Wplus_old);

    // Once we have the states on both sides of the face, we need to calculate the numerical flux.
    equations.numerical_normal_flux(fe_v.normal_vector(q), Wplus_old, Wminus_old, normal_fluxes_old[q]);

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
    if (fe_v.get_fe().has_support_on_face(i, face_no) == true)
    {
      double val = 0.;
      for (unsigned int q = 0; q < n_quadrature_points_face; ++q)
      {
        if (!is_primitive[i])
        {
          dealii::Tensor<1, dim> fe_v_value = (*fe_v_cell)[mag].value(i, q);

          val += this->parameters.time_step * (normal_fluxes_old[q][5] * fe_v_value[0] + normal_fluxes_old[q][6] * fe_v_value[1] + normal_fluxes_old[q][7] * fe_v_value[2]) * fe_v.JxW(q);
        }
        else
          val += this->parameters.time_step * normal_fluxes_old[q][component_ii[i]] * fe_v.shape_value(i, q) * fe_v.JxW(q);

        if (std::isnan(val))
        {
          equations.numerical_normal_flux(fe_v.normal_vector(q), Wplus_old, Wminus_old, normal_fluxes_old[q]);
          std::cout << "isnan: " << val << std::endl;
          std::cout << "i: " << i << ", ci: " << (!is_primitive[i] ? 1 : fe_v.get_fe().system_to_component_index(i).first) << std::endl;
          std::cout << "point: " << fe_v.quadrature_point(q)[0] << ", " << fe_v.quadrature_point(q)[1] << ", " << fe_v.quadrature_point(q)[2] << std::endl;
          std::cout << "normal: " << fe_v.normal_vector(q)[0] << ", " << fe_v.normal_vector(q)[1] << ", " << fe_v.normal_vector(q)[2] << std::endl;
          for (int j = 0; j < Equations<equationsType, dim>::n_components; j++)
            std::cout << "W+ [" << j << "]: " << (double)Wplus_old[j] << ", W- [" << j << "]: " << (double)Wminus_old[j] << ", F [" << j << "]: " << (double)normal_fluxes_old[q][j] << std::endl;
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
    typename Equations<equationsType, dim>::Postprocessor postprocessor(equations);
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
  typename Equations<equationsType, dim>::Postprocessor postprocessor(equations);

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
  const std::string filename_base = "solution-" + Utilities::int_to_string(output_file_number, 3);

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
  std::string filename = "solution-" + Utilities::int_to_string(output_file_number, 3) + ".vtk";
  std::ofstream output(filename.c_str());
  data_out.write_vtk(output);
#endif

  ++output_file_number;
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::setup_initial_solution()
{
  prev_solution.reinit(locally_relevant_dofs, mpi_communicator);
  lin_solution.reinit(locally_relevant_dofs, mpi_communicator);
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

  lin_solution = prev_solution;
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
  bool previous_bad_step = false;

#ifdef OUTPUT_BASE
  output_base();
  exit(1);
#endif

  while (time < parameters.final_time)
  {
    // Some output.
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      std::cout << "T: " << time << std::endl;
      if (initial_step)
        std::cout << "   Number of active cells:       " << triangulation.n_active_cells() << std::endl << "   Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl << std::endl;
      std::cout << "   NonLin Res" << std::endl << "   _____________________________________" << std::endl;
    }

    double res_norm_prev[2] = { 0., 0. };
    double res_norm_initial;
    bool bad_step = false;
    for (int linStep = 0; linStep < this->parameters.newton_max_iterations; linStep++)
    {
      // Assemble
      system_rhs = 0;
      if (linStep == 0 && initial_step)
        system_matrix = 0;
      assemble_system((linStep == 0 && initial_step));

      // Output matrix & rhs if required.
      if (parameters.output_matrix)
        output_matrix(system_matrix, "matrix", time_step, linStep);
      if (parameters.output_rhs)
        output_vector(system_rhs, "rhs", time_step, linStep);

      // Solve
      solve(current_unlimited_solution, (linStep == 0 && initial_step));

      // Output solution if required
      if (parameters.output_solution)
        output_vector(current_unlimited_solution, "current_unlimited_solution", time_step, linStep);

      // Postprocess if required
      if (parameters.polynomial_order_dg > 0 && parameters.limit_in_nonlin_loop)
        postprocess();
      else
        current_limited_solution = current_unlimited_solution;

      if (initial_step || newton_damping > (1. - 1.e-8))
      {
        newton_update = current_limited_solution;
        newton_update -= lin_solution;
        lin_solution = current_limited_solution;
      }
      else
      {
        newton_update = current_limited_solution;
        newton_update -= lin_solution;

        newton_update *= newton_damping;
        lin_solution += newton_update;
      }

      // Get the residual norm.
      res_norm = newton_update.linfty_norm();
      if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
        std::cout << "\tLin step #" << linStep << ", error: " << res_norm << std::endl;

      if (linStep == 0)
        res_norm_initial = res_norm;

      if ((res_norm < parameters.newton_residual_norm_threshold) || ((linStep > 1) && ((std::abs(res_norm_prev[1] - res_norm) / res_norm) < 1.e-4)))
      {
        bad_step = false;
        break;
      }
      else if ((linStep == this->parameters.newton_max_iterations - 1) || ((linStep > 1) && (((res_norm - res_norm_prev[0]) / res_norm_prev[0]) > .1)))
      {
        parameters.cfl_constant *= this->parameters.decrease_factor;
        time -= parameters.time_step;
        parameters.time_step *= this->parameters.decrease_factor;
        time += parameters.time_step;
        if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
        {
          std::cout << "\t\tWorse CFL coefficient: " << parameters.cfl_constant << std::endl;
        }
        this->lin_solution = this->prev_solution;
        bad_step = true;
        break;
      }
      else
      {
        res_norm_prev[1] = res_norm_prev[0];
        res_norm_prev[0] = res_norm;
      }
    }

    if (!bad_step && !previous_bad_step)
    {
      if (!initial_step)
      {
        parameters.cfl_constant = std::min(this->parameters.initial_and_max_cfl_constant, parameters.cfl_constant * parameters.increase_factor);
        if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
          std::cout << "\t\tBetter CFL coefficient: " << parameters.cfl_constant << std::endl;
      }
      move_time_step_handle_outputs();
    }

    previous_bad_step = bad_step;
  }
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::move_time_step_handle_outputs()
{
  if (parameters.polynomial_order_dg > 0 && !parameters.limit_in_nonlin_loop)
  {
    postprocess();
    lin_solution = current_limited_solution;
  }

  this->prev_solution = this->lin_solution;

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
bool Problem<equationsType, dim>::is_periodic_boundary(int boundary_id) const
{
  for (int pb = 0; pb < this->parameters.periodic_boundaries.size(); pb++)
    if (this->parameters.periodic_boundaries[pb][0] == boundary_id || this->parameters.periodic_boundaries[pb][1] == boundary_id)
      return true;
  return false;
}

template class Problem<EquationsTypeMhd, 3>;

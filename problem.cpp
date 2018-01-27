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
  initial_quadrature(parameters.initial_quadrature_order),
  verbose_cout(std::cout, false),
  initial_step(true),
  last_output_time(0.), last_snapshot_time(0.), time(0.),
  time_step(0),
  mag(dim + 2)
{
  update_flags = update_values | update_JxW_values | update_gradients;
  face_update_flags = update_values | update_JxW_values | update_normal_vectors | update_q_points;
  neighbor_face_update_flags = update_values | update_q_points;
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
  // Number of DOFs per cell - we assume uniform polynomial order (we will only do h-adaptivity)
  const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

  // DOF indices both on the currently assembled element and the neighbor.
  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
  std::vector<types::global_dof_index> dof_indices_neighbor(dofs_per_cell);

  // What values we need for the assembly.
  const UpdateFlags update_flags = update_values;

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
    std::vector<unsigned int> lambda_indices_to_multiply[Equations<equationsType, dim>::n_components];
    std::vector<unsigned int> lambda_indices_to_multiply_all_B_components;

    cell->get_dof_indices(dof_indices);

    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      if (!is_primitive[i])
        lambda_indices_to_multiply_all_B_components.push_back(dof_indices[i]);
      else
      {
        if (!u_c_set[component_ii[i]])
        {
          u_c[component_ii[i]] = current_unlimited_solution(dof_indices[i]);
          u_c_set[component_ii[i]] = true;
        }
        else
          lambda_indices_to_multiply[component_ii[i]].push_back(dof_indices[i]);
      }
    }

    if (parameters.debug)
      std::cout << "cell: " << ++cell_count << " - center: " << cell->center() << ", values: " << u_c[0] << ", " << u_c[1] << ", " << u_c[2] << ", " << u_c[3] << ", " << u_c[4] << std::endl;

    double alpha_e[Equations<equationsType, dim>::n_components];
    for (int i = 0; i < Equations<equationsType, dim>::n_components; i++)
      alpha_e[i] = 1.;

    // For all vertices -> v_i
    for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i)
    {
      unsigned int v_i = cell->vertex_index(i);

      bool is_boundary_vertex = false;
      // For all faces, such that the face contains the vertex
      for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
      {
        if (cell->at_boundary(face_no) && !boundary_conditions.should_limit_this_boundary_id(cell->face(face_no)->boundary_id()))
        {
          TriaIterator<TriaAccessor<dim - 1, dim, dim> > face = cell->face(face_no);
          for (unsigned int face_i = 0; face_i < GeometryInfo<dim>::vertices_per_face; ++face_i)
          {
            if (face->vertex_index(face_i) == v_i)
            {
              is_boundary_vertex = true;
              break;
            }
          }
          if (is_boundary_vertex)
            break;
        }
      }
      if (is_boundary_vertex)
        continue;

      std::set<unsigned int> visited_faces;

      // (!!!) Find out u_i
      Vector<double> u_i(Equations<equationsType, dim>::n_components);
      VectorTools::point_value(dof_handler, current_unlimited_solution, cell->center() + (1. - NEGLIGIBLE) * (cell->vertex(i) - cell->center()), u_i);

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
      // For all faces, such that the face contains the vertex
      for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
      {
        if (visited_faces.find(cell->face_index(face_no)) != visited_faces.end())
          continue;

        visited_faces.insert(cell->face_index(face_no));

        // Look at the right neighbor (h-adaptivity)
        // (!!!) Assuming there is no division at this point (no adaptivity)
        if (cell->at_boundary(face_no))
          continue;

        TriaIterator<TriaAccessor<dim - 1, dim, dim> > face = cell->face(face_no);
        bool is_relevant_face = false;
        for (unsigned int face_i = 0; face_i < GeometryInfo<dim>::vertices_per_face; ++face_i)
        {
          if (face->vertex_index(face_i) == v_i)
            is_relevant_face = true;
        }
        if (is_relevant_face)
        {
          // Update u_i_min, u_i_max
          const typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face_no);
          neighbor->get_dof_indices(dof_indices_neighbor);

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

          // From the right neighbor, look at all faces, such that the face contains the vertex
          for (unsigned int neighbor_face_no = 0; neighbor_face_no < GeometryInfo<dim>::faces_per_cell; ++neighbor_face_no)
          {
            if (visited_faces.find(neighbor->face_index(neighbor_face_no)) != visited_faces.end())
              continue;

            visited_faces.insert(neighbor->face_index(neighbor_face_no));

            if (neighbor->at_boundary(neighbor_face_no))
              continue;

            // Look at the right neighbor's neighbor (h-adaptivity)
            // (!!!) Assuming there is no division at this point (no adaptivity)
            TriaIterator<TriaAccessor<dim - 1, dim, dim> > neighbor_face = neighbor->face(neighbor_face_no);
            bool is_neighbor_relevant_face = false;
            for (unsigned int neighbor_face_i = 0; neighbor_face_i < GeometryInfo<dim>::vertices_per_face; ++neighbor_face_i)
            {
              if (neighbor_face->vertex_index(neighbor_face_i) == v_i)
                is_neighbor_relevant_face = true;
            }
            if (is_neighbor_relevant_face)
            {
              // Update u_i_min, u_i_max
              const typename DoFHandler<dim>::cell_iterator neighbor_neighbor = neighbor->neighbor(neighbor_face_no);
              neighbor_neighbor->get_dof_indices(dof_indices_neighbor);

              bool u_i_extrema_set_neighbor[Equations<equationsType, dim>::n_components];
              for (int i = 0; i < Equations<equationsType, dim>::n_components; i++)
                u_i_extrema_set_neighbor[i] = false;
              for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                if (is_primitive[i])
                {
                  if (!u_i_extrema_set_neighbor[component_ii[i]])
                  {
                    u_i_min[component_ii[i]] = std::min(u_i_min[component_ii[i]], (double)current_unlimited_solution(dof_indices_neighbor[i]));
                    u_i_max[component_ii[i]] = std::max(u_i_max[component_ii[i]], (double)current_unlimited_solution(dof_indices_neighbor[i]));
                    u_i_extrema_set_neighbor[component_ii[i]] = true;
                  }
                }
              }
            }
          }
        }
      }

      if (!is_boundary_vertex)
      {
        // Based on u_i_min, u_i_max, u_i, get alpha_e
        for (int k = 0; k < Equations<equationsType, dim>::n_components; k++)
          if (std::abs((u_c[k] - u_i[k]) / u_c[k]) > NEGLIGIBLE)
          {
            alpha_e[k] = std::min(alpha_e[k], ((u_i[k] - u_c[k]) > 0.) ? std::min(1.0, (u_i_max[k] - u_c[k]) / (u_i[k] - u_c[k])) : std::min(1.0, (u_i_min[k] - u_c[k]) / (u_i[k] - u_c[k])));
            if (parameters.debug)
              std::cout << "\talpha_e[" << k << "]: " << alpha_e[k] << std::endl;
          }
      }
    }

    for (int k = 0; k < Equations<equationsType, dim>::n_components; k++)
      for (int i = 0; i < lambda_indices_to_multiply[k].size(); i++)
        current_limited_solution(lambda_indices_to_multiply[k][i]) *= alpha_e[k];

    double alpha_e_B = std::min(std::min(alpha_e[5], alpha_e[6]), alpha_e[7]);
    for (int i = 0; i < lambda_indices_to_multiply_all_B_components.size(); i++)
      current_limited_solution(lambda_indices_to_multiply_all_B_components[i]) *= alpha_e_B;
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

  // DOF indices both on the currently assembled element and the neighbor.
  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
  std::vector<types::global_dof_index> dof_indices_neighbor(dofs_per_cell);

  // DOF indices both on the currently assembled element and the neighbor.
  FEValues<dim> fe_v(mapping, fe, this->initial_step ? initial_quadrature : quadrature, update_flags);
  FEFaceValues<dim> fe_v_face(mapping, fe, face_quadrature, face_update_flags);
  FESubfaceValues<dim> fe_v_subface(mapping, fe, face_quadrature, face_update_flags);
  FEFaceValues<dim> fe_v_face_neighbor(mapping, fe, face_quadrature, neighbor_face_update_flags);
  FESubfaceValues<dim> fe_v_subface_neighbor(mapping, fe, face_quadrature, neighbor_face_update_flags);

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

    if (assemble_matrix)
      cell_matrix = 0;
    cell_rhs = 0;

    fe_v.reinit(cell);
    cell->get_dof_indices(dof_indices);

    if (parameters.debug)
      std::cout << "NEW CELL: " << ith_cell++ << std::endl;

    // Assemble the volumetric integrals.
    assemble_cell_term(fe_v, dof_indices, cell_matrix, cell_rhs, assemble_matrix);

    // Assemble the face integrals - ONLY if this is not the initial step (where we do the projection of the initial condition).
    if (!initial_step)
    {
      for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
      {
        // Boundary face - here we pass the boundary id
        if (cell->at_boundary(face_no))
        {
          bool is_periodic_boundary = false;
          for (int pb = 0; pb < this->parameters.periodic_boundaries.size(); pb++)
          {
            if (this->parameters.periodic_boundaries[pb][0] == cell->face(face_no)->boundary_id() || this->parameters.periodic_boundaries[pb][1] == cell->face(face_no)->boundary_id())
            {
              is_periodic_boundary = true;
              break;
            }
          }

          if (is_periodic_boundary)
          {
            const DealIIExtensions::FacePair<dim>&  face_pair = periodic_cell_map.find(std::make_pair(cell, face_no))->second;
            typename DoFHandler<dim>::active_cell_iterator neighbor(cell);
            neighbor = ((*(face_pair.cell[0])).active_cell_index() == (*cell).active_cell_index()) ? face_pair.cell[1] : face_pair.cell[0];
            const unsigned int neighbor_face = ((*(face_pair.cell[0])).active_cell_index() == (*cell).active_cell_index()) ? face_pair.face_idx[1] : face_pair.face_idx[0];

            neighbor->get_dof_indices(dof_indices_neighbor);

            fe_v_face.reinit(cell, face_no);
            fe_v_face_neighbor.reinit(neighbor, neighbor_face);

            assemble_face_term(face_no, fe_v_face, fe_v_face_neighbor, dof_indices, dof_indices_neighbor, false, cell->face(face_no)->boundary_id(), cell->face(face_no)->diameter(), cell_matrix, cell_rhs, assemble_matrix);
          }
          else
          {
            fe_v_face.reinit(cell, face_no);
            assemble_face_term(face_no, fe_v_face, fe_v_face, dof_indices, std::vector<types::global_dof_index>(), true, cell->face(face_no)->boundary_id(), cell->face(face_no)->diameter(), cell_matrix, cell_rhs, assemble_matrix);
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

              assemble_face_term(face_no, fe_v_subface, fe_v_face_neighbor, dof_indices, dof_indices_neighbor, false, numbers::invalid_unsigned_int, neighbor_child->face(neighbor2)->diameter(), cell_matrix, cell_rhs, assemble_matrix);
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

            assemble_face_term(face_no, fe_v_face, fe_v_face_neighbor, dof_indices, dof_indices_neighbor, false, numbers::invalid_unsigned_int, cell->face(face_no)->diameter(), cell_matrix, cell_rhs, assemble_matrix);
          }
          // Here the neighbor face fits exactly the current face of the current element, this is the 'easy' part.
          // This is the only face assembly case performed without adaptivity.
          else
          {
            const typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face_no);
            neighbor->get_dof_indices(dof_indices_neighbor);

            fe_v_face.reinit(cell, face_no);
            fe_v_face_neighbor.reinit(neighbor, cell->neighbor_of_neighbor(face_no));

            assemble_face_term(face_no, fe_v_face, fe_v_face_neighbor, dof_indices, dof_indices_neighbor, false, numbers::invalid_unsigned_int, cell->face(face_no)->diameter(), cell_matrix, cell_rhs, assemble_matrix);
          }
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
Problem<equationsType, dim>::assemble_cell_term(const FEValues<dim> &fe_v, const std::vector<types::global_dof_index>& dof_indices, FullMatrix<double>& cell_matrix, Vector<double>& cell_rhs, bool assemble_matrix)
{
  const unsigned int dofs_per_cell = fe_v.dofs_per_cell;
  const unsigned int n_q_points = fe_v.n_quadrature_points;

  Table<2, double> W_prev(n_q_points, Equations<equationsType, dim>::n_components);
  Table<2, double> W_lin(n_q_points, Equations<equationsType, dim>::n_components);

  if (initial_step)
  {
    std::vector<Vector<double> > initial_values(n_q_points, Vector<double>(Equations<equationsType, dim>::n_components));
    initial_condition.vector_value(fe_v.get_quadrature_points(), initial_values);
    for (unsigned int q = 0; q < n_q_points; ++q)
      for (unsigned int i = 0; i < equations.n_components; ++i)
        W_prev[q][i] = initial_values[q][i];
  }
  else
  {
#ifdef USE_HDIV_FOR_B
    const FEValuesExtractors::Vector mag(dim + 2);
#endif
    for (unsigned int q = 0; q < n_q_points; ++q)
    {
      for (unsigned int c = 0; c < Equations<equationsType, dim>::n_components; ++c)
      {
        W_prev[q][c] = 0.;
        W_lin[q][c] = 0.;
      }

      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
#ifdef USE_HDIV_FOR_B
        const unsigned int component_i = fe_v.get_fe().system_to_base_index(i).first.first;
        if (component_i == 1)
        {
          dealii::Tensor<1, dim> fe_v_value = fe_v[mag].value(i, q);

          for (unsigned int d = 0; d < dim; d++)
          {
            W_prev[q][5 + d] += prev_solution(dof_indices[i]) * fe_v_value[d];
            W_lin[q][5 + d] += lin_solution(dof_indices[i]) * fe_v_value[d];
          }
        }
        else
#endif
        {
          W_prev[q][component_ii[i]] += prev_solution(dof_indices[i]) * fe_v.shape_value(i, q);
          W_lin[q][component_ii[i]] += lin_solution(dof_indices[i]) * fe_v.shape_value(i, q);
        }
      }
    }
  }

  std::vector < std_cxx11::array <std_cxx11::array <double, dim>, Equations<equationsType, dim>::n_components > > flux_old(n_q_points);

  for (unsigned int q = 0; q < n_q_points; ++q)
  {
    if (!initial_step)
      equations.compute_flux_matrix(W_lin[q], flux_old[q]);

    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      cell_rhs(i) += fe_v.JxW(q) * W_prev[q][component_ii[i]] * fe_v.shape_value(i, q);
      if (!initial_step)
      {
        for (int d = 0; d < dim; d++)
          cell_rhs(i) += fe_v.JxW(q) * parameters.time_step * flux_old[q][component_ii[i]][d] * fe_v.shape_grad(i, q)[d];
      }

      if (assemble_matrix)
        for (unsigned int j = 0; j < dofs_per_cell; ++j)
        {
          if (component_ii[i] == component_ii[j])
            cell_matrix(i, j) += fe_v.JxW(q) * fe_v.shape_value(i, q) * fe_v.shape_value(j, q);
        }
    }
  }
}

template <EquationsType equationsType, int dim>
void
Problem<equationsType, dim>::assemble_face_term(const unsigned int face_no,
  const FEFaceValuesBase<dim> &fe_v,
  const FEFaceValuesBase<dim> &fe_v_neighbor,
  const std::vector<types::global_dof_index> &dof_indices,
  const std::vector<types::global_dof_index> &dof_indices_neighbor,
  const bool                   external_face,
  const unsigned int           boundary_id,
  const double                 face_diameter,
  FullMatrix<double>& cell_matrix, Vector<double>& cell_rhs, bool assemble_matrix)
{
  const unsigned int n_q_points = fe_v.n_quadrature_points;
  const unsigned int dofs_per_cell = fe_v.dofs_per_cell;

  Table<2, double> Wplus_old(n_q_points, Equations<equationsType, dim>::n_components), Wminus_old(n_q_points, Equations<equationsType, dim>::n_components);

  std::vector< std_cxx11::array < double, Equations<equationsType, dim>::n_components> > normal_fluxes_old(n_q_points);

  const FEValuesExtractors::Vector mag(dim + 2);

  if (parameters.debug)
    std::cout << "\tEdqe: " << face_no << std::endl;

  // This loop is preparation - calculate all states (Wplus on the current element side of the currently assembled face, Wminus on the other side).
  for (unsigned int q = 0; q < n_q_points; ++q)
  {
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      if (fe_v.get_fe().has_support_on_face(i, face_no) == true)
      {
#ifdef USE_HDIV_FOR_B
        const unsigned int component_i = fe_v.get_fe().system_to_base_index(i).first.first;
        if (component_i == 1)
        {
          dealii::Tensor<1, dim> fe_v_value = fe_v[mag].value(i, q);
          Wplus_old[q][5] += lin_solution(dof_indices[i]) * fe_v_value[0];
          Wplus_old[q][6] += lin_solution(dof_indices[i]) * fe_v_value[1];
          Wplus_old[q][7] += lin_solution(dof_indices[i]) * fe_v_value[2];
        }
        else
#endif
        {
          Wplus_old[q][component_ii[i]] += lin_solution(dof_indices[i]) * fe_v.shape_value(i, q);
        }

        if (!external_face)
        {
#ifdef USE_HDIV_FOR_B
          const unsigned int component_i_neighbor = fe_v_neighbor.get_fe().system_to_base_index(i).first.first;
          if (component_i_neighbor == 1)
          {
            Wminus_old[q][5] += lin_solution(dof_indices_neighbor[i]) * fe_v_neighbor[mag].value(i, q)[0];
            Wminus_old[q][6] += lin_solution(dof_indices_neighbor[i]) * fe_v_neighbor[mag].value(i, q)[1];
            Wminus_old[q][7] += lin_solution(dof_indices_neighbor[i]) * fe_v_neighbor[mag].value(i, q)[2];
          }
          else
#endif
          {
            Wminus_old[q][component_ii[i]] += lin_solution(dof_indices_neighbor[i]) * fe_v_neighbor.shape_value(i, q);
          }
        }
      }
    }

    if (external_face)
    {
      dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1> Wminus_old_q = Wminus_old[q];
      boundary_conditions.bc_vector_value(boundary_id, fe_v.quadrature_point(q), Wminus_old_q, Wplus_old[q]);
      for (unsigned int di = 0; di < this->equations.n_components; ++di)
        Wminus_old[q][di] = Wminus_old_q[di];
    }

    // Once we have the states on both sides of the face, we need to calculate the numerical flux.
    equations.numerical_normal_flux(fe_v.normal_vector(q), Wplus_old[q], Wminus_old[q], normal_fluxes_old[q]);

    // Some debugging outputs.
    if (parameters.debug)
    {
      std::cout << "\t\tpoint_i: " << q << std::endl;
      std::cout << "\t\tq: " << fe_v.quadrature_point(q) << ", n: " << fe_v.normal_vector(q)[0] << ", " << fe_v.normal_vector(q)[1] << ", " << fe_v.normal_vector(q)[2] << std::endl;
      std::cout << "\t\tWplus: ";
      for (unsigned int i = 0; i < Equations<equationsType, dim>::n_components; i++)
        std::cout << Wplus_old[q][i] << (i < 7 ? ", " : "");
      std::cout << std::endl;

      std::cout << "\t\tWminus: ";
      for (unsigned int i = 0; i < Equations<equationsType, dim>::n_components; i++)
        std::cout << Wminus_old[q][i] << (i < 7 ? ", " : "");
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
      const unsigned int component_i = fe_v.get_fe().system_to_base_index(i).first.first;
      double val = 0.;

      for (unsigned int q = 0; q < n_q_points; ++q)
      {
#ifdef USE_HDIV_FOR_B
        if (component_i == 1)
        {
          dealii::Tensor<1, dim> fe_v_value = fe_v[mag].value(i, q);

          val += this->parameters.time_step * (normal_fluxes_old[q][5] * fe_v_value[0] + normal_fluxes_old[q][6] * fe_v_value[1] + normal_fluxes_old[q][7] * fe_v_value[2]) * fe_v.JxW(q);
        }
        else
#endif
          val += this->parameters.time_step * normal_fluxes_old[q][component_ii[i]] * fe_v.shape_value(i, q) * fe_v.JxW(q);

        if (std::isnan(val))
        {
          equations.numerical_normal_flux(fe_v.normal_vector(q), Wplus_old[q], Wminus_old[q], normal_fluxes_old[q]);
          std::cout << "isnan: " << val << std::endl;
#ifdef USE_HDIV_FOR_B
          std::cout << "i: " << i << ", ci: " << (component_i == 1 ? 1 : fe_v.get_fe().system_to_component_index(i).first) << std::endl;
#else
          std::cout << "i: " << i << ", ci: " << fe_v.get_fe().system_to_component_index(i).first << std::endl;
#endif
          std::cout << "point: " << fe_v.quadrature_point(q)[0] << ", " << fe_v.quadrature_point(q)[1] << ", " << fe_v.quadrature_point(q)[2] << std::endl;
          std::cout << "normal: " << fe_v.normal_vector(q)[0] << ", " << fe_v.normal_vector(q)[1] << ", " << fe_v.normal_vector(q)[2] << std::endl;
          for (int j = 0; j < Equations<equationsType, dim>::n_components; j++)
            std::cout << "W+ [" << j << "]: " << (double)Wplus_old[q][j] << ", W- [" << j << "]: " << (double)Wminus_old[q][j] << ", F [" << j << "]: " << (double)normal_fluxes_old[q][j] << std::endl;
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

      if ((res_norm < parameters.newton_residual_norm_threshold) || ((linStep > 1) && ((std::abs(res_norm_prev[1] - res_norm) / res_norm) < 1.e-6)))
      {
        bad_step = false;
        break;
      }
      else if ((linStep == this->parameters.newton_max_iterations - 1) || ((linStep > 1) && (((res_norm - res_norm_prev[0]) / res_norm_prev[0]) > 1.)))
      {
        newton_damping *= this->parameters.decrease_factor;
        parameters.cfl_constant *= this->parameters.decrease_factor;
        time -= parameters.time_step;
        parameters.time_step *= this->parameters.decrease_factor;
        time += parameters.time_step;
        if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
        {
          std::cout << "\t\tWorse damping coefficient: " << newton_damping << std::endl;
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
        newton_damping = std::min(this->parameters.initial_and_max_newton_damping, newton_damping * parameters.increase_factor);
        parameters.cfl_constant = std::min(this->parameters.initial_and_max_cfl_constant, parameters.cfl_constant * parameters.increase_factor);
        if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
        {
          std::cout << "\t\tBetter damping coefficient: " << newton_damping << std::endl;
          std::cout << "\t\tBetter CFL coefficient: " << parameters.cfl_constant << std::endl;
        }
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

template class Problem<EquationsTypeMhd, 3>;

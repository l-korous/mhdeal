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
  // Creating the FE system - first spaces for density & momentum, then for the mag field, and lastly for the energy.
  fe(FE_DGT<dim>(parameters.polynomial_order_dg), 4,
    FE_RaviartThomas<dim>(1), 1,
    FE_DGT<dim>(parameters.polynomial_order_dg), 1),
  dof_handler(triangulation),
  quadrature(parameters.quadrature_order),
  face_quadrature(parameters.quadrature_order),
  verbose_cout(std::cout, false),
  initial_step(true),
  assemble_just_rhs(false)
{
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
  constraints.close();

#ifdef HAVE_MPI
  SparsityTools::distribute_sparsity_pattern(dsp, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);
#endif

  system_rhs.reinit(locally_owned_dofs, mpi_communicator);
  system_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);
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

  // DOF indices both on the currently assembled element and the neighbor.
  FEValues<dim> fe_v(mapping, fe, quadrature, update_flags);
  FEValues<dim> fe_v_neighbor(mapping, fe, quadrature, update_flags);

  // This is what we return.
  current_limited_solution = current_solution;

  int cell_count = 0;
  // Loop through all cells.
  for (typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
  {
    if (!cell->is_locally_owned())
      continue;

    bool u_c_set[5] = { false, false, false, false, false };
    double u_c[5];
    std::vector<unsigned int> lambda_indices_to_multiply[5];

    fe_v.reinit(cell);
    cell->get_dof_indices(dof_indices);

    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      const unsigned int component_i = fe_v.get_fe().system_to_base_index(i).first.first;
      if (component_i != 1)
      {
        const unsigned int component_ii = fe_v.get_fe().system_to_component_index(i).first;
        if (!u_c_set[component_ii == 7 ? 4 : component_ii])
        {
          u_c[component_ii == 7 ? 4 : component_ii] = current_limited_solution(dof_indices[i]);
          u_c_set[component_ii == 7 ? 4 : component_ii] = true;
        }
        else
        {
          lambda_indices_to_multiply[component_ii == 7 ? 4 : component_ii].push_back(dof_indices[i]);
        }
      }
    }

    if (parameters.debug_limiter)
      std::cout << "cell: " << ++cell_count << " - center: " << cell->center() << ", values: " << u_c[0] << ", " << u_c[1] << ", " << u_c[2] << ", " << u_c[3] << ", " << u_c[4] << std::endl;

    double alpha_e[5] = { 1., 1., 1., 1., 1. };

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
      double u_i[5];
      Vector<double> vec_to_retrieve_value(8);
      VectorTools::point_value(dof_handler, current_limited_solution, cell->center() + (1. - 1.e-12) * (cell->vertex(i) - cell->center()), vec_to_retrieve_value);
      u_i[0] = vec_to_retrieve_value[0];
      u_i[1] = vec_to_retrieve_value[1];
      u_i[2] = vec_to_retrieve_value[2];
      u_i[3] = vec_to_retrieve_value[3];
      u_i[4] = vec_to_retrieve_value[7];

      if (this->parameters.debug_limiter)
        std::cout << "\tv_i: " << cell->vertex(i) << ", values: " << u_i[0] << ", " << u_i[1] << ", " << u_i[2] << ", " << u_i[3] << ", " << u_i[4] << std::endl;
      // Init u_i_min, u_i_max
      double u_i_min[5];
      double u_i_max[5];
      for (int k = 0; k < 5; k++)
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
          fe_v_neighbor.reinit(neighbor);
          neighbor->get_dof_indices(dof_indices_neighbor);

          bool u_i_extrema_set[5] = { false, false, false, false, false };
          for (unsigned int dof = 0; dof < dofs_per_cell; ++dof)
          {
            const unsigned int component_i = fe_v_neighbor.get_fe().system_to_base_index(dof).first.first;
            if (component_i != 1)
            {
              const unsigned int component_ii = fe_v_neighbor.get_fe().system_to_component_index(dof).first;
              if (!u_i_extrema_set[component_ii == 7 ? 4 : component_ii])
              {
                if (this->parameters.debug_limiter)
                {
                  if ((double)current_limited_solution(dof_indices_neighbor[dof]) < u_i_min[component_ii == 7 ? 4 : component_ii])
                    std::cout << "\tdecreasing u_i_min to: " << (double)current_limited_solution(dof_indices_neighbor[dof]) << std::endl;
                  if ((double)current_limited_solution(dof_indices_neighbor[dof]) > u_i_max[component_ii == 7 ? 4 : component_ii])
                    std::cout << "\tincreasing u_i_max to: " << (double)current_limited_solution(dof_indices_neighbor[dof]) << std::endl;
                }
                u_i_min[component_ii == 7 ? 4 : component_ii] = std::min(u_i_min[component_ii == 7 ? 4 : component_ii], (double)current_limited_solution(dof_indices_neighbor[dof]));
                u_i_max[component_ii == 7 ? 4 : component_ii] = std::max(u_i_max[component_ii == 7 ? 4 : component_ii], (double)current_limited_solution(dof_indices_neighbor[dof]));
                u_i_extrema_set[component_ii == 7 ? 4 : component_ii] = true;
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
              fe_v_neighbor.reinit(neighbor_neighbor);
              neighbor_neighbor->get_dof_indices(dof_indices_neighbor);

              bool u_i_extrema_set[5] = { false, false, false, false, false };
              for (unsigned int dof = 0; dof < dofs_per_cell; ++dof)
              {
                const unsigned int component_i = fe_v_neighbor.get_fe().system_to_base_index(dof).first.first;
                if (component_i != 1)
                {
                  const unsigned int component_ii = fe_v_neighbor.get_fe().system_to_component_index(dof).first;
                  if (!u_i_extrema_set[component_ii == 7 ? 4 : component_ii])
                  {
                    u_i_min[component_ii == 7 ? 4 : component_ii] = std::min(u_i_min[component_ii == 7 ? 4 : component_ii], (double)current_limited_solution(dof_indices_neighbor[dof]));
                    u_i_max[component_ii == 7 ? 4 : component_ii] = std::max(u_i_max[component_ii == 7 ? 4 : component_ii], (double)current_limited_solution(dof_indices_neighbor[dof]));
                    u_i_extrema_set[component_ii == 7 ? 4 : component_ii] = true;
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
        for (int k = 0; k < 5; k++)
          if ((std::abs(u_c[k]) > 1e-8) && (std::abs((u_c[k] - u_i[k]) / u_c[k]) > 1e-6))
            alpha_e[k] = std::min(alpha_e[k], ((u_i[k] - u_c[k]) > 0.) ? std::min(1.0, (u_i_max[k] - u_c[k]) / (u_i[k] - u_c[k])) : std::min(1.0, (u_i_min[k] - u_c[k]) / (u_i[k] - u_c[k])));
      }
    }

    for (int k = 0; k < 5; k++)
      for (int i = 0; i < lambda_indices_to_multiply[k].size(); i++)
        current_limited_solution(lambda_indices_to_multiply[k][i]) *= alpha_e[k];
  }
}


template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::assemble_system()
{
  // Number of DOFs per cell - we assume uniform polynomial order (we will only do h-adaptivity)
  const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

  // DOF indices both on the currently assembled element and the neighbor.
  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
  std::vector<types::global_dof_index> dof_indices_neighbor(dofs_per_cell);

  // What values we need for the assembly.
  const UpdateFlags update_flags = update_values | update_q_points | update_JxW_values | update_gradients;
  const UpdateFlags face_update_flags = update_values | update_q_points | update_JxW_values | update_normal_vectors;
  const UpdateFlags neighbor_face_update_flags = update_q_points | update_values;

  // DOF indices both on the currently assembled element and the neighbor.
  FEValues<dim> fe_v(mapping, fe, quadrature, update_flags);
  FEFaceValues<dim> fe_v_face(mapping, fe, face_quadrature, face_update_flags);
  FESubfaceValues<dim> fe_v_subface(mapping, fe, face_quadrature, face_update_flags);
  FEFaceValues<dim> fe_v_face_neighbor(mapping, fe, face_quadrature, neighbor_face_update_flags);
  FESubfaceValues<dim> fe_v_subface_neighbor(mapping, fe, face_quadrature, neighbor_face_update_flags);

  // Local (cell) matrices and rhs - for the currently assembled element and the neighbor
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_matrix_neighbor(dofs_per_cell, dofs_per_cell);
  Vector<double> cell_rhs(dofs_per_cell);
  Vector<double> cell_rhs_neighbor(dofs_per_cell);

  // Loop through all cells.
  for (typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
  {
    // Only assemble what belongs to this process.
    if (!cell->is_locally_owned())
      continue;

    if (!assemble_just_rhs)
      cell_matrix = 0;
    cell_rhs = 0;

    fe_v.reinit(cell);
    cell->get_dof_indices(dof_indices);

    if (parameters.debug)
      std::cout << cell << std::endl;

    // Assemble the volumetric integrals.
    assemble_cell_term(fe_v, dof_indices, cell_matrix, cell_rhs);

    // Assemble the face integrals - ONLY if this is not the initial step (where we do the projection of the initial condition).
    if (!initial_step)
    {
      for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
      {
        if (!assemble_just_rhs)
          cell_matrix_neighbor = 0;
        cell_rhs_neighbor = 0;

        // Boundary face - here we pass the boundary id
        if (cell->at_boundary(face_no))
        {
          fe_v_face.reinit(cell, face_no);
          assemble_face_term(face_no, fe_v_face, fe_v_face, dof_indices, std::vector<types::global_dof_index>(), true, cell->face(face_no)->boundary_id(), cell->face(face_no)->diameter(), cell_matrix, cell_rhs, cell_matrix_neighbor, cell_rhs_neighbor);
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

              assemble_face_term(face_no, fe_v_subface, fe_v_face_neighbor, dof_indices, dof_indices_neighbor, false, numbers::invalid_unsigned_int, neighbor_child->face(neighbor2)->diameter(), cell_matrix, cell_rhs, cell_matrix_neighbor, cell_rhs_neighbor);

              constraints.distribute_local_to_global(cell_matrix_neighbor, cell_rhs_neighbor, dof_indices_neighbor, system_matrix, system_rhs);
              if (assemble_just_rhs)
                constraints.distribute_local_to_global(cell_rhs_neighbor, dof_indices_neighbor, system_rhs);
              else
                constraints.distribute_local_to_global(cell_matrix_neighbor, cell_rhs_neighbor, dof_indices_neighbor, system_matrix, system_rhs);
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

            assemble_face_term(face_no, fe_v_face, fe_v_face_neighbor, dof_indices, dof_indices_neighbor, false, numbers::invalid_unsigned_int, cell->face(face_no)->diameter(), cell_matrix, cell_rhs, cell_matrix_neighbor, cell_rhs_neighbor);

            if (assemble_just_rhs)
              constraints.distribute_local_to_global(cell_rhs_neighbor, dof_indices_neighbor, system_rhs);
            else
              constraints.distribute_local_to_global(cell_matrix_neighbor, cell_rhs_neighbor, dof_indices_neighbor, system_matrix, system_rhs);
          }
          // Here the neighbor face fits exactly the current face of the current element, this is the 'easy' part.
          // This is the only face assembly case performed without adaptivity.
          else
          {
            const typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face_no);
            neighbor->get_dof_indices(dof_indices_neighbor);

            fe_v_face.reinit(cell, face_no);
            fe_v_face_neighbor.reinit(neighbor, cell->neighbor_of_neighbor(face_no));

            assemble_face_term(face_no, fe_v_face, fe_v_face_neighbor, dof_indices, dof_indices_neighbor, false, numbers::invalid_unsigned_int, cell->face(face_no)->diameter(), cell_matrix, cell_rhs, cell_matrix_neighbor, cell_rhs_neighbor);

            if (assemble_just_rhs)
              constraints.distribute_local_to_global(cell_rhs_neighbor, dof_indices_neighbor, system_rhs);
            else
              constraints.distribute_local_to_global(cell_matrix_neighbor, cell_rhs_neighbor, dof_indices_neighbor, system_matrix, system_rhs);
          }
        }
      }
    }

    if (assemble_just_rhs)
      constraints.distribute_local_to_global(cell_rhs, dof_indices, system_rhs);
    else
      constraints.distribute_local_to_global(cell_matrix, cell_rhs, dof_indices, system_matrix, system_rhs);
  }

  if (!assemble_just_rhs)
    system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
}

template <EquationsType equationsType, int dim>
void
Problem<equationsType, dim>::assemble_cell_term(const FEValues<dim> &fe_v, const std::vector<types::global_dof_index>& dof_indices, FullMatrix<double>& cell_matrix, Vector<double>& cell_rhs)
{
  const unsigned int dofs_per_cell = fe_v.dofs_per_cell;
  const unsigned int n_q_points = fe_v.n_quadrature_points;

  // This is for the explicit case.
  if (parameters.theta < 1.)
  {
    Table<2, double> W_old(n_q_points, Equations<equationsType, dim>::n_components);
    Table<3, double> grad_W_old(n_q_points, Equations<equationsType, dim>::n_components, dim);

    for (unsigned int q = 0; q < n_q_points; ++q)
    {
      for (unsigned int c = 0; c < Equations<equationsType, dim>::n_components; ++c)
      {
        W_old[q][c] = 0;
        if (parameters.needs_gradients)
        {
          for (unsigned int d = 0; d < dim; ++d)
            grad_W_old[q][c][d] = 0;
        }
      }
    }

    const FEValuesExtractors::Vector mag(dim + 1);

    if (initial_step)
    {
      std::vector<Vector<double> > initial_values(n_q_points, Vector<double>(Equations<equationsType, dim>::n_components));
      initial_condition.vector_value(fe_v.get_quadrature_points(), initial_values);
      for (unsigned int q = 0; q < n_q_points; ++q)
        for (unsigned int i = 0; i < equations.n_components; ++i)
          W_old[q][i] = initial_values[q][i];
    }
    else
    {
      for (unsigned int q = 0; q < n_q_points; ++q)
      {
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          const unsigned int component_i = fe_v.get_fe().system_to_base_index(i).first.first;
          // component_i == 1 means that this is in fact the vector-valued FE space for the magnetic field and we need to calculate the value for all three components of this vector field together.
          if (component_i == 1)
          {
            W_old[q][4] += old_solution(dof_indices[i]) * fe_v[mag].value(i, q)[0];
            W_old[q][5] += old_solution(dof_indices[i]) * fe_v[mag].value(i, q)[1];
            W_old[q][6] += old_solution(dof_indices[i]) * fe_v[mag].value(i, q)[2];
          }
          // For the other components (spaces), we go by each component.
          else
          {
            const unsigned int component_ii = fe_v.get_fe().system_to_component_index(i).first;
            W_old[q][component_ii] += old_solution(dof_indices[i]) * fe_v.shape_value_component(i, q, component_ii);
          }

          if (parameters.theta > 0. && parameters.needs_gradients)
          {
            for (unsigned int d = 0; d < dim; d++)
              grad_W_old[q][component_i][d] += old_solution(dof_indices[i]) * fe_v.shape_grad_component(i, q, component_i)[d];
          }
        }
      }
    }

    std::vector < std_cxx11::array <std_cxx11::array <double, dim>, Equations<equationsType, dim>::n_components > > flux_old(n_q_points);
    std::vector < std_cxx11::array< double, Equations<equationsType, dim>::n_components> > forcing_old(n_q_points);
    // LK: This is for e.g. stabilization for Euler
    std::vector < std_cxx11::array <std_cxx11::array <double, dim>, Equations<equationsType, dim>::n_components > > jacobian_addition_old(n_q_points);

    for (unsigned int q = 0; q < n_q_points; ++q)
    {
      equations.compute_flux_matrix(W_old[q], flux_old[q]);

      if (parameters.debug)
      {
        std::cout << "point_i: " << q << std::endl;
        std::cout << "q: " << fe_v.quadrature_point(q) << ", n: " << fe_v.quadrature_point(q)[0] << ", " << fe_v.quadrature_point(q)[1] << ", " << fe_v.quadrature_point(q)[2] << std::endl;
        std::cout << "W: ";
        for (unsigned int i = 0; i < 8; i++)
          std::cout << W_old[q][i] << (i < 7 ? ", " : "");
        std::cout << std::endl;

        std::cout << "F[x]: ";
        for (unsigned int i = 0; i < 8; i++)
          std::cout << flux_old[q][i][0] << (i < 7 ? ", " : "");
        std::cout << std::endl;

        std::cout << "F[y]: ";
        for (unsigned int i = 0; i < 8; i++)
          std::cout << flux_old[q][i][1] << (i < 7 ? ", " : "");
        std::cout << std::endl;

        std::cout << "F[z]: ";
        for (unsigned int i = 0; i < 8; i++)
          std::cout << flux_old[q][i][2] << (i < 7 ? ", " : "");
        std::cout << std::endl;
      }
      equations.compute_forcing_vector(W_old[q], forcing_old[q]);
      if (parameters.needs_gradients)
        equations.compute_jacobian_addition(fe_v.get_cell()->diameter(), grad_W_old[q], jacobian_addition_old[q]);
    }

    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      const unsigned int component_i = fe_v.get_fe().system_to_base_index(i).first.first;

      double val = 0;

      for (unsigned int q = 0; q < n_q_points; ++q)
      {
        // component_i == 1 means that this is in fact the vector-valued FE space for the magnetic field and we need to calculate the value for all three components of this vector field together.
        if (component_i == 1)
        {
          if (parameters.is_stationary == false)
            val -= (1.0 / parameters.time_step)
            * (W_old[q][4] * fe_v[mag].value(i, q)[0] + W_old[q][5] * fe_v[mag].value(i, q)[1] + W_old[q][6] * fe_v[mag].value(i, q)[2])
            * fe_v.JxW(q);

          if (parameters.debug && q == 0)
          {
            std::cout << "DOF: " << i << " - COMP: " << component_i << std::endl;
          }

          if (!initial_step)
          {
            for (unsigned int d = 0; d < dim; d++)
              val -= (1.0 - parameters.theta) * (flux_old[q][4][d] * fe_v[mag].gradient(i, q)[0][d] + flux_old[q][5][d] * fe_v[mag].gradient(i, q)[1][d] + flux_old[q][6][d] * fe_v[mag].gradient(i, q)[2][d]) * fe_v.JxW(q);

            if (parameters.needs_gradients)
              for (unsigned int d = 0; d < dim; d++)
                val += (1.0 - parameters.theta) * (jacobian_addition_old[q][4][d] * fe_v[mag].gradient(i, q)[0][d] + jacobian_addition_old[q][5][d] * fe_v[mag].gradient(i, q)[1][d] + jacobian_addition_old[q][6][d] * fe_v[mag].gradient(i, q)[2][d]) * fe_v.JxW(q);

            if (parameters.needs_forcing)
              val -= (1.0 - parameters.theta) * (forcing_old[q][4] * fe_v[mag].value(i, q)[0] + forcing_old[q][5] * fe_v[mag].value(i, q)[1] + forcing_old[q][6] * fe_v[mag].value(i, q)[2]) * fe_v.JxW(q);
          }
        }
        // For the other components (spaces), we go by each component.
        else
        {
          const unsigned int component_ii = fe_v.get_fe().system_to_component_index(i).first;
          if (parameters.is_stationary == false)
            val -= (1.0 / parameters.time_step) * W_old[q][component_ii] * fe_v.shape_value_component(i, q, component_ii) * fe_v.JxW(q);

          if (parameters.debug && q == 0)
          {
            std::cout << "DOF: " << i << " - COMP: " << component_i << " - SUB: " << component_ii << std::endl;
          }

          if (!initial_step)
          {
            for (unsigned int d = 0; d < dim; d++)
              val -= (1.0 - parameters.theta) * flux_old[q][component_ii][d] * fe_v.shape_grad_component(i, q, component_ii)[d] * fe_v.JxW(q);

            if (parameters.needs_gradients)
              for (unsigned int d = 0; d < dim; d++)
                val += (1.0 - parameters.theta) * jacobian_addition_old[q][component_i][d] * fe_v.shape_grad_component(i, q, component_i)[d] * fe_v.JxW(q);

            if (parameters.needs_forcing)
              val -= (1.0 - parameters.theta) * forcing_old[q][component_i] * fe_v.shape_value_component(i, q, component_i) * fe_v.JxW(q);
          }
        }
      }

      if (std::isnan(val))
      {
        std::cout << "isnan: " << val << std::endl;
        std::cout << "i: " << i << ", ci: " << (component_i == 1 ? 1 : fe_v.get_fe().system_to_component_index(i).first) << std::endl;
        std::cout << "point: " << fe_v.quadrature_point(0)[0] << ", " << fe_v.quadrature_point(0)[1] << ", " << fe_v.quadrature_point(0)[2] << std::endl;
        for (int j = 0; j < 8; j++)
          std::cout << "W [" << j << "]: " << W_old[0][j] << ", F [" << j << "]: " << flux_old[0][j][0] << ", " << flux_old[0][j][1] << ", " << flux_old[0][j][2] << std::endl;
      }

      cell_rhs(i) -= val;
    }
  }

  // Now this is for the implicit case. This is different from the face term, as the time derivative needs to be there in both explicit and implicit cases.
  // So we do all the preparations, but only add the time-derivative contribution to the assembly if we in fact have the explicit case.
  {
    Table<2, Sacado::Fad::DFad<double> > W(n_q_points, Equations<equationsType, dim>::n_components);
    Table<3, Sacado::Fad::DFad<double> > grad_W(n_q_points, Equations<equationsType, dim>::n_components, dim);
    std::vector<double> residual_derivatives(dofs_per_cell);

    std::vector<Sacado::Fad::DFad<double> > independent_local_dof_values(dofs_per_cell);
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      independent_local_dof_values[i] = current_solution(dof_indices[i]);

    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      independent_local_dof_values[i].diff(i, dofs_per_cell);

    for (unsigned int q = 0; q < n_q_points; ++q)
    {
      for (unsigned int c = 0; c < Equations<equationsType, dim>::n_components; ++c)
      {
        W[q][c] = 0;
        if (parameters.needs_gradients)
        {
          for (unsigned int d = 0; d < dim; ++d)
            grad_W[q][c][d] = 0;
        }
      }
    }

    const FEValuesExtractors::Vector mag(dim + 1);

    for (unsigned int q = 0; q < n_q_points; ++q)
    {
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        const unsigned int component_i = fe_v.get_fe().system_to_base_index(i).first.first;

        // component_i == 1 means that this is in fact the vector-valued FE space for the magnetic field and we need to calculate the value for all three components of this vector field together.
        if (component_i == 1)
        {
          W[q][4] += independent_local_dof_values[i] * fe_v[mag].value(i, q)[0];
          W[q][5] += independent_local_dof_values[i] * fe_v[mag].value(i, q)[1];
          W[q][6] += independent_local_dof_values[i] * fe_v[mag].value(i, q)[2];
        }
        // For the other components (spaces), we go by each component.
        else
        {
          const unsigned int component_ii = fe_v.get_fe().system_to_component_index(i).first;
          W[q][component_ii] += independent_local_dof_values[i] * fe_v.shape_value_component(i, q, component_ii);
        }

        if (parameters.theta > 0. && parameters.needs_gradients)
        {
          for (unsigned int d = 0; d < dim; d++)
            grad_W[q][component_i][d] += independent_local_dof_values[i] * fe_v.shape_grad_component(i, q, component_i)[d];
        }
      }
    }

    std::vector < std_cxx11::array <std_cxx11::array <Sacado::Fad::DFad<double>, dim>, Equations<equationsType, dim>::n_components > > flux(n_q_points);
    std::vector < std_cxx11::array< Sacado::Fad::DFad<double>, Equations<equationsType, dim>::n_components> > forcing(n_q_points);
    // This is for e.g. stabilization for Euler, it is any addition to the jacobian that is not a part of the flux.
    std::vector < std_cxx11::array <std_cxx11::array <Sacado::Fad::DFad<double>, dim>, Equations<equationsType, dim>::n_components > > jacobian_addition(n_q_points);

    if (parameters.theta > 0.)
    {
      for (unsigned int q = 0; q < n_q_points; ++q)
      {
        equations.compute_flux_matrix(W[q], flux[q]);
        equations.compute_forcing_vector(W[q], forcing[q]);
        if (parameters.needs_gradients)
          equations.compute_jacobian_addition(fe_v.get_cell()->diameter(), grad_W[q], jacobian_addition[q]);
      }
    }

    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      const unsigned int component_i = fe_v.get_fe().system_to_base_index(i).first.first;
      Sacado::Fad::DFad<double> R_i = 0;

      for (unsigned int q = 0; q < n_q_points; ++q)
      {
        // component_i == 1 means that this is in fact the vector-valued FE space for the magnetic field and we need to calculate the value for all three components of this vector field together.
        if (component_i == 1)
        {
          if (parameters.is_stationary == false)
            R_i += (1.0 / parameters.time_step)
            * (W[q][4] * fe_v[mag].value(i, q)[0] + W[q][5] * fe_v[mag].value(i, q)[1] + W[q][6] * fe_v[mag].value(i, q)[2])
            * fe_v.JxW(q);
        }
        else
        {
          const unsigned int component_ii = fe_v.get_fe().system_to_component_index(i).first;
          if (parameters.is_stationary == false)
            R_i += (1.0 / parameters.time_step) * W[q][component_ii] * fe_v.shape_value_component(i, q, component_ii) * fe_v.JxW(q);
        }

        if (parameters.theta > 0. && !initial_step) {
          // component_i == 1 means that this is in fact the vector-valued FE space for the magnetic field and we need to calculate the value for all three components of this vector field together.
          if (component_i == 1)
          {
            for (unsigned int d = 0; d < dim; d++)
              R_i -= parameters.theta * (flux[q][4][d] * fe_v[mag].gradient(i, q)[0][d] + flux[q][5][d] * fe_v[mag].gradient(i, q)[1][d] + flux[q][6][d] * fe_v[mag].gradient(i, q)[2][d]) * fe_v.JxW(q);

            if (parameters.needs_gradients)
              for (unsigned int d = 0; d < dim; d++)
                R_i += parameters.theta * (jacobian_addition[q][4][d] * fe_v[mag].gradient(i, q)[0][d] + jacobian_addition[q][5][d] * fe_v[mag].gradient(i, q)[1][d] + jacobian_addition[q][6][d] * fe_v[mag].gradient(i, q)[2][d]) * fe_v.JxW(q);

            if (parameters.needs_forcing)
              R_i -= parameters.theta * (forcing[q][4] * fe_v[mag].value(i, q)[0] + forcing[q][5] * fe_v[mag].value(i, q)[1] + forcing[q][6] * fe_v[mag].value(i, q)[2]) * fe_v.JxW(q);
          }
          // For the other components (spaces), we go by each component.
          else
          {
            for (unsigned int d = 0; d < dim; d++)
              R_i -= parameters.theta * flux[q][component_i][d] * fe_v.shape_grad_component(i, q, component_i)[d] * fe_v.JxW(q);

            if (parameters.needs_gradients)
              for (unsigned int d = 0; d < dim; d++)
                R_i += parameters.theta * jacobian_addition[q][component_i][d] * fe_v.shape_grad_component(i, q, component_i)[d] * fe_v.JxW(q);

            if (parameters.needs_forcing)
              R_i -= parameters.theta * forcing[q][component_i] * fe_v.shape_value_component(i, q, component_i) * fe_v.JxW(q);
          }
        }
      }

      if (!assemble_just_rhs)
        for (unsigned int k = 0; k < dofs_per_cell; ++k)
          cell_matrix(i, k) += R_i.fastAccessDx(k);

      if (std::isnan(R_i.val()))
      {
        std::cout << "isnan: " << R_i.val() << std::endl;
        std::cout << "i: " << i << ", ci: " << (component_i == 1 ? 1 : fe_v.get_fe().system_to_component_index(i).first) << std::endl;
        std::cout << "point: " << fe_v.quadrature_point(0)[0] << ", " << fe_v.quadrature_point(0)[1] << ", " << fe_v.quadrature_point(0)[2] << std::endl;
        for (int j = 0; j < 8; j++)
          std::cout << "W [" << j << "]: " << W[0][j] << ", F [" << j << "]: " << flux[0][j][0] << ", " << flux[0][j][1] << ", " << flux[0][j][2] << std::endl;
      }

      cell_rhs(i) -= R_i.val();
    }
  }
}

template <EquationsType equationsType, int dim>
void
Problem<equationsType, dim>::assemble_face_term(const unsigned int           face_no,
  const FEFaceValuesBase<dim> &fe_v,
  const FEFaceValuesBase<dim> &fe_v_neighbor,
  const std::vector<types::global_dof_index> &dof_indices,
  const std::vector<types::global_dof_index> &dof_indices_neighbor,
  const bool                   external_face,
  const unsigned int           boundary_id,
  const double                 face_diameter,
  FullMatrix<double>& cell_matrix, Vector<double>& cell_rhs, FullMatrix<double>& cell_matrix_neighbor, Vector<double>& cell_rhs_neighbor)
{
  const unsigned int n_q_points = fe_v.n_quadrature_points;
  const unsigned int dofs_per_cell = fe_v.dofs_per_cell;

  // This is for the explicit case.
  if (parameters.theta < 1.)
  {
    Table<2, double> Wplus_old(n_q_points, Equations<equationsType, dim>::n_components), Wminus_old(n_q_points, Equations<equationsType, dim>::n_components);

    std::vector< std_cxx11::array < double, Equations<equationsType, dim>::n_components> > normal_fluxes_old(n_q_points);

    const FEValuesExtractors::Vector mag(dim + 1);

    if (parameters.debug)
      std::cout << "edqe: " << face_no << std::endl;

    // This loop is preparation - calculate all states (Wplus on the current element side of the currently assembled face, Wminus on the other side).
    for (unsigned int q = 0; q < n_q_points; ++q)
    {
      dealii::internal::TableBaseAccessors::Accessor<2, double, false, 1> Wminus_old_q = Wminus_old[q];
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        const unsigned int component_i = fe_v.get_fe().system_to_base_index(i).first.first;

        // component_i == 1 means that this is in fact the vector-valued FE space for the magnetic field and we need to calculate the value for all three components of this vector field together.
        if (component_i == 1)
        {
          Wplus_old[q][4] += old_solution(dof_indices[i]) * fe_v[mag].value(i, q)[0];
          Wplus_old[q][5] += old_solution(dof_indices[i]) * fe_v[mag].value(i, q)[1];
          Wplus_old[q][6] += old_solution(dof_indices[i]) * fe_v[mag].value(i, q)[2];
        }
        // For the other components (spaces), we go by each component.
        else
        {
          const unsigned int component_ii = fe_v.get_fe().system_to_component_index(i).first;
          Wplus_old[q][component_ii] += old_solution(dof_indices[i]) * fe_v.shape_value_component(i, q, component_ii);
        }

        if (!external_face)
        {
          const unsigned int component_i_neighbor = fe_v_neighbor.get_fe().system_to_base_index(i).first.first;

          // component_i_neighbor == 1 means that this is in fact the vector-valued FE space for the magnetic field and we need to calculate the value for all three components of this vector field together.
          if (component_i_neighbor == 1)
          {
            Wminus_old[q][4] += old_solution(dof_indices_neighbor[i]) * fe_v_neighbor[mag].value(i, q)[0];
            Wminus_old[q][5] += old_solution(dof_indices_neighbor[i]) * fe_v_neighbor[mag].value(i, q)[1];
            Wminus_old[q][6] += old_solution(dof_indices_neighbor[i]) * fe_v_neighbor[mag].value(i, q)[2];
          }
          // For the other components (spaces), we go by each component.
          else
          {
            const unsigned int component_ii_neighbor = fe_v_neighbor.get_fe().system_to_component_index(i).first;
            Wminus_old[q][component_ii_neighbor] += old_solution(dof_indices_neighbor[i]) * fe_v_neighbor.shape_value_component(i, q, component_ii_neighbor);
          }
        }
      }

      // Wminus (state vector on the other side of the currently assembled face) on the boundary corresponds to the (Dirichlet) values.
      if (external_face)
      {
        // Boundary values handling.
        // If we have an external face, we need to get the actual values from the provided BoundaryConditions class using its method bc_vector_value().
        Vector<double> boundary_values(Equations<equationsType, dim>::n_components);
        boundary_conditions.bc_vector_value(boundary_id, fe_v.quadrature_point(q), Wminus_old_q, Wplus_old[q]);
      }


      // Once we have the states on both sides of the face, we need to calculate the numerical flux.
      equations.numerical_normal_flux(fe_v.normal_vector(q), Wplus_old[q], Wminus_old[q], normal_fluxes_old[q]);

      // Some debugging outputs.
      if (parameters.debug)
      {
        std::cout << "point_i: " << q << std::endl;
        std::cout << "q: " << fe_v.quadrature_point(q) << ", n: " << fe_v.normal_vector(q)[0] << ", " << fe_v.normal_vector(q)[1] << ", " << fe_v.normal_vector(q)[2] << std::endl;
        std::cout << "Wplus: ";
        for (unsigned int i = 0; i < 8; i++)
          std::cout << Wplus_old[q][i] << (i < 7 ? ", " : "");
        std::cout << std::endl;

        std::cout << "Wminus: ";
        for (unsigned int i = 0; i < 8; i++)
          std::cout << Wminus_old[q][i] << (i < 7 ? ", " : "");
        std::cout << std::endl;

        std::cout << "Num F: ";
        for (unsigned int i = 0; i < 8; i++)
          std::cout << normal_fluxes_old[q][i] << (i < 7 ? ", " : "");
        std::cout << std::endl;
      }
    }

    // The actual contributions to the right-hand-side (this is the explicit case, no contribution to the matrix here).
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      if (fe_v.get_fe().has_support_on_face(i, face_no) == true)
      {
        const unsigned int component_i = fe_v.get_fe().system_to_base_index(i).first.first;
        double val = 0.;

        for (unsigned int q = 0; q < n_q_points; ++q)
        {
          // component_i == 1 means that this is in fact the vector-valued FE space for the magnetic field and we need to calculate the value for all three components of this vector field together.
          if (component_i == 1)
          {
            val += (1.0 - parameters.theta)
              * (normal_fluxes_old[q][4] * fe_v[mag].value(i, q)[0] + normal_fluxes_old[q][5] * fe_v[mag].value(i, q)[1] + normal_fluxes_old[q][6] * fe_v[mag].value(i, q)[2])
              * fe_v.JxW(q);
          }
          // For the other components (spaces), we go by each component.
          else
          {
            const unsigned int component_ii = fe_v.get_fe().system_to_component_index(i).first;
            val += (1.0 - parameters.theta) * normal_fluxes_old[q][component_ii] * fe_v.shape_value_component(i, q, component_ii) * fe_v.JxW(q);
          }

          // Some debugging outputs.
          if (std::isnan(val))
          {
            equations.numerical_normal_flux(fe_v.normal_vector(q), Wplus_old[q], Wminus_old[q], normal_fluxes_old[q]);
            std::cout << "isnan: " << val << std::endl;
            std::cout << "i: " << i << ", ci: " << (component_i == 1 ? 1 : fe_v.get_fe().system_to_component_index(i).first) << std::endl;
            std::cout << "point: " << fe_v.quadrature_point(q)[0] << ", " << fe_v.quadrature_point(q)[1] << ", " << fe_v.quadrature_point(q)[2] << std::endl;
            std::cout << "normal: " << fe_v.normal_vector(q)[0] << ", " << fe_v.normal_vector(q)[1] << ", " << fe_v.normal_vector(q)[2] << std::endl;
            for (int j = 0; j < 8; j++)
              std::cout << "W+ [" << j << "]: " << Wplus_old[q][j] << ", W- [" << j << "]: " << Wminus_old[q][j] << ", F [" << j << "]: " << normal_fluxes_old[q][j] << std::endl;
          }
        }

        cell_rhs(i) -= val;
      }
    }
  }

  // Now this is for the implicit case.
  if (parameters.theta > 0.)
  {
    std::vector<Sacado::Fad::DFad<double> > independent_local_dof_values(dofs_per_cell), independent_neighbor_dof_values(external_face == false ? dofs_per_cell : 0);
    const unsigned int n_independent_variables = (external_face == false ? 2 * dofs_per_cell : dofs_per_cell);

    Table<2, Sacado::Fad::DFad<double> > Wplus(n_q_points, Equations<equationsType, dim>::n_components), Wminus(n_q_points, Equations<equationsType, dim>::n_components);

    std::vector< std_cxx11::array < Sacado::Fad::DFad<double>, Equations<equationsType, dim>::n_components> > normal_fluxes(n_q_points);

    for (unsigned int i = 0; i < dofs_per_cell; i++)
    {
      independent_local_dof_values[i] = current_solution(dof_indices[i]);
      independent_local_dof_values[i].diff(i, n_independent_variables);
    }

    // We need the neighbor values if this is an internal face.
    if (external_face == false)
    {
      for (unsigned int i = 0; i < dofs_per_cell; i++)
      {
        independent_neighbor_dof_values[i] = current_solution(dof_indices_neighbor[i]);
        independent_neighbor_dof_values[i].diff(i + dofs_per_cell, n_independent_variables);
      }
    }

    const FEValuesExtractors::Vector mag(dim + 1);

    for (unsigned int q = 0; q < n_q_points; ++q)
    {
      dealii::internal::TableBaseAccessors::Accessor<2, Sacado::Fad::DFad<double>, false, 1> Wminus_q = Wminus[q];
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        const unsigned int component_i = fe_v.get_fe().system_to_base_index(i).first.first;

        // component_i == 1 means that this is in fact the vector-valued FE space for the magnetic field and we need to calculate the value for all three components of this vector field together.
        if (component_i == 1)
        {
          Wplus[q][4] += independent_local_dof_values[i] * fe_v[mag].value(i, q)[0];
          Wplus[q][5] += independent_local_dof_values[i] * fe_v[mag].value(i, q)[1];
          Wplus[q][6] += independent_local_dof_values[i] * fe_v[mag].value(i, q)[2];
        }
        // For the other components (spaces), we go by each component.
        else
        {
          const unsigned int component_ii = fe_v.get_fe().system_to_component_index(i).first;
          Wplus[q][component_ii] += independent_local_dof_values[i] * fe_v.shape_value_component(i, q, component_ii);
        }

        if (!external_face)
        {
          const unsigned int component_i_neighbor = fe_v_neighbor.get_fe().system_to_base_index(i).first.first;

          // component_i_neighbor == 1 means that this is in fact the vector-valued FE space for the magnetic field and we need to calculate the value for all three components of this vector field together.
          if (component_i_neighbor == 1)
          {
            Wminus[q][4] += independent_neighbor_dof_values[i] * fe_v_neighbor[mag].value(i, q)[0];
            Wminus[q][5] += independent_neighbor_dof_values[i] * fe_v_neighbor[mag].value(i, q)[1];
            Wminus[q][6] += independent_neighbor_dof_values[i] * fe_v_neighbor[mag].value(i, q)[2];
          }
          // For the other components (spaces), we go by each component.
          else
          {
            const unsigned int component_ii_neighbor = fe_v_neighbor.get_fe().system_to_component_index(i).first;
            Wminus[q][component_ii_neighbor] += independent_neighbor_dof_values[i] * fe_v_neighbor.shape_value_component(i, q, component_ii_neighbor);
          }
        }
      }
      if (external_face)
      {
        // Boundary values handling.
        // If we have an external face, we need to get the actual values from the provided BoundaryConditions class using its method bc_vector_value().
        Vector<double> boundary_values(Equations<equationsType, dim>::n_components);
        boundary_conditions.bc_vector_value(boundary_id, fe_v.quadrature_point(q), Wminus_q, Wplus[q]);
      }

      equations.numerical_normal_flux(fe_v.normal_vector(q), Wplus[q], Wminus[q], normal_fluxes[q]);
    }

    std::vector<double> residual_derivatives(dofs_per_cell);

    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      if (fe_v.get_fe().has_support_on_face(i, face_no) == true)
      {
        const unsigned int component_i = fe_v.get_fe().system_to_base_index(i).first.first;
        Sacado::Fad::DFad<double> R_i = 0;

        for (unsigned int q = 0; q < n_q_points; ++q)
        {
          // component_i == 1 means that this is in fact the vector-valued FE space for the magnetic field and we need to calculate the value for all three components of this vector field together.
          if (component_i == 1)
          {
            R_i += (1.0 - parameters.theta)
              * (normal_fluxes[q][4] * fe_v[mag].value(i, q)[0] + normal_fluxes[q][5] * fe_v[mag].value(i, q)[1] + normal_fluxes[q][6] * fe_v[mag].value(i, q)[2])
              * fe_v.JxW(q);
          }
          // For the other components (spaces), we go by each component.
          else
          {
            const unsigned int component_ii = fe_v.get_fe().system_to_component_index(i).first;
            R_i += (1.0 - parameters.theta) * normal_fluxes[q][component_ii] * fe_v.shape_value_component(i, q, component_ii) * fe_v.JxW(q);
          }

          if (std::isnan(R_i.val()))
          {
            equations.numerical_normal_flux(fe_v.normal_vector(q), Wplus[q], Wminus[q], normal_fluxes[q]);
            std::cout << "isnan: " << R_i.val() << std::endl;
            std::cout << "i: " << i << ", ci: " << (component_i == 1 ? 1 : fe_v.get_fe().system_to_component_index(i).first) << std::endl;
            std::cout << "point: " << fe_v.quadrature_point(q)[0] << ", " << fe_v.quadrature_point(q)[1] << ", " << fe_v.quadrature_point(q)[2] << std::endl;
            std::cout << "normal: " << fe_v.normal_vector(q)[0] << ", " << fe_v.normal_vector(q)[1] << ", " << fe_v.normal_vector(q)[2] << std::endl;
            for (int j = 0; j < 8; j++)
              std::cout << "W+ [" << j << "]: " << Wplus[q][j] << ", W- [" << j << "]: " << Wminus[q][j] << ", F [" << j << "]: " << normal_fluxes[q][j] << std::endl;
          }
        }

        if (!assemble_just_rhs)
          for (unsigned int k = 0; k < dofs_per_cell; ++k)
            cell_matrix(i, k) += R_i.fastAccessDx(k);

        if (!assemble_just_rhs)
        {
          // We only add contribution to the matrix entries corresponding to a neighbor cell if there is any - and there is none there if we are dealing with an external face.
          if (external_face == false)
          {
            for (unsigned int k = 0; k < dofs_per_cell; ++k)
              cell_matrix_neighbor(i, k) += R_i.fastAccessDx(dofs_per_cell + k);
          }
        }
        cell_rhs(i) -= R_i.val();
      }
    }
  }
}

template <EquationsType equationsType, int dim>
void
Problem<equationsType, dim>::solve(TrilinosWrappers::MPI::Vector &newton_update)
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

    AztecOO solver;
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
    newton_update = completely_distributed_solution;
  }
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::output_results(const char* prefix) const
{
  typename Equations<equationsType, dim>::Postprocessor postprocessor(equations);

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);

  // Solution components.
  data_out.add_data_vector(current_solution, equations.component_names(), DataOut<dim>::type_dof_data, equations.component_interpretation());

  // Derived quantities.
  data_out.add_data_vector(current_solution, postprocessor);

#ifdef HAVE_MPI
  // Subdomains.
  Vector<float> subdomain(triangulation.n_active_cells());
  for (unsigned int i = 0; i < subdomain.size(); ++i)
    subdomain(i) = triangulation.locally_owned_subdomain();
  data_out.add_data_vector(subdomain, "subdomain");
#endif

  data_out.build_patches(2);

  static unsigned int output_file_number = 0;

#ifdef HAVE_MPI
  const std::string filename_base = std::string(prefix) + "solution-" + Utilities::int_to_string(output_file_number, 3);

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
  std::string filename = std::string(prefix) + "solution-" + Utilities::int_to_string(output_file_number, 3) + ".vtk";
  std::ofstream output(filename.c_str());
  data_out.write_vtk(output);
#endif

  ++output_file_number;
}


template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::setup_initial_solution()
{
  old_solution.reinit(locally_relevant_dofs, mpi_communicator);
  current_solution.reinit(locally_relevant_dofs, mpi_communicator);
  current_limited_solution.reinit(locally_owned_dofs, mpi_communicator);
  current_unlimited_solution.reinit(locally_relevant_dofs, mpi_communicator);
  newton_initial_guess.reinit(locally_owned_dofs, mpi_communicator);

  old_solution = 0;
  current_solution = old_solution;
  current_unlimited_solution = old_solution;
  newton_initial_guess = old_solution;
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
void Problem<equationsType, dim>::run()
{
  // Preparations.
  setup_system();
  setup_initial_solution();
  TrilinosWrappers::MPI::Vector newton_update(locally_relevant_dofs, mpi_communicator);

  // More preparations.
  double time = 0;
  int time_step = 0;
  double next_output = time + parameters.output_step;

  // Time loop.
  while (time < parameters.final_time)
  {
    // Some output.
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      std::cout << "T=" << time << std::endl << "   Number of active cells:       " << triangulation.n_active_cells() << std::endl
        << "   Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl << std::endl;
      std::cout << "   NonLin Res" << std::endl << "   _____________________________________" << std::endl;
    }

    // Preparation for the Newton loop.
    unsigned int newton_iter = 0;
    current_solution = newton_initial_guess;
    while (true)
    {
      if (!(initial_step && (newton_iter == 0)))
      {
        assemble_just_rhs = true;
        system_rhs = 0;
        assemble_system();
        // L2 norm of the residual - to check Newton convergence
        const double res_norm = system_rhs.l2_norm();
        if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
        {
          if (std::fabs(res_norm) < parameters.newton_residual_norm_threshold)
            std::printf("   %-16.3e (converged)\n\n", res_norm);
          else
            std::printf("   %-16.3e\n", res_norm);

          if (std::isnan(res_norm))
          {
            output_results();
            std::this_thread::sleep_for(std::chrono::milliseconds(10000));
            exit(1);
          }
        }

        // If we are below threshold for the L2 residual norm, break the loop and go to next time step
        if (std::fabs(res_norm) < parameters.newton_residual_norm_threshold)
          break;
      }

      // Assemble.
      current_solution = current_limited_solution;
      assemble_just_rhs = false;
      system_matrix = 0;
      system_rhs = 0;
      assemble_system();

      // Outputs of algebraic stuff (optional - possible to be set in the Parameters class).
      if (parameters.output_matrix)
        output_matrix(system_matrix, "matrix", time_step, newton_iter);
      if (parameters.output_rhs)
        output_vector(system_rhs, "rhs", time_step, newton_iter);
      if (parameters.output_solution)
        output_vector(newton_update, "newton_update", time_step, newton_iter);

      solve(newton_update);

      if (parameters.theta > 0.)
        newton_update *= parameters.newton_damping;

      // Update the unlimited solution, and make solution equal.
      current_solution += newton_update;
      current_unlimited_solution = current_solution;

      // Postprocess, and store into limited solution (keep current_solution intact)
      if (parameters.polynomial_order_dg > 0)
        postprocess();

      ++newton_iter;
      AssertThrow(newton_iter <= parameters.newton_max_iterations, ExcMessage("No convergence in nonlinear solver"));
    }

    // Make current_solution point to the limited one.
    current_solution = current_limited_solution;

    // Output solution vector optionally and move on to the next time step.
    if (parameters.output_solution)
      output_vector(newton_update, "current_solution", time_step);

    ++time_step;
    time += parameters.time_step;

    if (parameters.output_step < 0)
      output_results();
    else if (time >= next_output)
    {
      output_results();
      next_output += parameters.output_step;
    }

    newton_initial_guess = current_solution;
    old_solution = current_solution;
    initial_step = false;
  }
}

template class Problem<EquationsTypeMhd, 3>;

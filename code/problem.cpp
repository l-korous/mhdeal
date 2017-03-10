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
  fe(FE_DGQ<dim>(parameters.polynomial_order_dg), 4,
    FE_RaviartThomas<dim>(1), 1,
    FE_DGQ<dim>(parameters.polynomial_order_dg), 1),
  dof_handler(triangulation),
  quadrature(2 * std::max(parameters.polynomial_order_dg, parameters.polynomial_order_hdiv) + 3),
  face_quadrature(2 * std::max(parameters.polynomial_order_dg, parameters.polynomial_order_hdiv) + 3),
  verbose_cout(std::cout, false)
{
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::setup_system()
{
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

  system_rhs.reinit(locally_owned_dofs, mpi_communicator);

#ifdef HAVE_MPI
  SparsityTools::distribute_sparsity_pattern(dsp, dof_handler.n_locally_owned_dofs_per_processor(), mpi_communicator, locally_relevant_dofs);
#endif

  system_matrix.reinit(locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::assemble_system()
{
  const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
  std::vector<types::global_dof_index> dof_indices_neighbor(dofs_per_cell);

  const UpdateFlags update_flags = update_values | update_gradients | update_q_points | update_JxW_values;
  const UpdateFlags face_update_flags = update_values | update_q_points | update_JxW_values | update_normal_vectors;
  const UpdateFlags neighbor_face_update_flags = update_q_points | update_values;

  FEValues<dim> fe_v(mapping, fe, quadrature, update_flags);
  FEFaceValues<dim> fe_v_face(mapping, fe, face_quadrature, face_update_flags);
  FESubfaceValues<dim> fe_v_subface(mapping, fe, face_quadrature, face_update_flags);
  FEFaceValues<dim> fe_v_face_neighbor(mapping, fe, face_quadrature, neighbor_face_update_flags);
  FESubfaceValues<dim> fe_v_subface_neighbor(mapping, fe, face_quadrature, neighbor_face_update_flags);

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_matrix_neighbor(dofs_per_cell, dofs_per_cell);
  Vector<double> cell_rhs(dofs_per_cell);
  Vector<double> cell_rhs_neighbor(dofs_per_cell);

  for (typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
  {
    if (!cell->is_locally_owned())
      continue;

    cell_matrix = 0;
    cell_rhs = 0;

    fe_v.reinit(cell);
    cell->get_dof_indices(dof_indices);

    assemble_cell_term(fe_v, dof_indices, cell_matrix, cell_rhs);

    if (!parameters.initial_step)
    {
      for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
      {
        cell_matrix_neighbor = 0;
        cell_rhs_neighbor = 0;

        if (cell->at_boundary(face_no))
        {
          fe_v_face.reinit(cell, face_no);
          assemble_face_term(face_no, fe_v_face, fe_v_face, dof_indices, std::vector<types::global_dof_index>(), true, cell->face(face_no)->boundary_id(), cell->face(face_no)->diameter(), cell_matrix, cell_rhs, cell_matrix_neighbor, cell_rhs_neighbor);
        }
        else
        {
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
            }
          }
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

            constraints.distribute_local_to_global(cell_matrix_neighbor, cell_rhs_neighbor, dof_indices_neighbor, system_matrix, system_rhs);
          }
          else
          {
            const typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face_no);
            neighbor->get_dof_indices(dof_indices_neighbor);

            fe_v_face.reinit(cell, face_no);
            fe_v_face_neighbor.reinit(neighbor, cell->neighbor_of_neighbor(face_no));

            assemble_face_term(face_no, fe_v_face, fe_v_face_neighbor, dof_indices, dof_indices_neighbor, false, numbers::invalid_unsigned_int, cell->face(face_no)->diameter(), cell_matrix, cell_rhs, cell_matrix_neighbor, cell_rhs_neighbor);

            constraints.distribute_local_to_global(cell_matrix_neighbor, cell_rhs_neighbor, dof_indices_neighbor, system_matrix, system_rhs);
          }
        }
      }
    }

    constraints.distribute_local_to_global(cell_matrix, cell_rhs, dof_indices, system_matrix, system_rhs);
  }

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
        if (this->parameters.needs_gradients)
        {
          for (unsigned int d = 0; d < dim; ++d)
            grad_W_old[q][c][d] = 0;
        }
      }
    }

    const FEValuesExtractors::Vector mag(dim + 1);

    if (parameters.initial_step)
    {
      for (unsigned int q = 0; q < n_q_points; ++q)
      {
        Point<dim> p = fe_v.quadrature_point(q);
        for (unsigned int i = 0; i < this->equations.n_components; ++i)
          W_old[q][i] = initial_condition.value(p, i);
      }
    }
    else
    {
      for (unsigned int q = 0; q < n_q_points; ++q)
      {
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          const unsigned int component_i = fe_v.get_fe().system_to_base_index(i).first.first;

          if (component_i == 1)
          {
            W_old[q][4] += old_solution(dof_indices[i]) * fe_v[mag].value(i, q)[0];
            W_old[q][5] += old_solution(dof_indices[i]) * fe_v[mag].value(i, q)[1];
            W_old[q][6] += old_solution(dof_indices[i]) * fe_v[mag].value(i, q)[2];
          }
          else
          {
            const unsigned int component_ii = fe_v.get_fe().system_to_component_index(i).first;
            W_old[q][component_ii] += old_solution(dof_indices[i]) * fe_v.shape_value_component(i, q, component_ii);
          }

          if (parameters.theta > 0. && this->parameters.needs_gradients)
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
      equations.compute_forcing_vector(W_old[q], forcing_old[q]);
      if (this->parameters.needs_gradients)
        equations.compute_jacobian_addition(fe_v.get_cell()->diameter(), grad_W_old[q], jacobian_addition_old[q]);
    }

    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      const unsigned int component_i = fe_v.get_fe().system_to_base_index(i).first.first;

      double val = 0;

      for (unsigned int q = 0; q < n_q_points; ++q)
      {
        if (parameters.is_stationary == false)
        {
          if (component_i == 1)
          {
            val -= (1.0 / parameters.time_step)
              * (W_old[q][4] * fe_v[mag].value(i, q)[0] + W_old[q][5] * fe_v[mag].value(i, q)[1] + W_old[q][6] * fe_v[mag].value(i, q)[2])
              * fe_v.JxW(q);

            if (!parameters.initial_step)
              for (unsigned int d = 0; d < dim; d++)
                val -= (1.0 - parameters.theta) * (flux_old[q][4][d] * fe_v[mag].gradient(i, q)[0][d] + flux_old[q][5][d] * fe_v[mag].gradient(i, q)[1][d] + flux_old[q][6][d] * fe_v[mag].gradient(i, q)[2][d]) * fe_v.JxW(q);
          }
          else
          {
            const unsigned int component_ii = fe_v.get_fe().system_to_component_index(i).first;
            val -= (1.0 / parameters.time_step) * W_old[q][component_ii] * fe_v.shape_value_component(i, q, component_ii) * fe_v.JxW(q);

            if (!parameters.initial_step)
              for (unsigned int d = 0; d < dim; d++)
                val -= (1.0 - parameters.theta) * flux_old[q][component_ii][d] * fe_v.shape_grad_component(i, q, component_ii)[d] * fe_v.JxW(q);
          }
        }

        if (!parameters.initial_step)
          if (this->parameters.needs_gradients)
          {
            for (unsigned int d = 0; d < dim; d++)
              val += (1.0 - parameters.theta) * jacobian_addition_old[q][component_i][d] * fe_v.shape_grad_component(i, q, component_i)[d] * fe_v.JxW(q);
          }

        if (!parameters.initial_step)
          val -= (1.0 - parameters.theta) * forcing_old[q][component_i] * fe_v.shape_value_component(i, q, component_i) * fe_v.JxW(q);
      }

      cell_rhs(i) -= val;
    }
  }

  // Now this is for the implicit case. This is different from the face term, as the time derivative needs to be there in both explicit and implicit cases.
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
        if (this->parameters.needs_gradients)
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

        if (component_i == 1)
        {
          W[q][4] += independent_local_dof_values[i] * fe_v[mag].value(i, q)[0];
          W[q][5] += independent_local_dof_values[i] * fe_v[mag].value(i, q)[1];
          W[q][6] += independent_local_dof_values[i] * fe_v[mag].value(i, q)[2];
        }
        else
        {
          const unsigned int component_ii = fe_v.get_fe().system_to_component_index(i).first;
          W[q][component_ii] += independent_local_dof_values[i] * fe_v.shape_value_component(i, q, component_ii);
        }

        if (parameters.theta > 0. && this->parameters.needs_gradients)
        {
          for (unsigned int d = 0; d < dim; d++)
            grad_W[q][component_i][d] += independent_local_dof_values[i] * fe_v.shape_grad_component(i, q, component_i)[d];
        }
      }
    }

    std::vector < std_cxx11::array <std_cxx11::array <Sacado::Fad::DFad<double>, dim>, Equations<equationsType, dim>::n_components > > flux(n_q_points);
    std::vector < std_cxx11::array< Sacado::Fad::DFad<double>, Equations<equationsType, dim>::n_components> > forcing(n_q_points);
    // This is for e.g. stabilization for Euler
    std::vector < std_cxx11::array <std_cxx11::array <Sacado::Fad::DFad<double>, dim>, Equations<equationsType, dim>::n_components > > jacobian_addition(n_q_points);

    if (this->parameters.theta > 0.)
    {
      for (unsigned int q = 0; q < n_q_points; ++q)
      {
        equations.compute_flux_matrix(W[q], flux[q]);
        equations.compute_forcing_vector(W[q], forcing[q]);
        if (this->parameters.needs_gradients)
          equations.compute_jacobian_addition(fe_v.get_cell()->diameter(), grad_W[q], jacobian_addition[q]);
      }
    }

    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      const unsigned int component_i = fe_v.get_fe().system_to_base_index(i).first.first;
      Sacado::Fad::DFad<double> R_i = 0;

      for (unsigned int q = 0; q < n_q_points; ++q)
      {
        if (parameters.is_stationary == false)
        {
          if (component_i == 1)
          {
            R_i += (1.0 / parameters.time_step)
              * (W[q][4] * fe_v[mag].value(i, q)[0] + W[q][5] * fe_v[mag].value(i, q)[1] + W[q][6] * fe_v[mag].value(i, q)[2])
              * fe_v.JxW(q);
          }
          else
          {
            const unsigned int component_ii = fe_v.get_fe().system_to_component_index(i).first;
            R_i += (1.0 / parameters.time_step) * W[q][component_ii] * fe_v.shape_value_component(i, q, component_ii) * fe_v.JxW(q);
          }
        }

        if (this->parameters.theta > 0.) {
          for (unsigned int d = 0; d < dim; d++)
            R_i -= parameters.theta * flux[q][component_i][d] * fe_v.shape_grad_component(i, q, component_i)[d] * fe_v.JxW(q);

          if (this->parameters.needs_gradients)
          {
            for (unsigned int d = 0; d < dim; d++)
              R_i += parameters.theta * jacobian_addition[q][component_i][d] * fe_v.shape_grad_component(i, q, component_i)[d] * fe_v.JxW(q);
          }

          R_i -= parameters.theta * forcing[q][component_i] * fe_v.shape_value_component(i, q, component_i) * fe_v.JxW(q);
        }
      }

      for (unsigned int k = 0; k < dofs_per_cell; ++k)
        cell_matrix(i, k) += R_i.fastAccessDx(k);

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

  // Boundary values handling.
  std::vector<Vector<double> > boundary_values(n_q_points, Vector<double>(Equations<equationsType, dim>::n_components));
  if (external_face)
    boundary_conditions.bc_vector_value(boundary_id, fe_v.get_quadrature_points(), boundary_values);

  // This is for the explicit case.
  if (parameters.theta < 1.)
  {
    Table<2, double> Wplus_old(n_q_points, Equations<equationsType, dim>::n_components), Wminus_old(n_q_points, Equations<equationsType, dim>::n_components);

    std::vector< std_cxx11::array < double, Equations<equationsType, dim>::n_components> > normal_fluxes_old(n_q_points);

    const FEValuesExtractors::Vector mag(dim + 1);

    for (unsigned int q = 0; q < n_q_points; ++q)
    {
      if (parameters.initial_step)
      {
        Point<dim> p = fe_v.quadrature_point(q);
        for (unsigned int i = 0; i < this->equations.n_components; ++i)
          Wplus_old[q][i] = initial_condition.value(p, i);
        if (!external_face)
        {
          Point<dim> pn = fe_v_neighbor.quadrature_point(q);
          for (unsigned int i = 0; i < this->equations.n_components; ++i)
            Wminus_old[q][i] = initial_condition.value(pn, i);
        }
      }
      else
      {
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          const unsigned int component_i = fe_v.get_fe().system_to_base_index(i).first.first;

          if (component_i == 1)
          {
            Wplus_old[q][4] += old_solution(dof_indices[i]) * fe_v[mag].value(i, q)[0];
            Wplus_old[q][5] += old_solution(dof_indices[i]) * fe_v[mag].value(i, q)[1];
            Wplus_old[q][6] += old_solution(dof_indices[i]) * fe_v[mag].value(i, q)[2];
          }
          else
          {
            const unsigned int component_ii = fe_v.get_fe().system_to_component_index(i).first;
            Wplus_old[q][component_ii] += old_solution(dof_indices[i]) * fe_v.shape_value_component(i, q, component_ii);
          }

          if (!external_face)
          {
            const unsigned int component_i_neighbor = fe_v_neighbor.get_fe().system_to_base_index(i).first.first;

            if (component_i_neighbor == 1)
            {
              Wminus_old[q][4] += old_solution(dof_indices_neighbor[i]) * fe_v_neighbor[mag].value(i, q)[0];
              Wminus_old[q][5] += old_solution(dof_indices_neighbor[i]) * fe_v_neighbor[mag].value(i, q)[1];
              Wminus_old[q][6] += old_solution(dof_indices_neighbor[i]) * fe_v_neighbor[mag].value(i, q)[2];
            }
            else
            {
              const unsigned int component_ii_neighbor = fe_v_neighbor.get_fe().system_to_component_index(i).first;
              Wminus_old[q][component_ii_neighbor] += old_solution(dof_indices_neighbor[i]) * fe_v_neighbor.shape_value_component(i, q, component_ii_neighbor);
            }
          }
        }
      }
      if (external_face)
        equations.compute_Wminus(boundary_conditions.kind[boundary_id], fe_v.normal_vector(q), Wplus_old[q], boundary_values[q], Wminus_old[q]);

      equations.numerical_normal_flux(fe_v.normal_vector(q), Wplus_old[q], Wminus_old[q], normal_fluxes_old[q]);
    }

    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      if (fe_v.get_fe().has_support_on_face(i, face_no) == true)
      {
        const unsigned int component_i = fe_v.get_fe().system_to_base_index(i).first.first;
        double val = 0.;

        for (unsigned int q = 0; q < n_q_points; ++q)
        {
          if (component_i == 1)
          {
            val += (1.0 - parameters.theta)
              * (normal_fluxes_old[q][4] * fe_v[mag].value(i, q)[0] + normal_fluxes_old[q][5] * fe_v[mag].value(i, q)[1] + normal_fluxes_old[q][6] * fe_v[mag].value(i, q)[2])
              * fe_v.JxW(q);
          }
          else
          {
            const unsigned int component_ii = fe_v.get_fe().system_to_component_index(i).first;
            val += (1.0 - parameters.theta) * normal_fluxes_old[q][component_ii] * fe_v.shape_value_component(i, q, component_ii) * fe_v.JxW(q);
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
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        const unsigned int component_i = fe_v.get_fe().system_to_base_index(i).first.first;

        if (component_i == 1)
        {
          Wplus[q][4] += independent_local_dof_values[i] * fe_v[mag].value(i, q)[0];
          Wplus[q][5] += independent_local_dof_values[i] * fe_v[mag].value(i, q)[1];
          Wplus[q][6] += independent_local_dof_values[i] * fe_v[mag].value(i, q)[2];
        }
        else
        {
          const unsigned int component_ii = fe_v.get_fe().system_to_component_index(i).first;
          Wplus[q][component_ii] += independent_local_dof_values[i] * fe_v.shape_value_component(i, q, component_ii);
        }

        if (!external_face)
        {
          const unsigned int component_i_neighbor = fe_v_neighbor.get_fe().system_to_base_index(i).first.first;

          if (component_i_neighbor == 1)
          {
            Wminus[q][4] += independent_neighbor_dof_values[i] * fe_v_neighbor[mag].value(i, q)[0];
            Wminus[q][5] += independent_neighbor_dof_values[i] * fe_v_neighbor[mag].value(i, q)[1];
            Wminus[q][6] += independent_neighbor_dof_values[i] * fe_v_neighbor[mag].value(i, q)[2];
          }
          else
          {
            const unsigned int component_ii_neighbor = fe_v_neighbor.get_fe().system_to_component_index(i).first;
            Wminus[q][component_ii_neighbor] += independent_neighbor_dof_values[i] * fe_v_neighbor.shape_value_component(i, q, component_ii_neighbor);
          }
        }
      }
      if (external_face)
        equations.compute_Wminus(boundary_conditions.kind[boundary_id], fe_v.normal_vector(q), Wplus[q], boundary_values[q], Wminus[q]);

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
          if (component_i == 1)
          {
            R_i += (1.0 - parameters.theta)
              * (normal_fluxes[q][4] * fe_v[mag].value(i, q)[0] + normal_fluxes[q][5] * fe_v[mag].value(i, q)[1] + normal_fluxes[q][6] * fe_v[mag].value(i, q)[2])
              * fe_v.JxW(q);
          }
          else
          {
            const unsigned int component_ii = fe_v.get_fe().system_to_component_index(i).first;
            R_i += (1.0 - parameters.theta) * normal_fluxes[q][component_ii] * fe_v.shape_value_component(i, q, component_ii) * fe_v.JxW(q);
          }
        }

        for (unsigned int k = 0; k < dofs_per_cell; ++k)
          cell_matrix(i, k) += R_i.fastAccessDx(k);

        if (external_face == false)
        {
          for (unsigned int k = 0; k < dofs_per_cell; ++k)
            cell_matrix_neighbor(i, k) += R_i.fastAccessDx(dofs_per_cell + k);
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
void Problem<equationsType, dim>::output_results() const
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

  data_out.build_patches(std::min(4., std::ceil(std::pow(1000. / std::max(parameters.refinements[1], parameters.refinements[2]), .333))));

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
    data_out.write_visit_record(visit_master_output, filenames);
}
#else
  std::string filename = "solution-" + Utilities::int_to_string(output_file_number, 3) + ".vtk";
  std::ofstream output(filename.c_str());
  data_out.write_vtk(output);
#endif

  ++output_file_number;
}


template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::process_initial_condition()
{
  old_solution = 0;
  current_solution = old_solution;
  newton_initial_guess = old_solution;
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::run()
{
  setup_system();

  old_solution.reinit(locally_relevant_dofs, mpi_communicator);
  current_solution.reinit(locally_relevant_dofs, mpi_communicator);
  newton_initial_guess.reinit(locally_owned_dofs, mpi_communicator);
  TrilinosWrappers::MPI::Vector newton_update(locally_relevant_dofs, mpi_communicator);

  process_initial_condition();

  double time = 0;
  int time_step = 0;
  double next_output = time + parameters.output_step;

  while (time < parameters.final_time)
  {
    if (time > this->parameters.initialization_time) {
      this->parameters.time_step = this->parameters.time_step_after_initialization;
      this->parameters.theta = this->parameters.theta_after_initialization;
    }

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      std::cout << "T=" << time << std::endl << "   Number of active cells:       " << triangulation.n_active_cells() << std::endl
        << "   Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl << std::endl;
      std::cout << "   NonLin Res" << std::endl << "   _____________________________________" << std::endl;
    }

    unsigned int nonlin_iter = 0;

    current_solution = newton_initial_guess;
    while (true)
    {
      system_matrix = 0;
      system_rhs = 0;
      assemble_system();

      if (parameters.output_matrix)
      {
        std::ofstream m;
        std::stringstream ssm;
        ssm << time_step << "-" << nonlin_iter << "-" << Utilities::MPI::this_mpi_process(mpi_communicator) << ".matrix";
        m.open(ssm.str());
        system_matrix.print(m);
        m.close();
      }
      if (parameters.output_rhs)
      {
        std::ofstream r;
        std::stringstream ssr;
        ssr << time_step << "-" << nonlin_iter << "-" << Utilities::MPI::this_mpi_process(mpi_communicator) << ".rhs";
        r.open(ssr.str());
        system_rhs.print(r, 10, false, false);
        r.close();
      }
      const double res_norm = system_rhs.l2_norm();

      if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
        if (std::fabs(res_norm) < parameters.nonlinear_residual_norm_threshold)
          std::printf("   %-16.3e (converged)\n\n", res_norm);
        else
          std::printf("   %-16.3e\n", res_norm);
      }

      if (std::fabs(res_norm) < parameters.nonlinear_residual_norm_threshold)
        break;
      else
      {
        newton_update = 0;
        solve(newton_update);
        if (parameters.output_solution)
        {
          std::ofstream n;
          std::stringstream ssn;
          ssn << time_step << "-" << nonlin_iter << "-" << Utilities::MPI::this_mpi_process(mpi_communicator) << ".newton_update";
          n.open(ssn.str());
          newton_update.print(n, 10, false, false);
          n.close();
        }
        if (parameters.theta > 0.)
          newton_update *= parameters.newton_damping;

        current_solution += newton_update;

        if (std::isnan(res_norm))
        {
          output_results();
          exit(1);
        }
      }

      ++nonlin_iter;
      AssertThrow(nonlin_iter <= parameters.max_nonlinear_iterations, ExcMessage("No convergence in nonlinear solver"));
    }

    if (parameters.output_solution)
    {
      std::ofstream s;
      std::stringstream sss;
      sss << time_step << "-" << Utilities::MPI::this_mpi_process(mpi_communicator) << ".current_solution";
      s.open(sss.str());
      current_solution.print(s, 10, false, false);
      s.close();
    }

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
    this->parameters.initial_step = false;
  }
}

template class Problem<EquationsTypeEuler, 2>;
template class Problem<EquationsTypeEuler, 3>;

template class Problem<EquationsTypeMhd, 3>;

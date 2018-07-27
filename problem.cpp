#include "problem.h"

template <EquationsType equationsType, int dim>
Problem<equationsType, dim>::Problem(Parameters<dim>& parameters, Equations<equationsType, dim>& equations,
#ifdef HAVE_MPI
  parallel::distributed::Triangulation<dim>& triangulation,
#else
  Triangulation<dim>& triangulation,
#endif
  InitialCondition<equationsType, dim>& initial_condition, BoundaryCondition<equationsType, dim>& boundary_conditions) :
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
  last_output_time(0.), time(0.),
  time_step_number(0),
  mag(dim + 2),
  update_flags(update_values | update_JxW_values | update_gradients),
  face_update_flags(update_values | update_JxW_values | update_normal_vectors | update_q_points | update_gradients),
  neighbor_face_update_flags(update_values | update_q_points),
  fe_v_cell(mapping, fe, quadrature, update_flags),
  fe_v_face(mapping, fe, face_quadrature, face_update_flags),
  fe_v_face_neighbor(mapping, fe, face_quadrature, neighbor_face_update_flags),
  fe_v_subface(mapping, fe, face_quadrature, face_update_flags),
  fe_v_subface_neighbor(mapping, fe, face_quadrature, neighbor_face_update_flags),
  adaptivity(0),
  solver(new AztecOO()),
  reset_after_refinement(true)
{
  n_quadrature_points_cell = quadrature.get_points().size();
  fluxes_old.resize(n_quadrature_points_cell);
  W_prev.resize(n_quadrature_points_cell);
  n_quadrature_points_face = face_quadrature.get_points().size();
  Wplus_old.resize(n_quadrature_points_face);
  Wgrad_plus_old.resize(n_quadrature_points_face);
  Wminus_old.resize(n_quadrature_points_face);
  normal_fluxes_old.resize(n_quadrature_points_face);

  if (parameters.num_flux_type == parameters.hlld)
    this->numFlux = new NumFluxHLLD<equationsType, dim>(this->parameters);
  else if (parameters.num_flux_type == parameters.lax_friedrich)
    this->numFlux = new NumFluxLaxFriedrich<equationsType, dim>(this->parameters);

  if (parameters.slope_limiter == parameters.vertexBased)
    this->slopeLimiter = new VertexBasedSlopeLimiter<equationsType, dim>(parameters, mapping, fe, dof_handler, dofs_per_cell, triangulation, dof_indices, component_ii, is_primitive);
  else if (parameters.slope_limiter == parameters.barthJespersen)
    this->slopeLimiter = new BarthJespersenSlopeLimiter<equationsType, dim>(parameters, mapping, fe, dof_handler, dofs_per_cell, triangulation, dof_indices, component_ii, is_primitive);
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::set_adaptivity(Adaptivity<dim>* adaptivity)
{
  this->adaptivity = adaptivity;
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
  DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
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
  for (cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
  {
    // Only assemble what belongs to this process.
    if (!cell->is_locally_owned())
      continue;

    fe_v_cell.reinit(cell);

    if (assemble_matrix)
      cell_matrix = 0;
    cell_rhs = 0;

    cell->get_dof_indices(dof_indices);

    if (parameters.debug & parameters.DetailSteps)
      LOGL(2, "Cell: " << ith_cell);
    ith_cell++;

    // Assemble the volumetric integrals.
    assemble_cell_term(cell_matrix, cell_rhs, assemble_matrix);

    // Assemble the face integrals, only after the first (projection) step
    //if (time_step_number > 0)
    {
      for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
      {
        if (parameters.debug & parameters.DetailSteps)
          LOG(3, "Face: " << face_no);
        // Boundary face - here we pass the boundary id
        if (cell->at_boundary(face_no) && !(this->parameters.is_periodic_boundary(cell->face(face_no)->boundary_id())))
        {
          if (parameters.debug & parameters.DetailSteps)
            LOGL(1, " - boundary");
          fe_v_face.reinit(cell, face_no);
          assemble_face_term(face_no, fe_v_face, fe_v_face, true, cell->face(face_no)->boundary_id(), cell_rhs);
        }
        else
        {
          // Here the neighbor face is more split than the current one (has children with respect to the current face of the current element), we need to assemble sub-face by sub-face
          // Not performed if there is no adaptivity involved.
          if (cell->neighbor_or_periodic_neighbor(face_no)->has_children())
          {
            int n_children = cell->face(face_no)->number_of_children();
            unsigned int neighbor2;
            if (this->parameters.is_periodic_boundary(cell->face(face_no)->boundary_id()))
              neighbor2 = cell->periodic_neighbor_of_periodic_neighbor(face_no);
            else
              neighbor2 = cell->neighbor_of_neighbor(face_no);

            if (parameters.debug & parameters.DetailSteps)
              LOGL(1, " - neighbor more split, " << n_children << " children");

            for (unsigned int subface_no = 0; subface_no < n_children; ++subface_no)
            {
              const typename DoFHandler<dim>::active_cell_iterator neighbor_child =
                (this->parameters.is_periodic_boundary(cell->face(face_no)->boundary_id()) ?
                  cell->periodic_neighbor_child_on_subface(face_no, subface_no) :
                  cell->neighbor_child_on_subface(face_no, subface_no));

              fe_v_subface.reinit(cell, face_no, subface_no);
              fe_v_face_neighbor.reinit(neighbor_child, neighbor2);
              neighbor_child->get_dof_indices(dof_indices_neighbor);

              assemble_face_term(face_no, fe_v_subface, fe_v_face_neighbor, false, numbers::invalid_unsigned_int, cell_rhs);
            }
          }
          // Here the neighbor face is less split than the current one, there is some transformation needed.
          // Not performed if there is no adaptivity involved.
          else if (cell->neighbor_or_periodic_neighbor(face_no)->level() != cell->level())
          {
            if (parameters.debug & parameters.DetailSteps)
              LOGL(1, " - neighbor less split");
            const typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor_or_periodic_neighbor(face_no);
            Assert(neighbor->level() == cell->level() - 1, ExcInternalError());
            neighbor->get_dof_indices(dof_indices_neighbor);

            const std::pair<unsigned int, unsigned int> faceno_subfaceno =
              (this->parameters.is_periodic_boundary(cell->face(face_no)->boundary_id()) ?
                cell->periodic_neighbor_of_coarser_periodic_neighbor(face_no) :
                cell->neighbor_of_coarser_neighbor(face_no));

            const unsigned int neighbor_face_no = faceno_subfaceno.first, neighbor_subface_no = faceno_subfaceno.second;

            fe_v_face.reinit(cell, face_no);
            fe_v_subface_neighbor.reinit(neighbor, neighbor_face_no, neighbor_subface_no);

            assemble_face_term(face_no, fe_v_face, fe_v_subface_neighbor, false, numbers::invalid_unsigned_int, cell_rhs);
          }
          // Here the neighbor face fits exactly the current face of the current element, this is the 'easy' part.
          // This is the only face assembly case performed without adaptivity.
          else
          {
            if (parameters.debug & parameters.DetailSteps)
              LOGL(1, " - neighbor equally split");
            const typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor_or_periodic_neighbor(face_no);
            neighbor->get_dof_indices(dof_indices_neighbor);

            fe_v_face.reinit(cell, face_no);
            const unsigned int neighbor2 =
              (this->parameters.is_periodic_boundary(cell->face(face_no)->boundary_id()) ?
                cell->periodic_neighbor_of_periodic_neighbor(face_no) :
                cell->neighbor_of_neighbor(face_no));

            fe_v_face_neighbor.reinit(neighbor, neighbor2);
            assemble_face_term(face_no, fe_v_face, fe_v_face_neighbor, false, numbers::invalid_unsigned_int, cell_rhs);
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
Problem<equationsType, dim>::assemble_cell_term(FullMatrix<double>& cell_matrix, Vector<double>& cell_rhs, bool assemble_matrix)
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
              val += fe_v_cell.JxW(q) * fe_v_cell.shape_value(i, q) * fe_v_cell.shape_value(j, q);
          }
          else if (!is_primitive[j] && (component_ii[i] >= 5))
          {
            for (unsigned int q = 0; q < n_quadrature_points_cell; ++q)
            {
              Tensor<1, dim> fe_v_value = fe_v_cell[mag].value(j, q);
              val += fe_v_cell.JxW(q) * fe_v_cell.shape_value(i, q) * fe_v_value[component_ii[i] - 5];
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
              val += fe_v_cell.JxW(q) * fe_v_cell.shape_value(j, q) * fe_v_value[component_ii[j] - 5];
            }
          }
          else if (!is_primitive[j])
          {
            for (unsigned int q = 0; q < n_quadrature_points_cell; ++q)
            {
              Tensor<1, dim> fe_v_value_i = fe_v_cell[mag].value(i, q);
              Tensor<1, dim> fe_v_value_j = fe_v_cell[mag].value(j, q);
              val += fe_v_cell.JxW(q) * fe_v_value_i * fe_v_value_j;
            }
          }
        }

        cell_matrix(i, j) += val;
      }
    }
  }

  if (time_step_number == 0)
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
          Tensor<1, dim> fe_v_value = fe_v_cell[mag].value(i, q);
          for (unsigned int d = 0; d < dim; d++)
            W_prev[q][5 + d] += prev_solution(dof_indices[i]) * fe_v_value[d];
        }
        else
          W_prev[q][component_ii[i]] += prev_solution(dof_indices[i]) * fe_v_cell.shape_value(i, q);
      }
    }
  }

  for (unsigned int i = 0; i < dofs_per_cell; ++i)
  {
    double val = 0.;
    for (unsigned int q = 0; q < n_quadrature_points_cell; ++q)
    {
      if (is_primitive[i])
        val += fe_v_cell.JxW(q) * W_prev[q][component_ii[i]] * fe_v_cell.shape_value(i, q);
      else
      {
        Tensor<1, dim> fe_v_value = fe_v_cell[mag].value(i, q);
        for (unsigned int d = 0; d < dim; d++)
          val += fe_v_cell.JxW(q) * W_prev[q][5 + d] * fe_v_value[d];
      }
    }
    cell_rhs(i) += val;
  }

  //if (time_step_number > 0)
  {
    for (unsigned int q = 0; q < n_quadrature_points_cell; ++q)
      equations.compute_flux_matrix(W_prev[q], fluxes_old[q], this->parameters);

    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      double val = 0.;
      for (unsigned int q = 0; q < n_quadrature_points_cell; ++q)
      {
        if (is_primitive[i])
        {
          if (!basis_fn_is_constant[i])
            for (int d = 0; d < dim; d++)
              val += fe_v_cell.JxW(q) * parameters.current_time_step_length * fluxes_old[q][component_ii[i]][d] * fe_v_cell.shape_grad(i, q)[d];
        }
        else
        {
          Tensor<2, dim> fe_v_grad = fe_v_cell[mag].gradient(i, q);
          for (unsigned int d = 0; d < dim; d++)
            for (int e = 0; e < dim; e++)
              val += fe_v_cell.JxW(q) * parameters.current_time_step_length * fluxes_old[q][5 + d][e] * fe_v_grad[d][e];
        }
      }
      if (std::isnan(val))
      {
        LOG(0, ": isnan: " << val);
        LOG(0, ": i: " << i << ", ci: " << (!is_primitive[i] ? 1 : fe_v_cell.get_fe().system_to_component_index(i).first));
        LOG(0, ": point: " << fe_v_cell.quadrature_point(0)[0] << ", " << fe_v_cell.quadrature_point(0)[1] << ", " << fe_v_cell.quadrature_point(0)[2]);
        for (int j = 0; j < Equations<equationsType, dim>::n_components; j++)
          LOGL(1, ": W [" << j << "]: " << (double)W_prev[0][j]);
        exit(1);
      }
      cell_rhs(i) += val;
    }
  }
}

template <EquationsType equationsType, int dim>
void
Problem<equationsType, dim>::assemble_face_term(const unsigned int face_no, const FEFaceValuesBase<dim> &fe_v, const FEFaceValuesBase<dim> &fe_v_neighbor,
  const bool external_face, const unsigned int boundary_id, Vector<double>& cell_rhs)
{
  // This loop is preparation - calculate all states (Wplus on the current element side of the currently assembled face, Wminus on the other side).
  if (time_step_number == 0)
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
      {
        Wplus_old[q][c] = Wminus_old[q][c] = 0.;
        for (int d = 0; d < dim; d++)
          Wgrad_plus_old[q][c][d] = 0.;
      }
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        if (fe_v.get_fe().has_support_on_face(i, face_no))
        {
          if (!is_primitive[i])
          {
            Tensor<1, dim> fe_v_value = fe_v_cell[mag].value(i, q);
            Tensor<2, dim> fe_v_grad = fe_v_cell[mag].gradient(i, q);
            for (int d = 0; d < dim; d++)
            {
              Wplus_old[q][5 + d] += prev_solution(dof_indices[i]) * fe_v_value[d];
              for (int e = 0; e < dim; e++)
                Wgrad_plus_old[q][5 + d][e] += prev_solution(dof_indices[i]) * fe_v_grad[d][e];
            }
          }
          else
          {
            Wplus_old[q][component_ii[i]] += prev_solution(dof_indices[i]) * fe_v.shape_value(i, q);
            for (int d = 0; d < dim; d++)
              Wgrad_plus_old[q][component_ii[i]][d] += prev_solution(dof_indices[i]) * fe_v.shape_grad(i, q)[d];
          }
          if (!external_face)
          {
            if (!is_primitive[i])
            {
              Tensor<1, dim> fe_v_value_neighbor = fe_v_neighbor[mag].value(i, q);
              for (int d = 0; d < dim; d++)
                Wminus_old[q][5 + d] += prev_solution(dof_indices_neighbor[i]) * fe_v_value_neighbor[d];
            }
            else
              Wminus_old[q][component_ii[i]] += prev_solution(dof_indices_neighbor[i]) * fe_v_neighbor.shape_value(i, q);
          }
        }
      }
    }
  }

  for (unsigned int q = 0; q < n_quadrature_points_face; ++q)
  {
    if (external_face)
      boundary_conditions.bc_vector_value(boundary_id, fe_v.quadrature_point(q), fe_v.normal_vector(q), Wminus_old[q], Wgrad_plus_old[q], Wplus_old[q], this->time, this->cell);
    
    // Once we have the states on both sides of the face, we need to calculate the numerical flux.
    this->numFlux->numerical_normal_flux(fe_v.normal_vector(q), Wplus_old[q], Wminus_old[q], normal_fluxes_old[q], max_signal_speed);

    // Some debugging outputs.
    if ((parameters.debug & parameters.Assembling) || (parameters.debug & parameters.NumFlux))
    {
      LOG(0, "point_i: " << q);

      LOG(1, "q: " << fe_v.quadrature_point(q) << ", n: " << fe_v.normal_vector(q)[0] << ", " << fe_v.normal_vector(q)[1] << ", " << fe_v.normal_vector(q)[2]);
      LOG(1, "Wplus: ");
      for (unsigned int i = 0; i < Equations<equationsType, dim>::n_components; i++)
        LOG(0, Wplus_old[q][i] << (i < 7 ? ", " : ""));

      LOG(1, "Wminus: ");
      for (unsigned int i = 0; i < Equations<equationsType, dim>::n_components; i++)
        LOG(0, Wminus_old[q][i] << (i < 7 ? ", " : ""));

      LOG(1, "Num F: ");
      for (unsigned int i = 0; i < Equations<equationsType, dim>::n_components; i++)
      {
        LOG(0, normal_fluxes_old[q][i] << (i < 7 ? ", " : ""));
        if ((i + 1) == Equations<equationsType, dim>::n_components)
          std::cout << std::endl;
      }
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
          val += this->parameters.current_time_step_length * (normal_fluxes_old[q][5] * fe_v_value[0] + normal_fluxes_old[q][6] * fe_v_value[1] + normal_fluxes_old[q][7] * fe_v_value[2]) * fe_v.JxW(q);
        }
        else
          val += this->parameters.current_time_step_length * normal_fluxes_old[q][component_ii[i]] * fe_v.shape_value(i, q) * fe_v.JxW(q);

        if (std::isnan(val))
        {
          numFlux->numerical_normal_flux(fe_v.normal_vector(q), Wplus_old[q], Wminus_old[q], normal_fluxes_old[q], max_signal_speed);
          LOG(0, ": isnan: " << val);
          LOG(0, ": i: " << i << ", ci: " << (!is_primitive[i] ? 1 : fe_v.get_fe().system_to_component_index(i).first));
          LOG(0, ": point: " << fe_v.quadrature_point(q)[0] << ", " << fe_v.quadrature_point(q)[1] << ", " << fe_v.quadrature_point(q)[2]);
          LOG(0, ": normal: " << fe_v.normal_vector(q)[0] << ", " << fe_v.normal_vector(q)[1] << ", " << fe_v.normal_vector(q)[2]);
          for (int j = 0; j < Equations<equationsType, dim>::n_components; j++)
            LOGL(1, ": W+ [" << j << "]: " << (double)Wplus_old[q][j] << ", W- [" << j << "]: " << (double)Wminus_old[q][j] << ", F [" << j << "]: " << (double)normal_fluxes_old[q][j]);
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
    if (this->reset_after_refinement)
    {
      delete solver;
      solver = new AztecOO();
    }

    dealii::LinearAlgebraTrilinos::MPI::Vector completely_distributed_solution(locally_owned_dofs, mpi_communicator);

    Epetra_Vector x(View, system_matrix.trilinos_matrix().DomainMap(), completely_distributed_solution.begin());
    Epetra_Vector b(View, system_matrix.trilinos_matrix().RangeMap(), system_rhs.begin());

    solver->SetAztecOption(AZ_output, (parameters.output == Parameters<dim>::quiet_solver ? AZ_none : AZ_all));
    solver->SetAztecOption(AZ_solver, AZ_gmres);
    solver->SetRHS(&b);
    solver->SetLHS(&x);

    solver->SetAztecOption(AZ_precond, AZ_dom_decomp);
    solver->SetAztecOption(AZ_subdomain_solve, AZ_ilut);
    solver->SetAztecOption(AZ_overlap, 0);
    solver->SetAztecOption(AZ_reorder, 0);
    solver->SetAztecParam(AZ_drop, parameters.ilut_drop);
    solver->SetAztecParam(AZ_ilut_fill, parameters.ilut_fill);
    solver->SetAztecParam(AZ_athresh, parameters.ilut_atol);
    solver->SetAztecParam(AZ_rthresh, parameters.ilut_rtol);

    solver->SetUserMatrix(const_cast<Epetra_CrsMatrix *> (&system_matrix.trilinos_matrix()));

    solver->Iterate(parameters.max_iterations, parameters.linear_residual);

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
void Problem<equationsType, dim>::output_results(bool use_prev_solution) const
{
  typename Equations<equationsType, dim>::Postprocessor postprocessor(parameters);

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);

  // Solution components.
  data_out.add_data_vector(use_prev_solution ? prev_solution : current_limited_solution, equations.component_names(), DataOut<dim>::type_dof_data, equations.component_interpretation());

  // Derived quantities.
  data_out.add_data_vector(use_prev_solution ? prev_solution : current_limited_solution, postprocessor);

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
  const std::string filename_base = (parameters.output_file_prefix.length() > 0 ? parameters.output_file_prefix : (use_prev_solution ? "prev_solution" : "solution")) + "-" + Utilities::int_to_string(output_file_number, 3);

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
  std::string filename = (parameters.output_file_prefix.length() > 0 ? parameters.output_file_prefix : (use_prev_solution ? "prev_solution" : "solution")) + "-" + Utilities::int_to_string(output_file_number, 3) + ".vtk";
  std::ofstream output(filename.c_str());
  data_out.write_vtk(output);
#endif

  ++output_file_number;
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::output_matrix(TrilinosWrappers::SparseMatrix& mat, const char* suffix) const
{
  std::ofstream m;
  std::stringstream ssm;
  ssm << time_step_number << "-" << Utilities::MPI::this_mpi_process(mpi_communicator) << "." << suffix;
  m.open(ssm.str());
  mat.print(m);
  m.close();
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::output_vector(TrilinosWrappers::MPI::Vector& vec, const char* suffix) const
{
  std::ofstream n;
  std::stringstream ssn;
  ssn << time_step_number << "-" << Utilities::MPI::this_mpi_process(mpi_communicator) << "." << suffix;
  n.open(ssn.str());
  vec.print(n, 10, false, false);
  n.close();
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::run()
{
  // Preparations.
  setup_system();
  current_limited_solution.reinit(locally_owned_dofs, mpi_communicator);
  current_unlimited_solution.reinit(locally_relevant_dofs, mpi_communicator);
  prev_solution.reinit(locally_relevant_dofs, mpi_communicator);

#ifdef OUTPUT_BASE
  output_base();
  exit(1);
#endif

  int adaptivity_step = 0;
  while (time < parameters.final_time)
  {
    // Some output.
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      LOGL(0, "Step: " << time_step_number << ", T: " << time << ", adaptivity step: " << (this->adaptivity ? adaptivity_step++ : 0));
      LOGL(0, "- number of active cells:       " << triangulation.n_global_active_cells() << std::endl << " Number of degrees of freedom: " << dof_handler.n_dofs());
    }

    // Assemble
    if(this->parameters.debug | this->parameters.BasicSteps)
      LOGL(1, "Assembling...")
    system_rhs = 0;
    if (reset_after_refinement)
      system_matrix = 0;
    assemble_system(this->reset_after_refinement);

    // Output matrix & rhs if required.
    if (parameters.output_matrix)
      output_matrix(system_matrix, "matrix");
    if (parameters.output_rhs)
      output_vector(system_rhs, "rhs");

    // Solve
    if (this->parameters.debug | this->parameters.BasicSteps)
      LOGL(1, "Solving...")
    solve();

    // Postprocess if required
    if ((this->time >= this->parameters.start_limiting_at) && parameters.limit && parameters.polynomial_order_dg > 0)
    {
      if (this->parameters.debug | this->parameters.BasicSteps)
        LOGL(1, "Postprocessing...")
        postprocess();
    }
    else
      current_limited_solution = current_unlimited_solution;

    move_time_step_handle_outputs();
  }
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::perform_reset_after_refinement()
{
  this->reset_after_refinement = true;
  this->slopeLimiter->flush_cache();
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::move_time_step_handle_outputs()
{
  if (parameters.output_solution)
    output_vector(current_limited_solution, "solution");

  if ((parameters.output_step < 0) || (time - last_output_time >= parameters.output_step))
  {
    output_results();
    last_output_time = time;
  }

  if (time_step_number > 0)
  {
    calculate_cfl_condition();
    double global_cfl_time_step = Utilities::MPI::min(this->cfl_time_step, mpi_communicator);
    parameters.current_time_step_length = global_cfl_time_step;
  }

  if (this->adaptivity)
  {
    // refine mesh
    // we use the unlimited solution here for two reasons:
    // - it has ghost elements
    // - it is a useful indicator where to refine
    LOGL(1, "Refining...")
    bool refined = this->adaptivity->refine_mesh(time_step_number, time, current_unlimited_solution, dof_handler, triangulation, mapping);
    if (refined)
    {
      // transfer solution
#ifdef HAVE_MPI
      parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector> soltrans(dof_handler);
#else
      SolutionTransfer<dim, TrilinosWrappers::MPI::Vector> soltrans(dof_handler);
#endif

      if (time_step_number > 0)
        soltrans.prepare_for_coarsening_and_refinement(prev_solution);

      // Refine the current triangulation.
      triangulation.execute_coarsening_and_refinement();

      // reinit structures, periodicity, etc.
      this->setup_system();

      current_limited_solution.reinit(locally_owned_dofs, mpi_communicator);
      current_unlimited_solution.reinit(locally_relevant_dofs, mpi_communicator);

      // Now interpolate the solution
      if (time_step_number > 0)
      {
        TrilinosWrappers::MPI::Vector interpolated_solution;
        interpolated_solution.reinit(locally_owned_dofs, mpi_communicator);
#ifdef HAVE_MPI
        soltrans.interpolate(interpolated_solution);
#else
        soltrans.interpolate(prev_solution, interpolated_solution);
#endif
        prev_solution.reinit(locally_relevant_dofs, mpi_communicator);
        this->prev_solution = interpolated_solution;
      }
      else
        prev_solution.reinit(locally_relevant_dofs, mpi_communicator);

      this->perform_reset_after_refinement();
  }
    else
    {
      this->reset_after_refinement = false;
      this->prev_solution = this->current_limited_solution;
      ++time_step_number;
      time += parameters.current_time_step_length;
    }
}
  else
  {
    this->reset_after_refinement = false;
    this->prev_solution = this->current_limited_solution;
    ++time_step_number;
    time += parameters.current_time_step_length;
  }
}

template class Problem<EquationsTypeMhd, 3>;

#include "problem.h"

template <EquationsType equationsType, int dim>
Problem<equationsType, dim>::Problem(Parameters<dim>& parameters, Equations<equationsType, dim>& equations, Triangulation<dim>& triangulation, InitialCondition<equationsType, dim>& initial_condition, BoundaryConditions<equationsType, dim>& boundary_conditions) :
  parameters(parameters),
  equations(equations),
  triangulation(triangulation),
  initial_condition(initial_condition),
  boundary_conditions(boundary_conditions),
  mapping(),
  fe(FE_DGQ<dim>(parameters.polynomial_order), Equations<equationsType, dim>::n_components),
  dof_handler(triangulation),
  quadrature(2 * parameters.polynomial_order + 1),
  face_quadrature(2 * parameters.polynomial_order + 1),
  verbose_cout(std::cout, false)
{
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::setup_system()
{
  DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
  DoFTools::make_flux_sparsity_pattern(dof_handler, dsp);

  system_matrix.reinit(dsp);
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::assemble_system()
{
  const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
  std::vector<types::global_dof_index> dof_indices_neighbor(dofs_per_cell);

  const UpdateFlags update_flags = update_values | update_gradients | update_q_points | update_JxW_values;
  const UpdateFlags face_update_flags = update_values | update_q_points | update_JxW_values | update_normal_vectors;
  const UpdateFlags neighbor_face_update_flags = update_values;

  FEValues<dim> fe_v(mapping, fe, quadrature, update_flags);
  FEFaceValues<dim> fe_v_face(mapping, fe, face_quadrature, face_update_flags);
  FESubfaceValues<dim> fe_v_subface(mapping, fe, face_quadrature, face_update_flags);
  FEFaceValues<dim> fe_v_face_neighbor(mapping, fe, face_quadrature, neighbor_face_update_flags);
  FESubfaceValues<dim> fe_v_subface_neighbor(mapping, fe, face_quadrature, neighbor_face_update_flags);

  for (typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
  {
    fe_v.reinit(cell);
    cell->get_dof_indices(dof_indices);

    assemble_cell_term(fe_v, dof_indices);

    for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
    {
      if (cell->at_boundary(face_no))
      {
        fe_v_face.reinit(cell, face_no);
        assemble_face_term(face_no, fe_v_face, fe_v_face, dof_indices, std::vector<types::global_dof_index>(), true,
          cell->face(face_no)->boundary_id(), cell->face(face_no)->diameter());
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

            assemble_face_term(face_no, fe_v_subface, fe_v_face_neighbor, dof_indices, dof_indices_neighbor, false,
              numbers::invalid_unsigned_int, neighbor_child->face(neighbor2)->diameter());
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

          assemble_face_term(face_no, fe_v_face, fe_v_face_neighbor, dof_indices, dof_indices_neighbor, false,
            numbers::invalid_unsigned_int, cell->face(face_no)->diameter());
        }
        else
        {
          const typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face_no);
          neighbor->get_dof_indices(dof_indices_neighbor);

          fe_v_face.reinit(cell, face_no);
          fe_v_face_neighbor.reinit(neighbor, cell->neighbor_of_neighbor(face_no));

          assemble_face_term(face_no, fe_v_face, fe_v_face_neighbor, dof_indices, dof_indices_neighbor, false,
            numbers::invalid_unsigned_int, cell->face(face_no)->diameter());
        }
      }
    }
  }
}

template <EquationsType equationsType, int dim>
void
Problem<equationsType, dim>::assemble_cell_term(const FEValues<dim> &fe_v, const std::vector<types::global_dof_index> &dof_indices)
{
  const unsigned int dofs_per_cell = fe_v.dofs_per_cell;
  const unsigned int n_q_points = fe_v.n_quadrature_points;

  Table<2, Sacado::Fad::DFad<double> > W(n_q_points, Equations<equationsType, dim>::n_components);

  Table<2, double> W_old(n_q_points, Equations<equationsType, dim>::n_components);

  Table<3, Sacado::Fad::DFad<double> > grad_W(n_q_points, Equations<equationsType, dim>::n_components, dim);

  Table<3, double> grad_W_old(n_q_points, Equations<equationsType, dim>::n_components, dim);

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
      W_old[q][c] = 0;
      for (unsigned int d = 0; d < dim; ++d)
      {
        grad_W[q][c][d] = 0;
        grad_W_old[q][c][d] = 0;
      }
    }
  }

  for (unsigned int q = 0; q < n_q_points; ++q)
  {
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      const unsigned int c = fe_v.get_fe().system_to_component_index(i).first;

      W[q][c] += independent_local_dof_values[i] * fe_v.shape_value_component(i, q, c);
      W_old[q][c] += old_solution(dof_indices[i]) * fe_v.shape_value_component(i, q, c);

      for (unsigned int d = 0; d < dim; d++)
      {
        grad_W[q][c][d] += independent_local_dof_values[i] * fe_v.shape_grad_component(i, q, c)[d];
        grad_W_old[q][c][d] += old_solution(dof_indices[i]) * fe_v.shape_grad_component(i, q, c)[d];
      }
    }
  }

  std::vector < std_cxx11::array <std_cxx11::array <Sacado::Fad::DFad<double>, dim>, Equations<equationsType, dim>::n_components > > flux(n_q_points);
  std::vector < std_cxx11::array< Sacado::Fad::DFad<double>, Equations<equationsType, dim>::n_components> > forcing(n_q_points);

  std::vector < std_cxx11::array <std_cxx11::array <double, dim>, Equations<equationsType, dim>::n_components > > flux_old(n_q_points);
  std::vector < std_cxx11::array< double, Equations<equationsType, dim>::n_components> > forcing_old(n_q_points);

  // LK: This is for e.g. stabilization for Euler
  std::vector < std_cxx11::array <std_cxx11::array <Sacado::Fad::DFad<double>, dim>, Equations<equationsType, dim>::n_components > > jacobian_addition(n_q_points);
  std::vector < std_cxx11::array <std_cxx11::array <double, dim>, Equations<equationsType, dim>::n_components > > jacobian_addition_old(n_q_points);

  for (unsigned int q = 0; q < n_q_points; ++q)
  {
    equations.compute_flux_matrix(W[q], flux[q]);
    equations.compute_forcing_vector(W[q], forcing[q]);
    equations.compute_jacobian_addition(fe_v.get_cell()->diameter(), grad_W[q], jacobian_addition[q]);

    equations.compute_flux_matrix(W_old[q], flux_old[q]);
    equations.compute_forcing_vector(W_old[q], forcing_old[q]);
    equations.compute_jacobian_addition(fe_v.get_cell()->diameter(), grad_W_old[q], jacobian_addition_old[q]);
  }

  for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
  {
    Sacado::Fad::DFad<double> R_i = 0;

    const unsigned int component_i = fe_v.get_fe().system_to_component_index(i).first;

    for (unsigned int point = 0; point < fe_v.n_quadrature_points; ++point)
    {
      if (parameters.is_stationary == false)
        R_i += 1.0 / parameters.time_step * (W[point][component_i] - W_old[point][component_i]) * fe_v.shape_value_component(i, point, component_i) * fe_v.JxW(point);

      for (unsigned int d = 0; d < dim; d++)
        R_i -= (parameters.theta * flux[point][component_i][d] + (1.0 - parameters.theta) * flux_old[point][component_i][d])
        * fe_v.shape_grad_component(i, point, component_i)[d] * fe_v.JxW(point);

      for (unsigned int d = 0; d < dim; d++)
        R_i += (parameters.theta * jacobian_addition[point][component_i][d] + (1.0 - parameters.theta) * jacobian_addition_old[point][component_i][d])
        * fe_v.shape_grad_component(i, point, component_i)[d] * fe_v.JxW(point);

      R_i -= (parameters.theta  * forcing[point][component_i] + (1.0 - parameters.theta) * forcing_old[point][component_i])
        * fe_v.shape_value_component(i, point, component_i) * fe_v.JxW(point);
    }

    for (unsigned int k = 0; k < dofs_per_cell; ++k)
      residual_derivatives[k] = R_i.fastAccessDx(k);

    system_matrix.add(dof_indices[i], dof_indices, residual_derivatives);
    right_hand_side(dof_indices[i]) -= R_i.val();
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
  const double                 face_diameter)
{
  const unsigned int n_q_points = fe_v.n_quadrature_points;
  const unsigned int dofs_per_cell = fe_v.dofs_per_cell;

  std::vector<Sacado::Fad::DFad<double> > independent_local_dof_values(dofs_per_cell), independent_neighbor_dof_values(external_face == false ? dofs_per_cell : 0);

  const unsigned int n_independent_variables = (external_face == false ? 2 * dofs_per_cell : dofs_per_cell);

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

  Table<2, Sacado::Fad::DFad<double> > Wplus(n_q_points, Equations<equationsType, dim>::n_components), Wminus(n_q_points, Equations<equationsType, dim>::n_components);
  Table<2, double> Wplus_old(n_q_points, Equations<equationsType, dim>::n_components), Wminus_old(n_q_points, Equations<equationsType, dim>::n_components);

  for (unsigned int q = 0; q < n_q_points; ++q)
  {
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      const unsigned int component_i = fe_v.get_fe().system_to_component_index(i).first;
      Wplus[q][component_i] += independent_local_dof_values[i] * fe_v.shape_value_component(i, q, component_i);
      Wplus_old[q][component_i] += old_solution(dof_indices[i]) * fe_v.shape_value_component(i, q, component_i);
    }
  }

  if (external_face == false)
  {
    for (unsigned int q = 0; q < n_q_points; ++q)
    {
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        const unsigned int component_i = fe_v_neighbor.get_fe().system_to_component_index(i).first;
        Wminus[q][component_i] += independent_neighbor_dof_values[i] * fe_v_neighbor.shape_value_component(i, q, component_i);
        Wminus_old[q][component_i] += old_solution(dof_indices_neighbor[i])* fe_v_neighbor.shape_value_component(i, q, component_i);
      }
    }
  }
  else
  {
    std::vector<Vector<double> > boundary_values(n_q_points, Vector<double>(Equations<equationsType, dim>::n_components));
    boundary_conditions.bc_vector_value(boundary_id, fe_v.get_quadrature_points(), boundary_values);

    for (unsigned int q = 0; q < n_q_points; q++)
    {
      equations.compute_Wminus(boundary_conditions.kind[boundary_id], fe_v.normal_vector(q), Wplus[q], boundary_values[q], Wminus[q]);
      equations.compute_Wminus(boundary_conditions.kind[boundary_id], fe_v.normal_vector(q), Wplus_old[q], boundary_values[q], Wminus_old[q]);
    }
  }

  std::vector< std_cxx11::array < Sacado::Fad::DFad<double>, Equations<equationsType, dim>::n_components> > normal_fluxes(n_q_points);
  std::vector< std_cxx11::array < double, Equations<equationsType, dim>::n_components> > normal_fluxes_old(n_q_points);

  for (unsigned int q = 0; q < n_q_points; ++q)
  {
    equations.numerical_normal_flux(fe_v.normal_vector(q), Wplus[q], Wminus[q], normal_fluxes[q]);
    equations.numerical_normal_flux(fe_v.normal_vector(q), Wplus_old[q], Wminus_old[q], normal_fluxes_old[q]);
  }

  std::vector<double> residual_derivatives(dofs_per_cell);

  for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
  {
    if (fe_v.get_fe().has_support_on_face(i, face_no) == true)
    {
      Sacado::Fad::DFad<double> R_i = 0;

      for (unsigned int point = 0; point < n_q_points; ++point)
      {
        const unsigned int component_i = fe_v.get_fe().system_to_component_index(i).first;

        R_i += (parameters.theta * normal_fluxes[point][component_i] + (1.0 - parameters.theta) * normal_fluxes_old[point][component_i])
          * fe_v.shape_value_component(i, point, component_i) * fe_v.JxW(point);
      }

      for (unsigned int k = 0; k < dofs_per_cell; ++k)
        residual_derivatives[k] = R_i.fastAccessDx(k);
      system_matrix.add(dof_indices[i], dof_indices, residual_derivatives);

      if (external_face == false)
      {
        for (unsigned int k = 0; k < dofs_per_cell; ++k)
          residual_derivatives[k] = R_i.fastAccessDx(dofs_per_cell + k);
        system_matrix.add(dof_indices[i], dof_indices_neighbor, residual_derivatives);
      }

      right_hand_side(dof_indices[i]) -= R_i.val();
    }
  }
}

template <EquationsType equationsType, int dim>
std::pair<unsigned int, double>
Problem<equationsType, dim>::solve(Vector<double> &newton_update)
{
  if (parameters.solver == Parameters<dim>::direct)
  {
    SolverControl solver_control(1, 0);
    TrilinosWrappers::SolverDirect::AdditionalData data(parameters.output == Parameters<dim>::verbose_solver);
    TrilinosWrappers::SolverDirect direct(solver_control, data);

    direct.solve(system_matrix, newton_update, right_hand_side);

    return std::pair<unsigned int, double>(solver_control.last_step(),
      solver_control.last_value());
  }

  if (parameters.solver == Parameters<dim>::gmres)
  {
    Epetra_Vector x(View, system_matrix.trilinos_matrix().DomainMap(),
      newton_update.begin());
    Epetra_Vector b(View, system_matrix.trilinos_matrix().RangeMap(),
      right_hand_side.begin());

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

    solver.SetUserMatrix(const_cast<Epetra_CrsMatrix *>
      (&system_matrix.trilinos_matrix()));

    solver.Iterate(parameters.max_iterations, parameters.linear_residual);

    return std::pair<unsigned int, double>(solver.NumIters(), solver.TrueResidual());
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

  data_out.build_patches();

  static unsigned int output_file_number = 0;
  std::string filename = "solution-" + Utilities::int_to_string(output_file_number, 3) + ".vtk";
  std::ofstream output(filename.c_str());
  data_out.write_vtk(output);

  ++output_file_number;
}

template <EquationsType equationsType, int dim>
void Problem<equationsType, dim>::run()
{
  dof_handler.clear();
  dof_handler.distribute_dofs(fe);

  // Size all of the fields.
  old_solution.reinit(dof_handler.n_dofs());
  current_solution.reinit(dof_handler.n_dofs());
  newton_initial_guess.reinit(dof_handler.n_dofs());
  right_hand_side.reinit(dof_handler.n_dofs());

  setup_system();

  VectorTools::interpolate(dof_handler, initial_condition, old_solution);
  current_solution = old_solution;
  newton_initial_guess = old_solution;

  output_results();

  Vector<double> newton_update(dof_handler.n_dofs());

  double time = 0;
  double next_output = time + parameters.output_step;

  newton_initial_guess = old_solution;
  while (time < parameters.final_time)
  {
    if (time > this->parameters.initialization_time) {
      this->parameters.time_step = this->parameters.time_step_after_initialization;
      this->parameters.theta = this->parameters.theta_after_initialization;
    }

    std::cout << "T=" << time << std::endl << "   Number of active cells:       " << triangulation.n_active_cells() << std::endl
      << "   Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl << std::endl;
    std::cout << "   NonLin Res     Lin Iter       Lin Res" << std::endl << "   _____________________________________" << std::endl;

    unsigned int nonlin_iter = 0;
    current_solution = newton_initial_guess;
    while (true)
    {
      system_matrix = 0;

      right_hand_side = 0;
      assemble_system();

      const double res_norm = right_hand_side.l2_norm();
      if (std::fabs(res_norm) < parameters.nonlinear_residual_norm_threshold)
      {
        std::printf("   %-16.3e (converged)\n\n", res_norm);
        break;
      }
      else
      {
        newton_update = 0;
        std::pair<unsigned int, double> convergence = solve(newton_update);
        current_solution += newton_update;
        std::printf("   %-16.3e %04d        %-5.2e\n", res_norm, convergence.first, convergence.second);
      }

      ++nonlin_iter;
      AssertThrow(nonlin_iter <= parameters.max_nonlinear_iterations, ExcMessage("No convergence in nonlinear solver"));
    }

    time += parameters.time_step;

    if (parameters.output_step < 0)
      output_results();
    else if (time >= next_output)
    {
      output_results();
      next_output += parameters.output_step;
    }

    newton_initial_guess = current_solution;

    // Tohle je extrapolace - udelam stejny krok do dalsiho reseni jaky jsem udelal do tohohle.
    newton_initial_guess.sadd(2.0, -1.0, old_solution);

    old_solution = current_solution;
  }
}

template class Problem<EquationsTypeEuler, 2>;
template class Problem<EquationsTypeEuler, 3>;

template class Problem<EquationsTypeMhd, 3>;
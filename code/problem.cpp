#include "problem.h"

template <int dim>
Problem<dim>::Problem(Parameters<dim>& parameters) : parameters(parameters), mapping(), fe(FE_DGQ<dim>(parameters.polynomial_order), Equations<dim>::n_components),
dof_handler(triangulation), quadrature(2 * parameters.polynomial_order + 1), face_quadrature(2 * parameters.polynomial_order + 1), verbose_cout(std::cout, false)
{
}

template <int dim>
void Problem<dim>::setup_system()
{
  DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
  DoFTools::make_flux_sparsity_pattern(dof_handler, dsp);

  system_matrix.reinit(dsp);
}

template <int dim>
void Problem<dim>::assemble_system()
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

        // The other possibility we have to care for is if the neighbor
        // is coarser than the current cell (in particular, because of
        // the usual restriction of only one hanging node per face, the
        // neighbor must be exactly one level coarser than the current
        // cell, something that we check with an assertion). Again, we
        // then integrate over this interface:
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

template <int dim>
void
Problem<dim>::assemble_cell_term(const FEValues<dim> &fe_v, const std::vector<types::global_dof_index> &dof_indices)
{
  const unsigned int dofs_per_cell = fe_v.dofs_per_cell;
  const unsigned int n_q_points = fe_v.n_quadrature_points;

  Table<2, Sacado::Fad::DFad<double> > W(n_q_points, Equations<dim>::n_components);

  Table<2, double> W_old(n_q_points, Equations<dim>::n_components);

  Table<3, Sacado::Fad::DFad<double> > grad_W(n_q_points, Equations<dim>::n_components, dim);

  Table<3, double> grad_W_old(n_q_points, Equations<dim>::n_components, dim);

  std::vector<double> residual_derivatives(dofs_per_cell);

  // Next, we have to define the independent variables that we will try to
  // determine by solving a Newton step. These independent variables are the
  // values of the local degrees of freedom which we extract here:
  std::vector<Sacado::Fad::DFad<double> > independent_local_dof_values(dofs_per_cell);
  for (unsigned int i = 0; i < dofs_per_cell; ++i)
    independent_local_dof_values[i] = current_solution(dof_indices[i]);

  // The next step incorporates all the magic: we declare a subset of the
  // autodifferentiation variables as independent degrees of freedom,
  // whereas all the other ones remain dependent functions. These are
  // precisely the local degrees of freedom just extracted. All calculations
  // that reference them (either directly or indirectly) will accumulate
  // sensitivities with respect to these variables.
  //
  // In order to mark the variables as independent, the following does the
  // trick, marking <code>independent_local_dof_values[i]</code> as the
  // $i$th independent variable out of a total of
  // <code>dofs_per_cell</code>:
  for (unsigned int i = 0; i < dofs_per_cell; ++i)
    independent_local_dof_values[i].diff(i, dofs_per_cell);

  // After all these declarations, let us actually compute something. First,
  // the values of <code>W</code>, <code>W_old</code>, <code>grad_W</code>
  // and <code>grad_W_old</code>, which we can compute from the local DoF values
  // by using the formula $W(x_q)=\sum_i \mathbf W_i \Phi_i(x_q)$, where
  // $\mathbf W_i$ is the $i$th entry of the (local part of the) solution
  // vector, and $\Phi_i(x_q)$ the value of the $i$th vector-valued shape
  // function evaluated at quadrature point $x_q$. The gradient can be
  // computed in a similar way.
  //
  // Ideally, we could compute this information using a call into something
  // like FEValues::get_function_values and FEValues::get_function_gradients,
  // but since (i) we would have to extend the FEValues class for this, and
  // (ii) we don't want to make the entire <code>old_solution</code> vector
  // fad types, only the local cell variables, we explicitly code the loop
  // above. Before this, we add another loop that initializes all the fad
  // variables to zero:
  for (unsigned int q = 0; q < n_q_points; ++q)
    for (unsigned int c = 0; c < Equations<dim>::n_components; ++c)
    {
      W[q][c] = 0;
      W_old[q][c] = 0;
      for (unsigned int d = 0; d < dim; ++d)
      {
        grad_W[q][c][d] = 0;
        grad_W_old[q][c][d] = 0;
      }
    }

  for (unsigned int q = 0; q < n_q_points; ++q)
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

  // Next, in order to compute the cell contributions, we need to evaluate
  // $\mathbf{F}({\mathbf w}^k_{n+1})$, $\mathbf{G}({\mathbf w}^k_{n+1})$ and
  // $\mathbf{F}({\mathbf w}_n)$, $\mathbf{G}({\mathbf w}_n)$ at all quadrature
  // points. To store these, we also need to allocate a bit of memory. Note
  // that we compute the flux matrices and right hand sides in terms of
  // autodifferentiation variables, so that the Jacobian contributions can
  // later easily be computed from it:

  std::vector < std_cxx11::array <std_cxx11::array <Sacado::Fad::DFad<double>, dim>, Equations<dim>::n_components > > flux(n_q_points);

  std::vector < std_cxx11::array <std_cxx11::array <double, dim>, Equations<dim>::n_components > > flux_old(n_q_points);

  std::vector < std_cxx11::array< Sacado::Fad::DFad<double>, Equations<dim>::n_components> > forcing(n_q_points);

  std::vector < std_cxx11::array< double, Equations<dim>::n_components> > forcing_old(n_q_points);

  for (unsigned int q = 0; q < n_q_points; ++q)
  {
    Equations<dim>::compute_flux_matrix(W_old[q], flux_old[q]);
    Equations<dim>::compute_forcing_vector(W_old[q], forcing_old[q]);
    Equations<dim>::compute_flux_matrix(W[q], flux[q]);
    Equations<dim>::compute_forcing_vector(W[q], forcing[q]);
  }

  // We now have all of the pieces in place, so perform the assembly.  We
  // have an outer loop through the components of the system, and an inner
  // loop over the quadrature points, where we accumulate contributions to
  // the $i$th residual $R_i$. The general formula for this residual is
  // given in the introduction and at the top of this function. We can,
  // however, simplify it a bit taking into account that the $i$th
  // (vector-valued) test function $\mathbf{z}_i$ has in reality only a
  // single nonzero component (more on this topic can be found in the @ref
  // vector_valued module). It will be represented by the variable
  // <code>component_i</code> below. With this, the residual term can be
  // re-written as
  // @f{eqnarray*}
  // R_i &=&
  // \left(\frac{(\mathbf{w}_{n+1} -
  // \mathbf{w}_n)_{\text{component\_i}}}{\delta
  // t},(\mathbf{z}_i)_{\text{component\_i}}\right)_K
  // \\ &-& \sum_{d=1}^{\text{dim}} \left(  \theta \mathbf{F}
  // ({\mathbf{w}^k_{n+1}})_{\text{component\_i},d} + (1-\theta)
  // \mathbf{F} ({\mathbf{w}_{n}})_{\text{component\_i},d}  ,
  // \frac{\partial(\mathbf{z}_i)_{\text{component\_i}}} {\partial
  // x_d}\right)_K
  // \\ &+& \sum_{d=1}^{\text{dim}} h^{\eta} \left( \theta \frac{\partial
  // (\mathbf{w}^k_{n+1})_{\text{component\_i}}}{\partial x_d} + (1-\theta)
  // \frac{\partial (\mathbf{w}_n)_{\text{component\_i}}}{\partial x_d} ,
  // \frac{\partial (\mathbf{z}_i)_{\text{component\_i}}}{\partial x_d} \right)_K
  // \\ &-& \left( \theta\mathbf{G}({\mathbf{w}^k_n+1} )_{\text{component\_i}} +
  // (1-\theta)\mathbf{G}({\mathbf{w}_n})_{\text{component\_i}} ,
  // (\mathbf{z}_i)_{\text{component\_i}} \right)_K ,
  // @f}
  // where integrals are
  // understood to be evaluated through summation over quadrature points.
  //
  // We initially sum all contributions of the residual in the positive
  // sense, so that we don't need to negative the Jacobian entries.  Then,
  // when we sum into the <code>right_hand_side</code> vector, we negate
  // this residual.
  for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
  {
    Sacado::Fad::DFad<double> R_i = 0;

    const unsigned int component_i = fe_v.get_fe().system_to_component_index(i).first;

    // The residual for each row (i) will be accumulating into this fad
    // variable.  At the end of the assembly for this row, we will query
    // for the sensitivities to this variable and add them into the
    // Jacobian.
    for (unsigned int point = 0; point < fe_v.n_quadrature_points; ++point)
    {
      if (parameters.is_stationary == false)
        R_i += 1.0 / parameters.time_step * (W[point][component_i] - W_old[point][component_i]) * fe_v.shape_value_component(i, point, component_i) * fe_v.JxW(point);

      for (unsigned int d = 0; d < dim; d++)
        R_i -= (parameters.theta * flux[point][component_i][d] + (1.0 - parameters.theta) * flux_old[point][component_i][d])
        * fe_v.shape_grad_component(i, point, component_i)[d] * fe_v.JxW(point);

      for (unsigned int d = 0; d<dim; d++)
        R_i += std::pow(fe_v.get_cell()->diameter(), 2.0)
        * (parameters.theta * grad_W[point][component_i][d] + (1.0 - parameters.theta) * grad_W_old[point][component_i][d])
        * fe_v.shape_grad_component(i, point, component_i)[d] * fe_v.JxW(point);

      R_i -= (parameters.theta  * forcing[point][component_i] + (1.0 - parameters.theta) * forcing_old[point][component_i])
        * fe_v.shape_value_component(i, point, component_i) * fe_v.JxW(point);
    }

    // At the end of the loop, we have to add the sensitivities to the
    // matrix and subtract the residual from the right hand side. Trilinos
    // FAD data type gives us access to the derivatives using
    // <code>R_i.fastAccessDx(k)</code>, so we store the data in a
    // temporary array. This information about the whole row of local dofs
    // is then added to the Trilinos matrix at once (which supports the
    // data types we have chosen).
    for (unsigned int k = 0; k < dofs_per_cell; ++k)
      residual_derivatives[k] = R_i.fastAccessDx(k);

    system_matrix.add(dof_indices[i], dof_indices, residual_derivatives);
    right_hand_side(dof_indices[i]) -= R_i.val();
  }
}

template <int dim>
void
Problem<dim>::assemble_face_term(const unsigned int           face_no,
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

  Table<2, Sacado::Fad::DFad<double> > Wplus(n_q_points, Equations<dim>::n_components), Wminus(n_q_points, Equations<dim>::n_components);
  Table<2, double> Wplus_old(n_q_points, Equations<dim>::n_components), Wminus_old(n_q_points, Equations<dim>::n_components);

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
  // On the other hand, if this is an external boundary face, then the
  // values of $\mathbf{W}^-$ will be either functions of $\mathbf{W}^+$, or they will be
  // prescribed, depending on the kind of boundary condition imposed here.
  //
  // To start the evaluation, let us ensure that the boundary id specified
  // for this boundary is one for which we actually have data in the
  // parameters object. Next, we evaluate the function object for the
  // inhomogeneity.  This is a bit tricky: a given boundary might have both
  // prescribed and implicit values.  If a particular component is not
  // prescribed, the values evaluate to zero and are ignored below.
  //
  // The rest is done by a function that actually knows the specifics of
  // Euler equation boundary conditions. Note that since we are using fad
  // variables here, sensitivities will be updated appropriately, a process
  // that would otherwise be tremendously complicated.
  else
  {
    Assert(boundary_id < Parameters<dim>::max_n_boundaries, ExcIndexRange(boundary_id, 0, Parameters<dim>::max_n_boundaries));

    std::vector<Vector<double> > boundary_values(n_q_points, Vector<double>(Equations<dim>::n_components));
    Parameters<dim>::bc_vector_value(boundary_id, fe_v.get_quadrature_points(), boundary_values);

    for (unsigned int q = 0; q < n_q_points; q++)
    {
      Equations<dim>::compute_Wminus(parameters.boundary_conditions[boundary_id].kind, fe_v.normal_vector(q), Wplus[q], boundary_values[q], Wminus[q]);
      Equations<dim>::compute_Wminus(parameters.boundary_conditions[boundary_id].kind, fe_v.normal_vector(q), Wplus_old[q], boundary_values[q], Wminus_old[q]);
    }
  }

  std::vector< std_cxx11::array < Sacado::Fad::DFad<double>, Equations<dim>::n_components> > normal_fluxes(n_q_points);
  std::vector< std_cxx11::array < double, Equations<dim>::n_components> > normal_fluxes_old(n_q_points);

  double alpha;

  switch (parameters.stabilization_kind)
  {
  case Parameters<dim>::constant_stabilization:
    alpha = parameters.stabilization_value;
    break;
  case Parameters<dim>::mesh_dependent_stabilization:
    alpha = face_diameter / (2.0*parameters.time_step);
    break;
  default:
    Assert(false, ExcNotImplemented());
    alpha = 1;
  }

  for (unsigned int q = 0; q < n_q_points; ++q)
  {
    Equations<dim>::numerical_normal_flux(fe_v.normal_vector(q), Wplus[q], Wminus[q], alpha, normal_fluxes[q]);
    Equations<dim>::numerical_normal_flux(fe_v.normal_vector(q), Wplus_old[q], Wminus_old[q], alpha, normal_fluxes_old[q]);
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

template <int dim>
std::pair<unsigned int, double>
Problem<dim>::solve(Vector<double> &newton_update)
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

  // Likewise, if we are to use an iterative solver, we use Aztec's GMRES
  // solver. We could use the Trilinos wrapper classes for iterative
  // solvers and preconditioners here as well, but we choose to use an
  // Aztec solver directly. For the given problem, Aztec's internal
  // preconditioner implementations are superior over the ones deal.II has
  // wrapper classes to, so we use ILU-T preconditioning within the
  // AztecOO solver and set a bunch of options that can be changed from
  // the parameter file.
  //
  // There are two more practicalities: Since we have built our right hand
  // side and solution vector as deal.II Vector objects (as opposed to the
  // matrix, which is a Trilinos object), we must hand the solvers
  // Trilinos Epetra vectors.  Luckily, they support the concept of a
  // 'view', so we just send in a pointer to our deal.II vectors. We have
  // to provide an Epetra_Map for the vector that sets the parallel
  // distribution, which is just a dummy object in serial. The easiest way
  // is to ask the matrix for its map, and we're going to be ready for
  // matrix-vector products with it.
  //
  // Secondly, the Aztec solver wants us to pass a Trilinos
  // Epetra_CrsMatrix in, not the deal.II wrapper class itself. So we
  // access to the actual Trilinos matrix in the Trilinos wrapper class by
  // the command trilinos_matrix(). Trilinos wants the matrix to be
  // non-constant, so we have to manually remove the constantness using a
  // const_cast.
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

    return std::pair<unsigned int, double>(solver.NumIters(),
      solver.TrueResidual());
  }
}

template <int dim>
void Problem<dim>::output_results() const
{
  typename Equations<dim>::Postprocessor postprocessor(parameters.schlieren_plot);

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);

  data_out.add_data_vector(current_solution, Equations<dim>::component_names(), DataOut<dim>::type_dof_data, Equations<dim>::component_interpretation());
  data_out.add_data_vector(current_solution, postprocessor);

  data_out.build_patches();

  static unsigned int output_file_number = 0;
  std::string filename = "solution-" + Utilities::int_to_string(output_file_number, 3) + ".vtk";
  std::ofstream output(filename.c_str());
  data_out.write_vtk(output);

  ++output_file_number;
}

template <int dim>
void Problem<dim>::run()
{
  {
    GridIn<dim> grid_in;
    grid_in.attach_triangulation(triangulation);

    std::ifstream input_file(parameters.mesh_filename.c_str());
    Assert(input_file, ExcFileNotOpen(parameters.mesh_filename.c_str()));

    grid_in.read_ucd(input_file);
  }

  dof_handler.clear();
  dof_handler.distribute_dofs(fe);

  // Size all of the fields.
  old_solution.reinit(dof_handler.n_dofs());
  current_solution.reinit(dof_handler.n_dofs());
  newton_initial_guess.reinit(dof_handler.n_dofs());
  right_hand_side.reinit(dof_handler.n_dofs());

  setup_system();

  VectorTools::interpolate(dof_handler, parameters.initial_condition, old_solution);
  current_solution = old_solution;
  newton_initial_guess = old_solution;

  output_results();

  Vector<double> newton_update(dof_handler.n_dofs());

  double time = 0;
  double next_output = time + parameters.output_step;

  newton_initial_guess = old_solution;
  while (time < parameters.final_time)
  {
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

template class Problem<2>;
template class Problem<3>;

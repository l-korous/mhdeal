#include "problem.h"

template <int dim>
Problem<dim>::Problem(Parameters<dim>& parameters, Equations<dim>& equations,
#ifdef HAVE_MPI
  parallel::distributed::Triangulation<dim>& triangulation
#else
  Triangulation<dim>& triangulation
#endif
) :
  mpi_communicator(MPI_COMM_WORLD),
  parameters(parameters),
  equations(equations),
  triangulation(triangulation),
  mapping(),
  fe(FE_DGT<dim>(parameters.polynomial_order), 1),
  dof_handler(triangulation),
  quadrature(parameters.quadrature_order),
  verbose_cout(std::cout, false)
{
}

template <int dim>
void Problem<dim>::setup_system()
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

template <int dim>
void Problem<dim>::assemble_system()
{
  // Number of DOFs per cell - we assume uniform polynomial order (we will only do h-adaptivity)
  const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

  // DOF indices both on the currently assembled element and the neighbor.
  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  // What values we need for the assembly.
  const UpdateFlags update_flags = update_values | update_q_points | update_JxW_values | update_gradients;

  // DOF indices both on the currently assembled element and the neighbor.
  FEValues<dim> fe_v(mapping, fe, quadrature, update_flags);

  // Local (cell) matrices and rhs - for the currently assembled element and the neighbor
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double> cell_rhs(dofs_per_cell);

  // Loop through all cells.
  for (typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
  {
    // Only assemble what belongs to this process.
    if (!cell->is_locally_owned())
      continue;

    cell_matrix = 0;
    cell_rhs = 0;

    fe_v.reinit(cell);
    cell->get_dof_indices(dof_indices);

    // Assemble the volumetric integrals.
    assemble_cell_term(fe_v, dof_indices, cell_matrix, cell_rhs);

    constraints.distribute_local_to_global(cell_matrix, cell_rhs, dof_indices, system_matrix, system_rhs);
  }

  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
}

template <int dim>
void
Problem<dim>::assemble_cell_term(const FEValues<dim> &fe_v, const std::vector<types::global_dof_index>& dof_indices, FullMatrix<double>& cell_matrix, Vector<double>& cell_rhs)
{
  const unsigned int dofs_per_cell = fe_v.dofs_per_cell;
  const unsigned int n_q_points = fe_v.n_quadrature_points;

    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      const unsigned int component_i = fe_v.get_fe().system_to_base_index(i).first.first;

      for (unsigned int j = 0; j < dofs_per_cell; ++j)
      {
        const unsigned int component_j = fe_v.get_fe().system_to_base_index(j).first.first;
        double matVal = 0;

        for (unsigned int q = 0; q < n_q_points; ++q)
          matVal += fe_v.shape_value_component(i, q, component_i) * fe_v.shape_value_component(j, q, component_j) * fe_v.JxW(q);

        if(i == j)
          cell_matrix(i, j) += 1.;
      }

      double rhsVal = 0;

      for (unsigned int q = 0; q < n_q_points; ++q)
        rhsVal += fe_v.shape_value_component(i, q, component_i) * fe_v.JxW(q);
      
      if (i == 2)
        cell_rhs(i) += 1.;
    }
}

template <int dim>
void
Problem<dim>::solve(TrilinosWrappers::MPI::Vector &newton_update)
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

template <int dim>
void Problem<dim>::output_results() const
{
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);

  // Solution components.
  data_out.add_data_vector(solution, equations.component_names(), DataOut<dim>::type_dof_data, equations.component_interpretation());

#ifdef HAVE_MPI
  // Subdomains.
  Vector<float> subdomain(triangulation.n_active_cells());
  for (unsigned int i = 0; i < subdomain.size(); ++i)
    subdomain(i) = triangulation.locally_owned_subdomain();
  data_out.add_data_vector(subdomain, "subdomain");
#endif

  data_out.build_patches(4);

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

template <int dim>
void Problem<dim>::run()
{
  // Preparations.
  setup_system();
  solution.reinit(locally_relevant_dofs, mpi_communicator);
  TrilinosWrappers::MPI::Vector newton_update(locally_relevant_dofs, mpi_communicator);

  system_matrix = 0;
  system_rhs = 0;
  assemble_system();

  // Outputs of algebraic stuff (optional - possible to be set in the Parameters class).
  if (parameters.output_matrix)
  {
    std::ofstream m;
    std::stringstream ssm;
    ssm << "matrix";
    m.open(ssm.str());
    system_matrix.print(m);
    m.close();
  }
  if (parameters.output_rhs)
  {
    std::ofstream r;
    std::stringstream ssr;
    ssr << "rhs";
    r.open(ssr.str());
    system_rhs.print(r, 10, false, false);
    r.close();
  }

  solve(solution);
  if (parameters.output_solution)
  {
    std::ofstream n;
    std::stringstream ssn;
    ssn << "solution";
    n.open(ssn.str());
    solution.print(n, 10, false, false);
    n.close();
  }

   output_results();
}

template class Problem<3>;

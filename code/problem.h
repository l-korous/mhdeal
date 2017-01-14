#include "util.h"
#include "equationsEuler.h"
#include "equationsMhd.h"
#include "parameters.h"
#include "initialCondition.h"
#include "boundaryConditions.h"

template <EquationsType equationsType, int dim>
class Problem
{
public:
  Problem(Parameters<dim>& parameters, Equations<equationsType, dim>& equations, Triangulation<dim>& triangulation, InitialCondition<equationsType, dim>& initial_condition, BoundaryConditions<equationsType, dim>& boundary_conditions);
  void run();

private:
  void setup_system();

  void assemble_system();
  void assemble_cell_term(const FEValues<dim> &fe_v, const std::vector<types::global_dof_index> &dofs);
  void assemble_face_term(const unsigned int face_no, const FEFaceValuesBase<dim> &fe_v,
    const FEFaceValuesBase<dim> &fe_v_neighbor,
    const std::vector<types::global_dof_index> &dofs,
    const std::vector<types::global_dof_index> &dofs_neighbor,
    const bool                       external_face,
    const unsigned int               boundary_id,
    const double                     face_diameter);

  std::pair<unsigned int, double> solve(Vector<double> &solution);

  void output_results() const;

  Triangulation<dim>& triangulation;

  const MappingQ1<dim> mapping;

  const FESystem<dim> fe;

  DoFHandler<dim> dof_handler;

  const QGauss<dim> quadrature;

  const QGauss<dim - 1> face_quadrature;

  Vector<double> old_solution;

  Vector<double> current_solution;

  Vector<double> newton_initial_guess;

  Equations<equationsType, dim>& equations;

  Vector<double> system_rhs;

  TrilinosWrappers::SparseMatrix system_matrix;

  Parameters<dim>& parameters;

  InitialCondition<equationsType, dim>& initial_condition;

  BoundaryConditions<equationsType, dim>& boundary_conditions;

  ConditionalOStream              verbose_cout;
};
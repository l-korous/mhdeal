#ifndef COMMON_H
#define COMMON_H

#include <fstream>
#include <vector>
#include <iostream>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/precondition_block.h>
#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/simple.h>
#include <deal.II/meshworker/loop.h>
#include <deal.II/base/point.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/lac/full_matrix.h>    
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/multithread_info.h>

using namespace dealii;

#include "definitions.h"

typedef unsigned int ui;
typedef double d;
typedef dealii::Tensor<1, DIM> dimVec;
typedef std::vector<dimVec> vecDimVec;
typedef dealii::Vector<d> vec;

#include "postProcessor.h"
#include "equationImplementation.h"
#include "mhdSolver.h"
#include "numericalFlux.h"
#include "boundaryConditions.h"
#include "initialSln.h"

#endif
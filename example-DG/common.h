#pragma region INCLUDES

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

#pragma endregion

#pragma region DIRECTIVES

#define DIM 3
#define DG_ORDER 1
#define INIT_REF_NUM 4
#define COMPONENT_COUNT 2

#pragma endregion

typedef dealii::Tensor<1, COMPONENT_COUNT> vec;
typedef dealii::Tensor<1, DIM> dimVec;
typedef std::vector<dimVec> vecDimVec;
typedef unsigned int ui;
typedef double d;

class Postprocessor : public DataPostprocessor < DIM >
{
public:
  Postprocessor();
  virtual void compute_derived_quantities_vector(const std::vector<Vector<d> > &uh, const std::vector<std::vector<Tensor<1, DIM> > > &duh, const std::vector<std::vector<Tensor<2, DIM> > > &dduh, const std::vector<Point<DIM> > &normals, const std::vector<Point<DIM> > &evaluation_points, const dealii::types::material_id mat_id, std::vector<Vector<d> > &computed_quantities) const;
  virtual std::vector<std::string> get_names() const;
  virtual std::vector < DataComponentInterpretation::DataComponentInterpretation > get_data_component_interpretation() const;
  virtual UpdateFlags get_needed_update_flags() const;
};

class MHDSolver
{
public:
  MHDSolver();
  void run();

private:
  void setup_system();
  void assemble_system();
  void solve(Vector<d> &solution);
  void solveOneStep(Vector<d> &solution);
  void outputResults() const;
  void add_markers(Triangulation<DIM>::cell_iterator cell);

  Triangulation<DIM>   triangulation;
  dealii::hp::FECollection<DIM> feCollection;
  dealii::hp::MappingCollection<DIM> mappingCollection;
  dealii::hp::QCollection<DIM> qCollection;
  dealii::hp::QCollection<DIM - 1> qCollectionFace;
  hp::DoFHandler<DIM>      dofHandler;
  ConstraintMatrix     hangingNodeConstraints;

  SparsityPattern      sparsityPattern;
  SparseMatrix<d> systemMatrix;

  Vector<d>       solution;
  Vector<d>       rightHandSide;

  typedef MeshWorker::DoFInfo<DIM> DoFInfo;
  typedef MeshWorker::IntegrationInfo<DIM> CellInfo;

  static void assembleVolumetric(DoFInfo &dinfo, CellInfo &info);
  static void assembleBoundaryEdge(DoFInfo &dinfo, CellInfo &info);
  static void assembleInternalEdge(DoFInfo &dinfo1, DoFInfo &dinfo2, CellInfo &info1, CellInfo &info2);
};

class EquationImplementation
{
public:

  static d matrixVolValue(ui component_i, ui component_j,
    d u_val, d v_val, vec Un_val, dimVec u_grad, dimVec v_grad,
    vecDimVec Un_grad, Point<DIM> quadPoint);

  static d matrixBoundaryEdgeValue(ui component_i, ui component_j,
    d u_val, d v_val, vec Un_val, dimVec u_grad, dimVec v_grad,
    vecDimVec Un_grad, vec U_bnd_val, Point<DIM> quadPoint, Point<DIM> normal);

  static d matrixInternalEdgeValue(ui component_i, ui component_j,
    d u_val, d v_val, vec Un_val, vec Un_valN, dimVec u_grad, dimVec v_grad,
    vecDimVec Un_grad, vecDimVec Un_gradN, bool u_N, bool v_N,
    Point<DIM> quadPoint, Point<DIM> normal);


  static d rhsVolValue(ui component_i,
    d v_val, vec Un_val, dimVec v_grad,
    vecDimVec Un_grad, Point<DIM> quadPoint);

  static d rhsBoundaryEdgeValue(ui component_i,
    d v_val, vec Un_val, dimVec v_grad,
    vecDimVec Un_grad, vec U_bnd_val, Point<DIM> quadPoint, Point<DIM> normal);

  static d rhsInternalEdgeValue(ui component_i,
    d v_val, vec Un_val, dimVec v_grad,
    vecDimVec Un_grad,
    d v_valN, vec Un_valN, dimVec v_gradN,
    vecDimVec Un_gradN, bool v_N, Point<DIM> quadPoint, Point<DIM> normal);
};

class NumFlux
{
public:
  virtual void calculate(vec U_L, vec U_R, Point<DIM> normal) = 0;
};

class NumFluxLaxFriedrichs
{
public:
  void calculate(vec U_L, vec U_R, Point<DIM> normal);
};
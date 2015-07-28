#pragma region INCLUDES

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
// Discontinuous space.
#include <deal.II/fe/fe_dgq.h>
// Solver.
#include <deal.II/lac/solver_richardson.h>
#include <deal.II/lac/precondition_block.h>

// MeshWorker.
#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/simple.h>
#include <deal.II/meshworker/loop.h>
#include <deal.II/base/point.h>
#include <vector>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/vector.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/lac/full_matrix.h>    
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_nothing.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>


#include <fstream>
#include <iostream>
#include <deal.II/numerics/solution_transfer.h>

#include <deal.II/base/work_stream.h>
#include <deal.II/base/multithread_info.h>

#include <iostream>
#include <fstream>

using namespace dealii;

#pragma endregion

#pragma region DIRECTIVES

#define DIM 3
#define DG_ORDER 1
#define INIT_REF_NUM 5
#define COMPONENT_COUNT 2

#pragma endregion

// @sect3{Equation data}
//
// First, we define a class describing the inhomogeneous boundary
// data. Since only its values are used, we implement value_list(), but
// leave all other functions of Function undefined.
class BoundaryValues : public Function < DIM >
{
public:
  BoundaryValues() {};
  virtual void value_list(const std::vector<Point<DIM> > &points,
    std::vector<double> &values,
    const unsigned int component = 0) const;
};

// Given the flow direction, the inflow boundary of the unit square
// $[0,1]^2$ are the right and the lower boundaries. We prescribe
// discontinuous boundary values 1 and 0 on the x-axis and value 0 on the
// right boundary. The values of this function on the outflow boundaries
// will not be used within the DG scheme.
void BoundaryValues::value_list(const std::vector<Point<DIM> > &points,
  std::vector<double> &values,
  const unsigned int) const
{
  Assert(values.size() == points.size(),
    ExcDimensionMismatch(values.size(), points.size()));

  for (unsigned int i = 0; i < values.size(); ++i)
  {
    if (points[i](0) < 0.5)
      values[i] = 1.;
    else
      values[i] = 0.;
  }
}

class FEProblem
{
public:
  FEProblem();
  void run();

private:
  void setup_system();
  void assemble_system();
  void solve(Vector<double> &solution);
  void output_results() const;
  void add_markers(Triangulation<DIM>::cell_iterator cell);

  Triangulation<DIM>   triangulation;
  dealii::hp::FECollection<DIM> feCollection;
  dealii::hp::MappingCollection<DIM> mappingCollection;
  dealii::hp::QCollection<DIM> qCollection;
  dealii::hp::QCollection<DIM - 1> qCollectionFace;
  hp::DoFHandler<DIM>      dof_handler;
  ConstraintMatrix     hanging_node_constraints;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double>       solution;
  Vector<double>       right_hand_side;

  typedef MeshWorker::DoFInfo<DIM> DoFInfo;
  typedef MeshWorker::IntegrationInfo<DIM> CellInfo;

  static void integrate_cell_term(DoFInfo &dinfo,
    CellInfo &info);
  static void integrate_boundary_term(DoFInfo &dinfo,
    CellInfo &info);
  static void integrate_face_term(DoFInfo &dinfo1,
    DoFInfo &dinfo2,
    CellInfo &info1,
    CellInfo &info2);
};


// We start with the constructor. The 1 in the constructor call of
// <code>fe</code> is the polynomial degree.
FEProblem::FEProblem()
  :
  dof_handler(triangulation)
{
  std::vector<const dealii::FiniteElement<DIM> *> fes;
  std::vector<unsigned int> multiplicities;

  // Spaces
  fes.push_back(new dealii::FE_DGQ<DIM>(DG_ORDER));
  multiplicities.push_back(COMPONENT_COUNT);

  feCollection.push_back(dealii::FESystem<DIM, DIM>(fes, multiplicities));

  mappingCollection.push_back(dealii::MappingQ<DIM>(1, true));

  qCollection.push_back(dealii::QGauss<DIM>(3 * DG_ORDER));
  qCollectionFace.push_back(dealii::QGauss<DIM - 1>(3 * DG_ORDER));
}

class Postprocessor : public DataPostprocessor < DIM >
{
public:
  Postprocessor();
  virtual void compute_derived_quantities_vector(const std::vector<Vector<double> > &uh, const std::vector<std::vector<Tensor<1, DIM> > > &duh, const std::vector<std::vector<Tensor<2, DIM> > > &dduh, const std::vector<Point<DIM> > &normals, const std::vector<Point<DIM> >                  &evaluation_points, const dealii::types::material_id mat_id, std::vector<Vector<double> >                    &computed_quantities) const;
  virtual std::vector<std::string> get_names() const;
  virtual std::vector < DataComponentInterpretation::DataComponentInterpretation > get_data_component_interpretation() const;
  virtual UpdateFlags get_needed_update_flags() const;
};

Postprocessor::Postprocessor() : DataPostprocessor<DIM>()
{}

void Postprocessor::compute_derived_quantities_vector(const std::vector<Vector<double> > &uh, const std::vector<std::vector<Tensor<1, DIM> > > &duh, const std::vector<std::vector<Tensor<2, DIM> > > &dduh, const std::vector<Point<DIM> > &normals, const std::vector<Point<DIM> > &evaluation_points, const dealii::types::material_id mat_id, std::vector<Vector<double> > &computed_quantities) const
{
  const unsigned int n_quadrature_points = uh.size();

  for (unsigned int q = 0; q < n_quadrature_points; ++q)
  {
    // Velocities
    for (unsigned int d = 0; d < COMPONENT_COUNT; ++d)
      computed_quantities[q](d) = uh[q](d);
  }
}

std::vector<std::string> Postprocessor::get_names() const
{
  std::vector<std::string> names;
  for (unsigned int d = 0; d < COMPONENT_COUNT; ++d)
  {
    std::stringstream ss;
    ss << "solution_";
    ss << std::to_string(d);
    names.push_back(ss.str());
  }
  return names;
}

std::vector<DataComponentInterpretation::DataComponentInterpretation> Postprocessor::get_data_component_interpretation() const
{
  std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretation;
  for (unsigned int d = 0; d < COMPONENT_COUNT; ++d)
    interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  return interpretation;
}

UpdateFlags Postprocessor::get_needed_update_flags() const
{
  return update_values | update_gradients | update_quadrature_points;
}

void FEProblem::setup_system()
{
  dof_handler.distribute_dofs(feCollection);

  dealii::DoFRenumbering::component_wise(dof_handler);

  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_flux_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);

  system_matrix.reinit(sparsity_pattern);
  right_hand_side.reinit(dof_handler.n_dofs());
  solution.reinit(dof_handler.n_dofs());
}

void FEProblem::assemble_system()
{
  MeshWorker::IntegrationInfoBox<DIM> info_box;

  info_box.initialize_gauss_quadrature(3 * DG_ORDER, 3 * DG_ORDER, 3 * DG_ORDER);

  // These are the types of values we need for integrating our system. They
  // are added to the flags used on cells, boundary and interior faces, as
  // well as interior neighbor faces, which is forced by the four @p true
  // values.
  info_box.initialize_update_flags();
  UpdateFlags update_flags = update_quadrature_points |
    update_values |
    update_gradients;
  info_box.add_update_flags(update_flags, true, true, true, true);

  // After preparing all data in <tt>info_box</tt>, we initialize the
  // FEValues objects in there.
  info_box.initialize(this->feCollection[0], this->mappingCollection[0]);

  // The object created so far helps us do the local integration on each
  // cell and face. Now, we need an object which receives the integrated
  // (local) data and forwards them to the assembler.
  MeshWorker::DoFInfo<DIM> dof_info(dof_handler);

  // Now, we have to create the assembler object and tell it, where to put
  // the local data. These will be our system matrix and the right hand
  // side.
  MeshWorker::Assembler::SystemSimple < SparseMatrix<double>, Vector<double> >
    assembler;
  assembler.initialize(system_matrix, right_hand_side);

  // Finally, the integration loop over all active cells (determined by the
  // first argument, which is an active iterator).
  //
  // As noted in the discussion when declaring the local integration
  // functions in the class declaration, the arguments expected by the
  // assembling integrator class are not actually function pointers. Rather,
  // they are objects that can be called like functions with a certain
  // number of arguments. Consequently, we could also pass objects with
  // appropriate operator() implementations here, or the result of std::bind
  // if the local integrators were, for example, non-static member
  // functions.
  MeshWorker::loop<DIM, DIM, MeshWorker::DoFInfo<DIM>, MeshWorker::IntegrationInfoBox<DIM> >
    (dof_handler.begin_active(), dof_handler.end(),
    dof_info, info_box,
    &FEProblem::integrate_cell_term,
    &FEProblem::integrate_boundary_term,
    &FEProblem::integrate_face_term,
    assembler);
}


// @sect4{The local integrators}

// These are the functions given to the MeshWorker::integration_loop()
// called just above. They compute the local contributions to the system
// matrix and right hand side on cells and faces.
void FEProblem::integrate_cell_term(DoFInfo &dinfo,
  CellInfo &info)
{
  // First, let us retrieve some of the objects used here from @p info. Note
  // that these objects can handle much more complex structures, thus the
  // access here looks more complicated than might seem necessary.
  const FEValuesBase<DIM> &fe_v = info.fe_values();
  FullMatrix<double> &local_matrix = dinfo.matrix(0).matrix;
  const std::vector<double> &JxW = fe_v.get_JxW_values();

  const unsigned int dofs_per_cell = info.finite_element().dofs_per_cell;
  std::vector<int> components(dofs_per_cell);

  for (unsigned int i = 0; i < dofs_per_cell; ++i)
  {
    components[i] = info.finite_element().system_to_component_index(i).first;
  }

  // With these objects, we continue local integration like always. First,
  // we loop over the quadrature points and compute the advection vector in
  // the current point.
  for (unsigned int point = 0; point < fe_v.n_quadrature_points; ++point)
  {
    Point<DIM> beta;
    beta(0) = -fe_v.quadrature_point(point)(1);
    beta(1) = fe_v.quadrature_point(point)(0);
    beta /= beta.norm();

    // We solve a homogeneous equation, thus no right hand side shows up
    // in the cell term.  What's left is integrating the matrix entries.
    for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
      for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j)
      {
        if (components[i] == components[j])
        {
          local_matrix(i, j) -= beta*fe_v.shape_grad(i, point)*
            fe_v.shape_value(j, point) *
            JxW[point];
        }
      }
  }
}

// Now the same for the boundary terms. Note that now we use FEValuesBase,
// the base class for both FEFaceValues and FESubfaceValues, in order to get
// access to normal vectors.
void FEProblem::integrate_boundary_term(DoFInfo &dinfo,
  CellInfo &info)
{
  const FEValuesBase<DIM> &fe_v = info.fe_values();
  FullMatrix<double> &local_matrix = dinfo.matrix(0).matrix;
  Vector<double> &local_vector = dinfo.vector(0).block(0);

  const std::vector<double> &JxW = fe_v.get_JxW_values();
  const std::vector<Point<DIM> > &normals = fe_v.get_normal_vectors();

  std::vector<double> g(fe_v.n_quadrature_points);

  static BoundaryValues boundary_function;
  boundary_function.value_list(fe_v.get_quadrature_points(), g);

  const unsigned int dofs_per_cell = info.finite_element().dofs_per_cell;
  std::vector<int> components(dofs_per_cell);

  for (unsigned int i = 0; i < dofs_per_cell; ++i)
  {
    components[i] = info.finite_element().system_to_component_index(i).first;
  }

  for (unsigned int point = 0; point < fe_v.n_quadrature_points; ++point)
  {
    Point<DIM> beta;
    beta(0) = -fe_v.quadrature_point(point)(1);
    beta(1) = fe_v.quadrature_point(point)(0);
    beta /= beta.norm();

    const double beta_n = beta * normals[point];
    if (beta_n > 0)
      for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
        for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j)
        {
          if (components[i] == components[j])
          {
            local_matrix(i, j) += beta_n *
              fe_v.shape_value(j, point) *
              fe_v.shape_value(i, point) *
              JxW[point];
          }
        }
    else
      for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
        local_vector(i) -= beta_n *
        g[point] *
        fe_v.shape_value(i, point) *
        JxW[point];
  }
}

// Finally, the interior face terms. The difference here is that we receive
// two info objects, one for each cell adjacent to the face and we assemble
// four matrices, one for each cell and two for coupling back and forth.
void FEProblem::integrate_face_term(DoFInfo &dinfo1,
  DoFInfo &dinfo2,
  CellInfo &info1,
  CellInfo &info2)
{
  // For quadrature points, weights, etc., we use the FEValuesBase object of
  // the first argument.
  const FEValuesBase<DIM> &fe_v = info1.fe_values();

  // For additional shape functions, we have to ask the neighbors
  // FEValuesBase.
  const FEValuesBase<DIM> &fe_v_neighbor = info2.fe_values();

  // Then we get references to the four local matrices. The letters u and v
  // refer to trial and test functions, respectively. The %numbers indicate
  // the cells provided by info1 and info2. By convention, the two matrices
  // in each info object refer to the test functions on the respective
  // cell. The first matrix contains the interior couplings of that cell,
  // while the second contains the couplings between cells.
  FullMatrix<double> &u1_v1_matrix = dinfo1.matrix(0, false).matrix;
  FullMatrix<double> &u2_v1_matrix = dinfo1.matrix(0, true).matrix;
  FullMatrix<double> &u1_v2_matrix = dinfo2.matrix(0, true).matrix;
  FullMatrix<double> &u2_v2_matrix = dinfo2.matrix(0, false).matrix;

  // Here, following the previous functions, we would have the local right
  // hand side vectors. Fortunately, the interface terms only involve the
  // solution and the right hand side does not receive any contributions.

  const std::vector<double> &JxW = fe_v.get_JxW_values();
  const std::vector<Point<DIM> > &normals = fe_v.get_normal_vectors();

  const unsigned int dofs_per_cell1 = info1.finite_element().dofs_per_cell;
  const unsigned int dofs_per_cell2 = info2.finite_element().dofs_per_cell;

  std::vector<int> components1(dofs_per_cell1);
  std::vector<int> components2(dofs_per_cell2);

  for (unsigned int i = 0; i < dofs_per_cell1; ++i)
  {
    components1[i] = info1.finite_element().system_to_component_index(i).first;
  }
  for (unsigned int i = 0; i < dofs_per_cell2; ++i)
  {
    components2[i] = info2.finite_element().system_to_component_index(i).first;
  }

  for (unsigned int point = 0; point < fe_v.n_quadrature_points; ++point)
  {
    Point<DIM> beta;
    beta(0) = -fe_v.quadrature_point(point)(1);
    beta(1) = fe_v.quadrature_point(point)(0);
    beta /= beta.norm();

    const double beta_n = beta * normals[point];
    if (beta_n > 0)
    {
      // This term we've already seen:
      for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
        for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j)
        {
          if (components1[i] == components1[j])
          {
            u1_v1_matrix(i, j) += beta_n *
              fe_v.shape_value(j, point) *
              fe_v.shape_value(i, point) *
              JxW[point];
          }
        }

      // We additionally assemble the term $(\beta\cdot n u,\hat
      // v)_{\partial \kappa_+}$,
      for (unsigned int k = 0; k < fe_v_neighbor.dofs_per_cell; ++k)
        for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j)
        {
          if (components1[j] == components2[k])
          {
            u1_v2_matrix(k, j) -= beta_n *
              fe_v.shape_value(j, point) *
              fe_v_neighbor.shape_value(k, point) *
              JxW[point];
          }
        }
    }
    else
    {
      // This one we've already seen, too:
      for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i)
        for (unsigned int l = 0; l < fe_v_neighbor.dofs_per_cell; ++l)
        {
          if (components1[i] == components2[l])
          {
            u2_v1_matrix(i, l) += beta_n *
              fe_v_neighbor.shape_value(l, point) *
              fe_v.shape_value(i, point) *
              JxW[point];
          }
        }

      // And this is another new one: $(\beta\cdot n \hat u,\hat
      // v)_{\partial \kappa_-}$:
      for (unsigned int k = 0; k < fe_v_neighbor.dofs_per_cell; ++k)
        for (unsigned int l = 0; l < fe_v_neighbor.dofs_per_cell; ++l)
        {
          if (components2[k] == components2[l])
          {
            u2_v2_matrix(k, l) -= beta_n *
              fe_v_neighbor.shape_value(l, point) *
              fe_v_neighbor.shape_value(k, point) *
              JxW[point];
          }
        }
    }
  }
}


// @sect3{All the rest}
//
// For this simple problem we use the simplest possible solver, called
// Richardson iteration, that represents a simple defect correction. This,
// in combination with a block SSOR preconditioner, that uses the special
// block matrix structure of system matrices arising from DG
// discretizations. The size of these blocks are the number of DoFs per
// cell. Here, we use a SSOR preconditioning as we have not renumbered the
// DoFs according to the flow field. If the DoFs are renumbered in the
// downstream direction of the flow, then a block Gauss-Seidel
// preconditioner (see the PreconditionBlockSOR class with relaxation=1)
// does a much better job.
void FEProblem::solve(Vector<double> &solution)
{
  SolverControl           solver_control(1000, 1e-12);
  SolverRichardson<>      solver(solver_control);

  // Here we create the preconditioner,
  PreconditionBlockSSOR<SparseMatrix<double> > preconditioner;

  // then assign the matrix to it and set the right block size:
  preconditioner.initialize(system_matrix, feCollection.max_dofs_per_cell());

  // After these preparations we are ready to start the linear solver.
  solver.solve(system_matrix, solution, right_hand_side,
    preconditioner);
}

void FEProblem::output_results() const
{
  Postprocessor postprocessor;
  DataOut<DIM, hp::DoFHandler<DIM> > data_out;
  data_out.attach_dof_handler(dof_handler);
  const DataOut<DIM, hp::DoFHandler<DIM> >::DataVectorType data_vector_type = DataOut<DIM, hp::DoFHandler<DIM> >::type_dof_data;
  data_out.add_data_vector(solution, postprocessor);
  data_out.build_patches();
  std::string filename = "solution.vtk";
  std::ofstream output(filename.c_str());
  data_out.write_vtk(output);
}


// The following <code>run</code> function is similar to previous examples.
void FEProblem::run()
{
  GridGenerator::hyper_cube(triangulation);

  triangulation.refine_global(INIT_REF_NUM);

  deallog << "Number of active cells:       "
    << triangulation.n_active_cells()
    << std::endl;

  setup_system();

  deallog << "Number of degrees of freedom: "
    << dof_handler.n_dofs()
    << std::endl;

  assemble_system();
  solve(solution);

  output_results();
}


// The following <code>main</code> function is similar to previous examples as
// well, and need not be commented on.
int main()
{
  try
  {
    FEProblem feProblem;
    feProblem.run();
  }
  catch (std::exception &exc)
  {
    std::cerr << std::endl << std::endl
      << "----------------------------------------------------"
      << std::endl;
    std::cerr << "Exception on processing: " << std::endl
      << exc.what() << std::endl
      << "Aborting!" << std::endl
      << "----------------------------------------------------"
      << std::endl;
    return 1;
  }
  catch (...)
  {
    std::cerr << std::endl << std::endl
      << "----------------------------------------------------"
      << std::endl;
    std::cerr << "Unknown exception!" << std::endl
      << "Aborting!" << std::endl
      << "----------------------------------------------------"
      << std::endl;
    return 1;
  };

  return 0;
}
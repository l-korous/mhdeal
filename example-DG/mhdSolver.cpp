#include "common.h"

typedef EquationImplementation Eq;

MHDSolver::MHDSolver()
  :
  dofHandler(triangulation)
{
  std::vector<const dealii::FiniteElement<DIM> *> fes;
  std::vector<ui> multiplicities;

  fes.push_back(new dealii::FE_DGQ<DIM>(DG_ORDER));
  multiplicities.push_back(COMPONENT_COUNT);
  feCollection.push_back(dealii::FESystem<DIM, DIM>(fes, multiplicities));

  mappingCollection.push_back(dealii::MappingQ<DIM>(1, true));

  qCollection.push_back(dealii::QGauss<DIM>(3 * DG_ORDER));
  qCollectionFace.push_back(dealii::QGauss<DIM - 1>(3 * DG_ORDER));
}

void MHDSolver::setup_system()
{
  dofHandler.distribute_dofs(feCollection);

  dealii::DoFRenumbering::component_wise(dofHandler);

  DynamicSparsityPattern dsp(dofHandler.n_dofs());
  DoFTools::make_flux_sparsity_pattern(dofHandler, dsp);
  sparsityPattern.copy_from(dsp);

  systemMatrix.reinit(sparsityPattern);
  rightHandSide.reinit(dofHandler.n_dofs());
  solution.reinit(dofHandler.n_dofs());
}

void MHDSolver::assemble_system()
{
  MeshWorker::IntegrationInfoBox<DIM> info_box;

  // \todo This is wrong probably.
  info_box.initialize_gauss_quadrature(3 * DG_ORDER, 3 * DG_ORDER, 3 * DG_ORDER);

  info_box.initialize_update_flags();
  UpdateFlags update_flags = update_quadrature_points | update_values | update_gradients;
  info_box.add_update_flags(update_flags, true, true, true, true);

  // \todo What about multiple FEs in feCollection?
  info_box.initialize(this->feCollection[0], this->mappingCollection[0]);

  // \todo This has to be done properly for hpDoFHandler (varying number of DOFs per cell)
  MeshWorker::DoFInfo<DIM> dof_info(dofHandler);

  MeshWorker::Assembler::SystemSimple < SparseMatrix<d>, Vector<d> > assembler;
  assembler.initialize(systemMatrix, rightHandSide);

  // \todo This comes from tutorial, it may need some adjustment.
  MeshWorker::loop<DIM, DIM, MeshWorker::DoFInfo<DIM>, MeshWorker::IntegrationInfoBox<DIM> >
    (dofHandler.begin_active(), dofHandler.end(), dof_info, info_box, &MHDSolver::assembleVolumetric,
    &MHDSolver::assembleBoundaryEdge, &MHDSolver::assembleInternalEdge, assembler);
}

void MHDSolver::assembleVolumetric(DoFInfo &dinfo,
  CellInfo &info)
{
  const FEValuesBase<DIM> &fe_v = info.fe_values();
  FullMatrix<d> &local_matrix = dinfo.matrix(0).matrix;
  const std::vector<d> &JxW = fe_v.get_JxW_values();

  const ui dofs_per_cell = info.finite_element().dofs_per_cell;

  // Components
  std::vector<int> components(dofs_per_cell);
  for (ui i = 0; i < dofs_per_cell; ++i)
    components[i] = info.finite_element().system_to_component_index(i).first;

  for (ui point = 0; point < fe_v.n_quadrature_points; ++point)
  {
    for (ui i = 0; i < fe_v.dofs_per_cell; ++i)
      for (ui j = 0; j < fe_v.dofs_per_cell; ++j)
      {
        local_matrix(i, j) += JxW[point] * Eq::matrixVolValue(components[j], components[i],
          fe_v.shape_value(j, point), fe_v.shape_value(i, point),
          vec(), fe_v.shape_grad(j, point), fe_v.shape_grad(i, point),
          vecDimVec(), fe_v.quadrature_point(point));
      }
  }
}

// Now the same for the boundary terms. Note that now we use FEValuesBase,
// the base class for both FEFaceValues and FESubfaceValues, in order to get
// access to normal vectors.
void MHDSolver::assembleBoundaryEdge(DoFInfo &dinfo,
  CellInfo &info)
{
  const FEValuesBase<DIM> &fe_v = info.fe_values();
  FullMatrix<d> &local_matrix = dinfo.matrix(0).matrix;
  Vector<d> &local_vector = dinfo.vector(0).block(0);

  const std::vector<d> &JxW = fe_v.get_JxW_values();
  const std::vector<Point<DIM> > &normals = fe_v.get_normal_vectors();

  const ui dofs_per_cell = info.finite_element().dofs_per_cell;
  std::vector<int> components(dofs_per_cell);

  for (ui i = 0; i < dofs_per_cell; ++i)
  {
    components[i] = info.finite_element().system_to_component_index(i).first;
  }

  for (ui point = 0; point < fe_v.n_quadrature_points; ++point)
  {
    for (ui i = 0; i < fe_v.dofs_per_cell; ++i)
    {
      for (ui j = 0; j < fe_v.dofs_per_cell; ++j)
      {
        local_matrix(i, j) += JxW[point] * Eq::matrixBoundaryEdgeValue(components[j], components[i],
          fe_v.shape_value(j, point), fe_v.shape_value(i, point),
          vec(), fe_v.shape_grad(j, point), fe_v.shape_grad(i, point),
          vecDimVec(), vec(), fe_v.quadrature_point(point), normals[point]);
      }
      local_vector(i) += JxW[point] * Eq::rhsBoundaryEdgeValue(components[i],
        fe_v.shape_value(i, point),
        vec(), fe_v.shape_grad(i, point),
        vecDimVec(), vec(), fe_v.quadrature_point(point), normals[point]);
    }
  }
}

// Finally, the interior face terms. The difference here is that we receive
// two info objects, one for each cell adjacent to the face and we assemble
// four matrices, one for each cell and two for coupling back and forth.
void MHDSolver::assembleInternalEdge(DoFInfo &dinfo1,
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
  FullMatrix<d> &u1_v1_matrix = dinfo1.matrix(0, false).matrix;
  FullMatrix<d> &u2_v1_matrix = dinfo1.matrix(0, true).matrix;
  FullMatrix<d> &u1_v2_matrix = dinfo2.matrix(0, true).matrix;
  FullMatrix<d> &u2_v2_matrix = dinfo2.matrix(0, false).matrix;

  // Here, following the previous functions, we would have the local right
  // hand side vectors. Fortunately, the interface terms only involve the
  // solution and the right hand side does not receive any contributions.

  const std::vector<d> &JxW = fe_v.get_JxW_values();
  const std::vector<Point<DIM> > &normals = fe_v.get_normal_vectors();

  const ui dofs_per_cell1 = info1.finite_element().dofs_per_cell;
  const ui dofs_per_cell2 = info2.finite_element().dofs_per_cell;

  std::vector<int> components1(dofs_per_cell1);
  std::vector<int> components2(dofs_per_cell2);

  for (ui i = 0; i < dofs_per_cell1; ++i)
  {
    components1[i] = info1.finite_element().system_to_component_index(i).first;
  }
  for (ui i = 0; i < dofs_per_cell2; ++i)
  {
    components2[i] = info2.finite_element().system_to_component_index(i).first;
  }

  for (ui point = 0; point < fe_v.n_quadrature_points; ++point)
  {
    for (ui i = 0; i < fe_v.dofs_per_cell; ++i)
    {
      for (ui j = 0; j < fe_v.dofs_per_cell; ++j)
      {
        u1_v1_matrix(i, j) += JxW[point] * Eq::matrixInternalEdgeValue(components1[j], components1[i],
          fe_v.shape_value(j, point), fe_v.shape_value(i, point),
          vec(), vec(), fe_v.shape_grad(j, point), fe_v.shape_grad(i, point),
          vecDimVec(), vecDimVec(), false, false, fe_v.quadrature_point(point), normals[point]);
      }
    }

    for (ui k = 0; k < fe_v_neighbor.dofs_per_cell; ++k)
    {
      for (ui j = 0; j < fe_v.dofs_per_cell; ++j)
      {
        u1_v2_matrix(k, j) += JxW[point] * Eq::matrixInternalEdgeValue(components1[j], components1[k],
          fe_v.shape_value(j, point), fe_v_neighbor.shape_value(k, point),
          vec(), vec(), fe_v.shape_grad(j, point), fe_v_neighbor.shape_grad(k, point),
          vecDimVec(), vecDimVec(), false, true, fe_v.quadrature_point(point), normals[point]);
      }
    }

    for (ui i = 0; i < fe_v.dofs_per_cell; ++i)
    {
      for (ui l = 0; l < fe_v_neighbor.dofs_per_cell; ++l)
      {
        u2_v1_matrix(i, l) += JxW[point] * Eq::matrixInternalEdgeValue(components1[l], components1[i],
          fe_v_neighbor.shape_value(l, point), fe_v.shape_value(i, point),
          vec(), vec(), fe_v_neighbor.shape_grad(l, point), fe_v.shape_grad(i, point),
          vecDimVec(), vecDimVec(), true, false, fe_v.quadrature_point(point), normals[point]);
      }
    }

    for (ui k = 0; k < fe_v_neighbor.dofs_per_cell; ++k)
    {
      for (ui l = 0; l < fe_v_neighbor.dofs_per_cell; ++l)
      {
        u2_v2_matrix(k, l) += JxW[point] * Eq::matrixInternalEdgeValue(components1[l], components1[k],
          fe_v_neighbor.shape_value(l, point), fe_v_neighbor.shape_value(k, point),
          vec(), vec(), fe_v_neighbor.shape_grad(l, point), fe_v_neighbor.shape_grad(k, point),
          vecDimVec(), vecDimVec(), true, true, fe_v.quadrature_point(point), normals[point]);
      }
    }
  }
}

void MHDSolver::solve(Vector<d> &solution)
{
  SolverControl           solver_control(1000, 1e-12);
  SolverRichardson<>      solver(solver_control);

  // Here we create the preconditioner,
  PreconditionBlockSSOR<SparseMatrix<d> > preconditioner;

  // then assign the matrix to it and set the right block size:
  preconditioner.initialize(systemMatrix, feCollection.max_dofs_per_cell());

  // After these preparations we are ready to start the linear solver.
  solver.solve(systemMatrix, solution, rightHandSide,
    preconditioner);
}

void MHDSolver::outputResults() const
{
  Postprocessor postprocessor;
  DataOut<DIM, hp::DoFHandler<DIM> > data_out;
  data_out.attach_dof_handler(dofHandler);
  const DataOut<DIM, hp::DoFHandler<DIM> >::DataVectorType data_vector_type = DataOut<DIM, hp::DoFHandler<DIM> >::type_dof_data;
  data_out.add_data_vector(solution, postprocessor);
  data_out.build_patches();
  std::string filename = "solution.vtk";
  std::ofstream output(filename.c_str());
  data_out.write_vtk(output);
}

void MHDSolver::run()
{
  GridGenerator::hyper_cube(triangulation);

  triangulation.refine_global(INIT_REF_NUM);

  deallog << "Number of active cells:       "
    << triangulation.n_active_cells()
    << std::endl;

  setup_system();

  deallog << "Number of degrees of freedom: "
    << dofHandler.n_dofs()
    << std::endl;

  assemble_system();
  solve(solution);

  outputResults();
}
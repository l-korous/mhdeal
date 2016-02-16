#include <fstream>
#include <vector>
#include <iostream>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/base/timer.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/solution_transfer.h>

#include "mhdSolver.h"
#include "equationImplementation.h"
#include "postProcessor.h"
#include "initialSln.h"
#include "boundaryConditions.h"
#include "numericalFlux.h"



using namespace dealii;
using namespace std;

typedef EquationImplementation Eq;

d MHDSolver::A[3][COMPONENT_COUNT][COMPONENT_COUNT];
Vector<d> MHDSolver::slnPrev;
Vector<d> MHDSolver::slnLin;
// For initial conditions.
Vector<d> MHDSolver::slnUtil;
NumFlux* numFlux;
DirichletBoundaryCondition bc;

MHDSolver::MHDSolver()
  :
  feSystem(dealii::FE_DGQ<DIM>(DG_ORDER), COMPONENT_COUNT),
  dofHandler(triangulation),
  mapping(),
  quad(2 * DG_ORDER),
  quadFace(2 * DG_ORDER)
{
  numFlux = new NumFluxHLLD();
  
  for(unsigned int i=0;i<COMPONENT_COUNT;i++){
    for(unsigned int j=0;j<COMPONENT_COUNT;j++){
      A[0][i][j]=A[1][i][j]=A[2][i][j]=0.0;
    }
  }
}

void MHDSolver::setup_system()
{
  dofHandler.distribute_dofs(feSystem);

  dealii::DoFRenumbering::component_wise(dofHandler);

  DynamicSparsityPattern dsp(dofHandler.n_dofs());
  DoFTools::make_flux_sparsity_pattern(dofHandler, dsp);
  sparsityPattern.copy_from(dsp);

  systemMatrix.reinit(sparsityPattern);
  rightHandSide.reinit(dofHandler.n_dofs());
  solution.reinit(dofHandler.n_dofs());
  slnPrev.reinit(dofHandler.n_dofs());
  slnLin.reinit(dofHandler.n_dofs());
  slnUtil.reinit(dofHandler.n_dofs());
}

void MHDSolver::assemble_system(bool firstIteration)
{
  MeshWorker::IntegrationInfoBox<DIM> info_box;

  // \todo This is wrong probably.
  info_box.initialize_gauss_quadrature(3 * DG_ORDER, 3 * DG_ORDER, 3 * DG_ORDER);

  AnyData solution_data;
  solution_data.add(&slnPrev, "solution");
  info_box.cell_selector.add("solution", true, true, false);
  info_box.boundary_selector.add("solution", true, false, false);
  info_box.face_selector.add("solution", true, false, false);
  
  solution_data.add(&slnLin, "linearization");
  info_box.cell_selector.add("linearization", true, false, false);
  info_box.boundary_selector.add("linearization", true, false, false);
  info_box.face_selector.add("linearization", true, false, false);

  info_box.initialize_update_flags();
  UpdateFlags update_flags = update_quadrature_points | update_values | update_gradients;
  info_box.add_update_flags(update_flags, true, true, true, true);

  // \todo What about multiple FEs in feCollection?
  info_box.initialize(feSystem, mapping, solution_data, rightHandSide);

  // \todo This has to be done properly for hpDoFHandler (varying number of DOFs per cell)
  MeshWorker::DoFInfo<DIM> dof_info(dofHandler);

  if (!firstIteration)
  {
    systemMatrix.reinit(sparsityPattern);
    rightHandSide.reinit(dofHandler.n_dofs());
  }
  MeshWorker::Assembler::SystemSimple < SparseMatrix<d>, Vector<d> > assembler;
  assembler.initialize(systemMatrix, rightHandSide);

  // \todo This comes from tutorial, it may need some adjustment.
  MeshWorker::loop<DIM, DIM, MeshWorker::DoFInfo<DIM>, MeshWorker::IntegrationInfoBox<DIM> >
    (dofHandler.begin_active(), dofHandler.end(), dof_info, info_box, &MHDSolver::assembleVolumetric,
    &MHDSolver::saveInfoBoundaryEdge, &MHDSolver::assembleInternalEdge, assembler);

  MeshWorker::loop<DIM, DIM, MeshWorker::DoFInfo<DIM>, MeshWorker::IntegrationInfoBox<DIM> >
    (dofHandler.begin_active(), dofHandler.end(), dof_info, info_box, &MHDSolver::assembleVolumetricEmpty,
    &MHDSolver::assembleBoundaryEdge, &MHDSolver::assembleInternalEdgeEmpty, assembler);

  periodicBdr.clear();
}

void MHDSolver::assembleVolumetric(DoFInfo &dinfo,
  CellInfo &info)
{
  double rhs[COMPONENT_COUNT];
  const FEValuesBase<DIM> &fe_v = info.fe_values();
  dealii::FullMatrix<d> &local_matrix = dinfo.matrix(0).matrix;
  Vector<d> &local_vector = dinfo.vector(0).block(0);
  const std::vector<d> &JxW = fe_v.get_JxW_values();

  const ui dofs_per_cell = info.finite_element().dofs_per_cell;

  std::vector<dealii::Vector<double> > prev_values(fe_v.n_quadrature_points, dealii::Vector<double>(COMPONENT_COUNT));
  std::vector<std::vector<Tensor<1,DIM> > > prev_grads(fe_v.n_quadrature_points, std::vector<Tensor<1,DIM> >(COMPONENT_COUNT));
  std::vector<dealii::Vector<double> > lin_values(fe_v.n_quadrature_points, dealii::Vector<double>(COMPONENT_COUNT));
  for (ui i = 0; i < COMPONENT_COUNT; i++)
    for (ui point = 0; point < fe_v.n_quadrature_points; ++point){
      prev_values[point][i] = info.values[0][i][point];
      prev_grads[point][i] = info.gradients[0][i][point];
      lin_values[point][i] = info.values[1][i][point];
    }

  // Components
  std::vector<int> components(dofs_per_cell);
  for (ui i = 0; i < dofs_per_cell; ++i)
    components[i] = info.finite_element().system_to_component_index(i).first;

  for (ui point = 0; point < fe_v.n_quadrature_points; ++point)
  {
    JacobiM(A,lin_values[point]);
    for (ui i = 0; i < fe_v.dofs_per_cell; ++i)
      for (ui j = 0; j < fe_v.dofs_per_cell; ++j){ // u.v + dt * sum_d A_d * u * dv/dx_d
        if (components[i]==components[j])
          local_matrix(i, j) = JxW[point]*fe_v.shape_value(i, point)*fe_v.shape_value(j, point);
        for(ui d=0;d<DIM;d++)
          local_matrix(i, j) += JxW[point]* DELTA_T*
           A[d][components[i]][components[j]]*fe_v.shape_value(i, point)*fe_v.shape_grad(j, point)[d];
      }

    JacobiM(A,prev_values[point]);
    for(ui i=0;i<COMPONENT_COUNT_T;i++){  //  sum_d dF_d(u_old)/dx_d
      rhs[i]=A[0][i][0]*prev_grads[point][0][0];
      for(ui d=1;d<DIM;d++)
        rhs[i]+=A[d][i][0]*prev_grads[point][0][d];
      for(ui d=0;d<DIM;d++)
        for(ui j=1;j<COMPONENT_COUNT;j++)
          rhs[i]+=A[d][i][j]*prev_grads[point][j][d];
    }
    for (ui i = 0; i < fe_v.dofs_per_cell; ++i){  // u_old - dt * sum_d dF_d(u_old)/dx_d
      local_vector(i) += JxW[point] * (prev_values[point][components[i]] - DELTA_T*
            rhs[components[i]] );
    }
  }
}

bool eqiv(double a, double b)
{
  return fabs(a-b) < 1e-10;
}

std::vector<dealii::Vector<double> > MHDSolver::findCorrespondingInfo(dealii::Point<DIM> myCenter)
{
  for(PeriodicBdrInfo info : periodicBdr)
  {
    if(myCenter[0] != info.center[0] && myCenter[1] == info.center[1] && myCenter[2] == info.center[2])
    {
//      cout << "for " << myCenter << " found " << info.center << endl;
      return info.prev_values;
    }
  }
  assert(false);
}

void MHDSolver::saveInfoBoundaryEdge(DoFInfo &dinfo, CellInfo &info)
{
}

void MHDSolver::assembleBoundaryEdge(DoFInfo &dinfo, CellInfo &info)
{
  const FEValuesBase<DIM> &fe_v = info.fe_values();
  dealii::FullMatrix<d> &local_matrix = dinfo.matrix(0).matrix;
  Vector<d> &local_vector = dinfo.vector(0).block(0);

  const std::vector<d> &JxW = fe_v.get_JxW_values();
  const std::vector<Point<DIM> > &normals = fe_v.get_normal_vectors();

  const ui dofs_per_cell = info.finite_element().dofs_per_cell;

  // Previous values.
  std::vector<dealii::Vector<double> > prev_values(fe_v.n_quadrature_points, dealii::Vector<double>(COMPONENT_COUNT));
  for (ui i = 0; i < COMPONENT_COUNT; i++)
    for (ui point = 0; point < fe_v.n_quadrature_points; ++point)
      prev_values[point][i] = info.values[0][i][point];

  std::vector<dealii::Vector<double> > prev_values_oposite_elem;//(fe_v.n_quadrature_points, dealii::Vector<double>(COMPONENT_COUNT));

  // Components.
  std::vector<int> components(dofs_per_cell);
  for (ui i = 0; i < dofs_per_cell; ++i)
    components[i] = info.finite_element().system_to_component_index(i).first;


  for (ui point = 0; point < fe_v.n_quadrature_points; ++point)
  {
    for (ui i = 0; i < fe_v.dofs_per_cell; ++i)
    {
      for (ui j = 0; j < fe_v.dofs_per_cell; ++j)
      {
          local_matrix(i, j) += JxW[point] * Eq::matrixBoundaryEdgeValue(components[j], components[i],
            fe_v.shape_value(j, point), fe_v.shape_value(i, point),
            prev_values[point], fe_v.shape_grad(j, point), fe_v.shape_grad(i, point),
            vecDimVec(), vec(), fe_v.quadrature_point(point), normals[point], numFlux, &bc, dinfo.face->boundary_id());
      }
        local_vector(i) += JxW[point] * Eq::rhsBoundaryEdgeValue(components[i],
          fe_v.shape_value(i, point), prev_values[point], fe_v.shape_grad(i, point),
          vecDimVec(), vec(), fe_v.quadrature_point(point), normals[point], numFlux, &bc, dinfo.face->boundary_id());
    }
  }
}


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
  dealii::FullMatrix<d> &u1_v1_matrix = dinfo1.matrix(0, false).matrix;
  dealii::FullMatrix<d> &u2_v1_matrix = dinfo1.matrix(0, true).matrix;
  dealii::FullMatrix<d> &u1_v2_matrix = dinfo2.matrix(0, true).matrix;
  dealii::FullMatrix<d> &u2_v2_matrix = dinfo2.matrix(0, false).matrix;

  Vector<d> &v1_vector = dinfo1.vector(0).block(0);
  Vector<d> &v2_vector = dinfo2.vector(0).block(0);

  // Here, following the previous functions, we would have the local right
  // hand side vectors. Fortunately, the interface terms only involve the
  // solution and the right hand side does not receive any contributions.

  const std::vector<d> &JxW = fe_v.get_JxW_values();
  const std::vector<Point<DIM> > &normals = fe_v.get_normal_vectors();

  const ui dofs_per_cell1 = info1.finite_element().dofs_per_cell;
  const ui dofs_per_cell2 = info2.finite_element().dofs_per_cell;

  // Previous values.
  std::vector<dealii::Vector<double> > prev_values1(fe_v.n_quadrature_points, dealii::Vector<double>(COMPONENT_COUNT));
  for (ui i = 0; i < COMPONENT_COUNT; i++)
    for (ui point = 0; point < fe_v.n_quadrature_points; ++point)
      prev_values1[point][i] = info1.values[1][i][point];

  std::vector<dealii::Vector<double> > prev_values2(fe_v_neighbor.n_quadrature_points, dealii::Vector<double>(COMPONENT_COUNT));
  for (ui i = 0; i < COMPONENT_COUNT; i++)
    for (ui point = 0; point < fe_v_neighbor.n_quadrature_points; ++point)
      prev_values2[point][i] = info2.values[1][i][point];

  // Components.
  std::vector<int> components1(dofs_per_cell1);
  for (ui i = 0; i < dofs_per_cell1; ++i)
    components1[i] = info1.finite_element().system_to_component_index(i).first;

  std::vector<int> components2(dofs_per_cell2);
  for (ui i = 0; i < dofs_per_cell2; ++i)
    components2[i] = info2.finite_element().system_to_component_index(i).first;

  for (ui point = 0; point < fe_v.n_quadrature_points; ++point)
  {
    for (ui i = 0; i < fe_v.dofs_per_cell; ++i)
    {
      for (ui j = 0; j < fe_v.dofs_per_cell; ++j)
      {
        u1_v1_matrix(i, j) += JxW[point] * Eq::matrixInternalEdgeValue(components1[j], components1[i],
          fe_v.shape_value(i, point), fe_v.shape_value(j, point),
          prev_values1[point], prev_values2[point], fe_v.shape_grad(j, point), fe_v.shape_grad(i, point),
          vecDimVec(), vecDimVec(), false, false, fe_v.quadrature_point(point), normals[point], numFlux);
      }
    }

    for (ui k = 0; k < fe_v_neighbor.dofs_per_cell; ++k)
    {
      for (ui j = 0; j < fe_v.dofs_per_cell; ++j)
      {
        u1_v2_matrix(k, j) += JxW[point] * Eq::matrixInternalEdgeValue(components1[j], components2[k],
          fe_v.shape_value(k, point), fe_v_neighbor.shape_value(j, point),
          prev_values1[point], prev_values2[point], fe_v.shape_grad(j, point), fe_v_neighbor.shape_grad(k, point),
          vecDimVec(), vecDimVec(), true, false, fe_v.quadrature_point(point), normals[point], numFlux);
      }
    }

    for (ui i = 0; i < fe_v.dofs_per_cell; ++i)
    {
      for (ui l = 0; l < fe_v_neighbor.dofs_per_cell; ++l)
      {
        u2_v1_matrix(i, l) += JxW[point] * Eq::matrixInternalEdgeValue(components2[l], components1[i],
          fe_v_neighbor.shape_value(i, point), fe_v.shape_value(l, point),
          prev_values1[point], prev_values2[point], fe_v_neighbor.shape_grad(l, point), fe_v.shape_grad(i, point),
          vecDimVec(), vecDimVec(), false, true, fe_v.quadrature_point(point), normals[point], numFlux);
      }
    }

    for (ui k = 0; k < fe_v_neighbor.dofs_per_cell; ++k)
    {
      for (ui l = 0; l < fe_v_neighbor.dofs_per_cell; ++l)
      {
        u2_v2_matrix(k, l) += JxW[point] * Eq::matrixInternalEdgeValue(components2[l], components2[k],
          fe_v_neighbor.shape_value(k, point), fe_v_neighbor.shape_value(l, point),
          prev_values1[point], prev_values2[point], fe_v_neighbor.shape_grad(l, point), fe_v_neighbor.shape_grad(k, point),
          vecDimVec(), vecDimVec(), true, true, fe_v.quadrature_point(point), normals[point], numFlux);
      }
    }

    for (ui i = 0; i < fe_v.dofs_per_cell; ++i)
    {
      v1_vector(i) += JxW[point] * Eq::rhsInternalEdgeValue(components1[i], fe_v.shape_value(i, point),
                      fe_v.shape_grad(i, point), false, prev_values1[point], vecDimVec(), prev_values2[point],
                      vecDimVec(), fe_v.quadrature_point(point), normals[point], numFlux);
    }

    for (ui l = 0; l < fe_v_neighbor.dofs_per_cell; ++l)
    {
      assert(fe_v.quadrature_point(point) == fe_v_neighbor.quadrature_point(point));
      v2_vector(l) += JxW[point] * Eq::rhsInternalEdgeValue(components2[l], fe_v_neighbor.shape_value(l, point),
                      fe_v_neighbor.shape_grad(l, point), true, prev_values1[point], vecDimVec(), prev_values2[point],
                      vecDimVec(), fe_v.quadrature_point(point), normals[point], numFlux);
    }
  }
}

void MHDSolver::solve(Vector<d> &solution)
{
  if (PRINT_ALGEBRA)
  {
    std::cout << "  Printing system... " << std::endl;

    std::string matrix_file = "Matrix.txt";
    std::string rhs_file = "Rhs.txt";

    std::ofstream matrix_out(matrix_file);
    std::ofstream rhs_out(rhs_file);

    matrix_out << std::fixed;
    rhs_out << std::fixed;
    matrix_out << std::setprecision(6);
    rhs_out << std::setprecision(6);
    systemMatrix.print(matrix_out);
    rightHandSide.print(rhs_out, 6, false, false);

    matrix_out.close();
    rhs_out.close();
  }
  /*
  SolverControl           solver_control(1000, 1e-12);
  SolverRichardson<>      solver(solver_control);

  // Here we create the preconditioner,
  PreconditionBlockSSOR<SparseMatrix<d> > preconditioner;

  // then assign the matrix to it and set the right block size:
  preconditioner.initialize(systemMatrix, feSystem.dofs_per_cell);

  // After these preparations we are ready to start the linear solver.
  solver.solve(systemMatrix, solution, rightHandSide, preconditioner);
  */

  dealii::SparseDirectUMFPACK solver;

  solver.initialize(systemMatrix);

  solver.Tvmult(solution, rightHandSide);
}

void MHDSolver::outputResults(ui timeStep, d currentTime) const
{
  Postprocessor postprocessor;
  DataOut<DIM, DoFHandler<DIM> > data_out;
  data_out.set_flags(DataOutBase::VtkFlags(std::numeric_limits<d>::min(), std::numeric_limits<d>::min(), false, DataOutBase::VtkFlags::no_compression));
  data_out.attach_dof_handler(dofHandler);
  const DataOut<DIM, DoFHandler<DIM> >::DataVectorType data_vector_type = DataOut<DIM, DoFHandler<DIM> >::type_dof_data;
  data_out.add_data_vector(slnPrev, postprocessor);
  data_out.build_patches(mapping);
  std::stringstream ss;
  ss << "solution-";
  ss << timeStep;
  ss << ".vtk";
  std::ofstream output(ss.str());
  data_out.write_vtk(output);
}

void MHDSolver::add_markers(Triangulation<DIM>::cell_iterator cell)
{
}

void MHDSolver::run()
{
  GridGenerator::subdivided_hyper_rectangle(triangulation, std::vector<ui>({ INIT_REF_NUM_X, INIT_REF_NUM_Y, INIT_REF_NUM_Z}), p1, p4);
  
  std::string tria_file = "Tria.vtk";
  std::ofstream tria_out(tria_file);
  GridOut().write_vtk(triangulation, tria_out);

  Triangulation<DIM>::cell_iterator
    cell = triangulation.begin(),
    endc = triangulation.end();
  for (; cell != endc; ++cell)
  {
    this->add_markers(cell);
  }

  deallog << "Number of active cells:       "
    << triangulation.n_active_cells()
    << std::endl;

  setup_system();

  deallog << "Number of degrees of freedom: "
    << dofHandler.n_dofs()
    << std::endl;

  // Initial sln.
  VectorFunctionFromScalarFunctionObject<DIM> initialSlnRho(InitialSlnRho::value, 0, COMPONENT_COUNT);
  VectorTools::interpolate(this->dofHandler, initialSlnRho, this->slnPrev);
   
  VectorFunctionFromScalarFunctionObject<DIM> InitialSlnMomentumX(InitialSlnMomentumX::value, 1, COMPONENT_COUNT);
  VectorTools::interpolate(this->dofHandler, InitialSlnMomentumX, this->slnUtil);
  this->slnPrev += this->slnUtil;

  VectorFunctionFromScalarFunctionObject<DIM> InitialSlnMomentumY(InitialSlnMomentumY::value, 2, COMPONENT_COUNT);
  VectorTools::interpolate(this->dofHandler, InitialSlnMomentumY, this->slnUtil);
  this->slnPrev += this->slnUtil;

  VectorFunctionFromScalarFunctionObject<DIM> InitialSlnMomentumZ(InitialSlnMomentumZ::value, 3, COMPONENT_COUNT);
  VectorTools::interpolate(this->dofHandler, InitialSlnMomentumZ, this->slnUtil);
  this->slnPrev += this->slnUtil;
  
  VectorFunctionFromScalarFunctionObject<DIM> InitialSlnB1(InitialSlnB1::value, 4, COMPONENT_COUNT);
  VectorTools::interpolate(this->dofHandler, InitialSlnB1, this->slnUtil);
  this->slnPrev += this->slnUtil;

  VectorFunctionFromScalarFunctionObject<DIM> InitialSlnB2(InitialSlnB2::value, 5, COMPONENT_COUNT);
  VectorTools::interpolate(this->dofHandler, InitialSlnB2, this->slnUtil);
  this->slnPrev += this->slnUtil;

  VectorFunctionFromScalarFunctionObject<DIM> InitialSlnB3(InitialSlnB3::value, 6, COMPONENT_COUNT);
  VectorTools::interpolate(this->dofHandler, InitialSlnB3, this->slnUtil);
  this->slnPrev += this->slnUtil;

  VectorFunctionFromScalarFunctionObject<DIM> initialSlnEnergy(InitialSlnEnergy::value, 7, COMPONENT_COUNT);
  VectorTools::interpolate(this->dofHandler, initialSlnEnergy, this->slnUtil);
  this->slnPrev += this->slnUtil;

  this->slnLin = this->slnPrev;
  d currentTime = 0.;
  for (ui timeStep = 0; currentTime < T_FINAL; timeStep++, currentTime += DELTA_T)
  {
    Timer timer;
    timer.start();
    for(ui l=0;l<8;l++){
      assemble_system(timeStep == 0);
      solve(solution);
      this->slnUtil= solution;
      this->slnUtil-=this->slnLin;
      this->slnLin = solution;
      std::cout << "   ln: " << l << " er: " << this->slnUtil.linfty_norm() << std::endl; // debug only
      if (this->slnUtil.linfty_norm()<1e-10) break;
    }
//     for(ui i=0;i<this->slnPrev.size();i+=COMPONENT_COUNT){ // display arrays
//       for(ui k=0;k<COMPONENT_COUNT;k++)
//         std::cout <<this->slnPrev[i+k]<<":"<<solution[i+k]<<" ";
//       std::cout << std::endl;
//     }
    //break;
    this->slnPrev = this->slnLin;
    outputResults(timeStep, currentTime);
    Eq::currentTime = currentTime;
    timer.stop();
    std::cout << "Time step #" << timeStep << " : " << timer.wall_time() << " s." << std::endl;
  }
}

std::vector<PeriodicBdrInfo> MHDSolver::periodicBdr;

void MHDSolver::JacobiM(double A[3][COMPONENT_COUNT][COMPONENT_COUNT], 
                              Vector<double> lv)
{
  double v[COMPONENT_COUNT], iRh, iRh2, Uk, p, gmmo, Ec1, Ec2, Ec3, E1, E2, E3;

  // using shorter notation for old solution
  // order of the variables is following: rho, v(3), B(3), U, J(3)
  for (unsigned int i = 0; i < COMPONENT_COUNT; i++)
    v[i] = lv[i];

  iRh = 1.0 / v[0];
  iRh2 = iRh*iRh;
  Uk = iRh*(v[1] * v[1] + v[2] * v[2] + v[3] * v[3]);
  p = GAMMA*(v[7] - (v[4] * v[4] + v[5] * v[5] + v[6] * v[6]) - Uk);
  gmmo = GAMMA - 1.0;
  Ec1 = (v[3] * v[5] - v[2] * v[6])*iRh;
  Ec2 = (v[1] * v[6] - v[3] * v[4])*iRh;
  Ec3 = (v[2] * v[4] - v[1] * v[5])*iRh;
  E1 = Ec1 + ETA*v[8];
  E2 = Ec2 + ETA*v[9];
  E3 = Ec3 + ETA*v[10];

  //A[0][0][0] = 0;
  A[0][0][1] = 1;
  //A[0][0][2] = 0;
  //A[0][0][3] = 0;
  //A[0][0][4] = 0;
  //A[0][0][5] = 0;
  //A[0][0][6] = 0;
  //A[0][0][7] = 0;
  //A[0][0][8] = 0;
  //A[0][0][9] = 0;
  //A[0][0][10] = 0;

  A[0][1][0] = -(v[1] * v[1] * iRh2) + gmmo*Uk*.5*iRh;
  A[0][1][1] = (2  - gmmo)*v[1]*iRh;
  A[0][1][2] = -gmmo*v[2]*iRh;
  A[0][1][3] = -gmmo*v[3]*iRh;
  A[0][1][4] = -GAMMA*v[4];
  A[0][1][5] = (1 - gmmo)*v[5];
  A[0][1][6] = (1 - gmmo)*v[6];
  A[0][1][7] = 0.5*gmmo;
  //A[0][1][8] = 0;
  //A[0][1][9] = 0;
  //A[0][1][10] = 0;

  A[0][2][0] = -v[1] * v[2]*iRh2;
  A[0][2][1] = v[2] * iRh;
  A[0][2][2] = v[1] * iRh;
  //A[0][2][3] = 0;
  A[0][2][4] = -v[5];
  A[0][2][5] = -v[4];
  //A[0][2][6] = 0;
  //A[0][2][7] = 0;
  //A[0][2][8] = 0;
  //A[0][2][9] = 0;
  //A[0][2][10] = 0;

  A[0][3][0] = -v[1] * v[3]*iRh2;
  A[0][3][1] = v[3] * iRh;
  //A[0][3][2] = 0;
  A[0][3][3] = v[1] * iRh;
  A[0][3][4] = -v[6];
  //A[0][3][5] = 0;
  A[0][3][6] = -v[4];
  //A[0][3][7] = 0;
  //A[0][3][8] = 0;
  //A[0][3][9] = 0;
  //A[0][3][10] = 0;

  //A[0][4][0] = 0;
  //A[0][4][1] = 0;
  //A[0][4][2] = 0;
  //A[0][4][3] = 0;
  //A[0][4][4] = 0;
  //A[0][4][5] = 0;
  //A[0][4][6] = 0;
  //A[0][4][7] = 0;
  //A[0][4][8] = 0;
  //A[0][4][9] = 0;
  //A[0][4][10] = 0;

  A[0][5][0] = Ec3*iRh;
  A[0][5][1] = v[5] * iRh;
  A[0][5][2] =-v[4] * iRh;
  //A[0][5][3] = 0;
  A[0][5][4] =-v[2] * iRh;
  A[0][5][5] = v[1] * iRh;
  //A[0][5][6] = 0;
  //A[0][5][7] = 0;
  //A[0][5][8] = 0;
  //A[0][5][9] = 0;
  A[0][5][10] = -ETA;

  A[0][6][0] = -Ec2*iRh;
  A[0][6][1] = v[6] * iRh;
  //A[0][6][2] = 0;
  A[0][6][3] = -v[4] * iRh;
  A[0][6][4] = -v[3] * iRh;
  //A[0][6][5] = 0;
  A[0][6][6] = v[1] * iRh;
  //A[0][6][7] = 0;
  //A[0][6][8] = 0;
  A[0][6][9] = ETA;
  //A[0][6][10] = 0;

  A[0][7][0] = 2 * iRh*(v[5]*Ec3 - v[6]*Ec2) + v[1] * gmmo*Uk*iRh2 - v[1]*(Uk + p)*iRh2;
  A[0][7][1] = 2 * (v[5] * v[5] + v[6] * v[6])*iRh - 2*v[1] * gmmo*v[1]*iRh2 + (Uk + p)*iRh;
  A[0][7][2] = -2 * v[4] * v[5] * iRh - v[1] * 2 * gmmo*v[2]*iRh2;
  A[0][7][3] = -2 * v[4] * v[6] * iRh - v[1] * 2 * gmmo*v[3]*iRh2;
  A[0][7][4] = -2 * GAMMA*v[4] * v[1] * iRh + 2 * (-v[5] * v[2] - v[6] * v[3])*iRh;
  A[0][7][5] = -2 * GAMMA*v[5] * v[1] * iRh + 2 * (v[5] * v[1]*iRh - E3);
  A[0][7][6] = -2 * GAMMA*v[6] * v[1] * iRh + 2 * (v[6] * v[1]*iRh + E2);
  A[0][7][7] = GAMMA*v[1]*iRh;
  //A[0][7][8] = 0;
  A[0][7][9] = 2 * ETA*v[6];
  A[0][7][10] = -2 * ETA*v[5];

  //A[0][8][0] = 0;
  //A[0][8][1] = 0;
  //A[0][8][2] = 0;
  //A[0][8][3] = 0;
  //A[0][8][4] = 0;
  //A[0][8][5] = 0;
  //A[0][8][6] = 0;
  //A[0][8][7] = 0;
  //A[0][8][8] = 0;
  //A[0][8][9] = 0;
  //A[0][8][10] = 0;

  //A[0][9][0] = 0;
  //A[0][9][1] = 0;
  //A[0][9][2] = 0;
  //A[0][9][3] = 0;
  //A[0][9][4] = 0;
  //A[0][9][5] = 0;
  A[0][9][6] = 1;
  //A[0][9][7] = 0;
  //A[0][9][8] = 0;
  //A[0][9][9] = 0;
  //A[0][9][10] = 0;

  //A[0][10][0] = 0;
  //A[0][10][1] = 0;
  //A[0][10][2] = 0;
  //A[0][10][3] = 0;
  //A[0][10][4] = 0;
  A[0][10][5] = -1;
  //A[0][10][6] = 0;
  //A[0][10][7] = 0;
  //A[0][10][8] = 0;
  //A[0][10][9] = 0;
  //A[0][10][10] = 0;

  if (DIM>1){
  //A[1][0][0] = 0;
  //A[1][0][1] = 0;
  A[1][0][2] = 1;
  //A[1][0][3] = 0;
  //A[1][0][4] = 0;
  //A[1][0][5] = 0;
  //A[1][0][6] = 0;
  //A[1][0][7] = 0;
  //A[1][0][8] = 0;
  //A[1][0][9] = 0;
  //A[1][0][10] = 0;

  A[1][1][0] = -v[1] * v[2]*iRh2;
  A[1][1][1] = v[2] * iRh;
  A[1][1][2] = v[1] * iRh;
  //A[1][1][3] = 0;
  A[1][1][4] = -v[5];
  A[1][1][5] = -v[4];
  //A[1][1][6] = 0;
  //A[1][1][7] = 0;
  //A[1][1][8] = 0;
  //A[1][1][9] = 0;
  //A[1][1][10] = 0;

  A[1][2][0] = -v[2] * v[2] * iRh2 + gmmo*Uk*.5*iRh;
  A[1][2][1] = -gmmo*v[1]*iRh;
  A[1][2][2] = (2 - gmmo)*v[2]*iRh;
  A[1][2][3] = -gmmo*v[3]*iRh;
  A[1][2][4] = (1 - gmmo)*v[4];
  A[1][2][5] = -GAMMA*v[5];
  A[1][2][6] = (1 - gmmo)*v[6];
  A[1][2][7] = 0.5*gmmo;
  //A[1][2][8] = 0;
  //A[1][2][9] = 0;
  //A[1][2][10] = 0;

  A[1][3][0] = -v[2] * v[3]*iRh2;
  //A[1][3][1] = 0;
  A[1][3][2] = v[3] * iRh;
  A[1][3][3] = v[2] * iRh;
  //A[1][3][4] = 0;
  A[1][3][5] = -v[6];
  A[1][3][6] = -v[5];
  //A[1][3][7] = 0;
  //A[1][3][8] = 0;
  //A[1][3][9] = 0;
  //A[1][3][10] = 0;

  A[1][4][0] = -Ec3*iRh;
  A[1][4][1] =-v[5] * iRh;
  A[1][4][2] = v[4] * iRh;
  //A[1][4][3] = 0;
  A[1][4][4] = v[2] * iRh;
  A[1][4][5] =-v[1] * iRh;
  //A[1][4][6] = 0;
  //A[1][4][7] = 0;
  //A[1][4][8] = 0;
  //A[1][4][9] = 0;
  A[1][4][10] = ETA;

  //A[1][5][0] = 0;
  //A[1][5][1] = 0;
  //A[1][5][2] = 0;
  //A[1][5][3] = 0;
  //A[1][5][4] = 0;
  //A[1][5][5] = 0;
  //A[1][5][6] = 0;
  //A[1][5][7] = 0;
  //A[1][5][8] = 0;
  //A[1][5][9] = 0;
  //A[1][5][10] = 0;

  A[1][6][0] = Ec1*iRh;
  //A[1][6][1] = 0;
  A[1][6][2] = v[6] * iRh;
  A[1][6][3] =-v[5] * iRh;
  //A[1][6][4] = 0;
  A[1][6][5] =-v[3] * iRh;
  A[1][6][6] = v[2] * iRh;
  //A[1][6][7] = 0;
  A[1][6][8] = -ETA;
  //A[1][6][9] = 0;
  //A[1][6][10] = 0;

  A[1][7][0] = 2 * iRh*(-v[4]*Ec3 + v[6]*Ec1) + v[2] * gmmo*Uk*iRh2 - v[2]*(Uk + p)*iRh2;
  A[1][7][1] = -2 * v[4] * v[5] * iRh - 2 * gmmo*v[1] * v[2] * iRh2;
  A[1][7][2] = 2 * (v[4] * v[4] + v[6] * v[6])*iRh - 2*v[2] * gmmo*v[2]*iRh2 + (Uk + p)*iRh;
  A[1][7][3] = -2 * v[5] * v[6] * iRh - 2 * gmmo*v[2] * v[3] * iRh2;
  A[1][7][4] = -2 * GAMMA*v[4] * v[2] * iRh + 2 * (v[4] * v[2]*iRh + E3);
  A[1][7][5] = -2 * GAMMA*v[5] * v[2] * iRh + 2 * (-v[4] * v[1] - v[6] * v[3])*iRh;
  A[1][7][6] = -2 * GAMMA*v[6] * v[2] * iRh + 2 * (v[6] * v[2]*iRh - E1);
  A[1][7][7] = GAMMA*v[2]*iRh;
  A[1][7][8] = -2 * ETA*v[6];
  //A[1][7][9] = 0;
  A[1][7][10] = 2 * ETA*v[4];

  //A[1][8][0] 0;
  //A[1][8][1] 0;
  //A[1][8][2] 0;
  //A[1][8][3] 0;
  //A[1][8][4] 0;
  //A[1][8][5] 0;
  A[1][8][6] = -1;
  //A[1][8][7] 0;
  //A[1][8][8] 0;
  //A[1][8][9] 0;
  //A[1][8][10] 0;

  //A[1][9][0] = 0;
  //A[1][9][1] = 0;
  //A[1][9][2] = 0;
  //A[1][9][3] = 0;
  //A[1][9][4] = 0;
  //A[1][9][5] = 0;
  //A[1][9][6] = 0;
  //A[1][9][7] = 0;
  //A[1][9][8] = 0;
  //A[1][9][9] = 0;
  //A[1][9][10] = 0;

  //A[1][10][0] - 0;
  //A[1][10][1] - 0;
  //A[1][10][2] - 0;
  //A[1][10][3] - 0;
  A[1][10][4] = 1;
  //A[1][10][5] - 0;
  //A[1][10][6] - 0;
  //A[1][10][7] - 0;
  //A[1][10][8] - 0;
  //A[1][10][9] - 0;
  //A[1][10][10] - 0;
  }

  if (DIM>2){
  //A[2][0][0] = 0;
  //A[2][0][1] = 0;
  //A[2][0][2] = 0;
  A[2][0][3] = 1;
  //A[2][0][4] = 0;
  //A[2][0][5] = 0;
  //A[2][0][6] = 0;
  //A[2][0][7] = 0;
  //A[2][0][8] = 0;
  //A[2][0][9] = 0;
  //A[2][0][10] = 0;

  A[2][1][0] = -v[1] * v[3]*iRh2;
  A[2][1][1] = v[3] * iRh;
  //A[2][1][2] = 0;
  A[2][1][3] = v[1] * iRh;
  A[2][1][4] = -v[6];
  //A[2][1][5] = 0;
  A[2][1][6] = -v[4];
  //A[2][1][7] = 0;
  //A[2][1][8] = 0;
  //A[2][1][9] = 0;
  //A[2][1][10] = 0;

  A[2][2][0] = -v[2] * v[3]*iRh2;
  //A[2][2][1] = 0;
  A[2][2][2] = v[3] * iRh;
  A[2][2][3] = v[2] * iRh;
  //A[2][2][4] = 0;
  A[2][2][5] = -v[6];
  A[2][2][6] = -v[5];
  //A[2][2][7] = 0;
  //A[2][2][8] = 0;
  //A[2][2][9] = 0;
  //A[2][2][10] = 0;

  A[2][3][0] = -v[3] * v[3] * iRh2 + (gmmo*Uk)*.5*iRh;
  A[2][3][1] = -gmmo*v[1]*iRh;
  A[2][3][2] = -gmmo*v[2]*iRh;
  A[2][3][3] = (2 - gmmo)*v[3]*iRh;
  A[2][3][4] = (1 - gmmo)*v[4];
  A[2][3][5] = (1 - gmmo)*v[5];
  A[2][3][6] = -GAMMA*v[6];
  A[2][3][7] = 0.5*gmmo;
  //A[2][3][8] = 0;
  //A[2][3][9] = 0;
  //A[2][3][10] = 0;

  A[2][4][0] = Ec2*iRh;
  A[2][4][1] =-v[6] * iRh;
  //A[2][4][2] = 0;
  A[2][4][3] = v[4] * iRh;
  A[2][4][4] = v[3] * iRh;
  //A[2][4][5] = 0;
  A[2][4][6] =-v[1] * iRh;
  //A[2][4][7] = 0;
  //A[2][4][8] = 0;
  A[2][4][9] = -ETA;
  //A[2][4][10] = 0;

  A[2][5][0] = -Ec1*iRh;
  //A[2][5][1] = 0;
  A[2][5][2] =-v[6] * iRh;
  A[2][5][3] = v[5] * iRh;
  //A[2][5][4] = 0;
  A[2][5][5] = v[3] * iRh;
  A[2][5][6] =-v[2] * iRh;
  //A[2][5][7] = 0;
  A[2][5][8] = ETA;
  //A[2][5][9] = 0;
  //A[2][5][10] = 0;

  //A[2][6][0] = 0;
  //A[2][6][1] = 0;
  //A[2][6][2] = 0;
  //A[2][6][3] = 0;
  //A[2][6][4] = 0;
  //A[2][6][5] = 0;
  //A[2][6][6] = 0;
  //A[2][6][7] = 0;
  //A[2][6][8] = 0;
  //A[2][6][9] = 0;
  //A[2][6][10] = 0;

  A[2][7][0] = 2 * iRh*(v[4]*Ec2 - v[5]*Ec1) + v[3] * gmmo*Uk*iRh2 - v[3]*(Uk + p)*iRh2;
  A[2][7][1] = -2*v[4]*v[6]*iRh - 2*gmmo*v[1]*v[3]*iRh2;
  A[2][7][2] = -2*v[5]*v[6]*iRh - 2*gmmo*v[2]*v[3]*iRh2;
  A[2][7][3] = 2 * (v[4] * v[4] + v[5] * v[5])*iRh + 2*v[3] * gmmo*v[3]*iRh2 + (Uk + p)*iRh;
  A[2][7][4] = -2 * GAMMA*v[4] * v[3] * iRh + 2 * (v[4] * v[3]*iRh - E2);
  A[2][7][5] = -2 * GAMMA*v[5] * v[3] * iRh + 2 * (v[5] * v[3]*iRh + E1);
  A[2][7][6] = -2 * GAMMA*v[6] * v[3] * iRh + 2 * (-v[4] * v[1] - v[5] * v[2])*iRh;
  A[2][7][7] = GAMMA*v[3]*iRh;
  A[2][7][8] = 2 * ETA*v[5];
  A[2][7][9] = -2 * ETA*v[4];
  //A[2][7][10] = 0;

  //A[2][8][0] = 0;
  //A[2][8][1] = 0;
  //A[2][8][2] = 0;
  //A[2][8][3] = 0;
  //A[2][8][4] = 0;
  A[2][8][5] = 1;
  //A[2][8][6] = 0;
  //A[2][8][7] = 0;
  //A[2][8][8] = 0;
  //A[2][8][9] = 0;
  //A[2][8][10] = 0;

  //A[2][9][0] = 0;
  //A[2][9][1] = 0;
  //A[2][9][2] = 0;
  //A[2][9][3] = 0;
  A[2][9][4] = -1;
  //A[2][9][5] = 0;
  //A[2][9][6] = 0;
  //A[2][9][7] = 0;
  //A[2][9][8] = 0;
  //A[2][9][9] = 0;
  //A[2][9][10] = 0;

  //A[2][10][0] = 0;
  //A[2][10][1] = 0;
  //A[2][10][2] = 0;
  //A[2][10][3] = 0;
  //A[2][10][4] = 0;
  //A[2][10][5] = 0;
  //A[2][10][6] = 0;
  //A[2][10][7] = 0;
  //A[2][10][8] = 0;
  //A[2][10][9] = 0;
  //A[2][10][10] = 0;
    }
}

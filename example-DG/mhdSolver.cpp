#include "common.h"

typedef EquationImplementation Eq;

Vector<d> MHDSolver::slnPrev;
NumFlux* numFlux;
DirichletBoundaryCondition bc;

MHDSolver::MHDSolver()
    :
    feSystem(dealii::FE_DGQ<COMPONENT_COUNT, DIM>(DG_ORDER), COMPONENT_COUNT),
    dofHandler(triangulation),
    mapping(),
    quad(2 * DG_ORDER),
    quadFace(2 * DG_ORDER)
{
    numFlux = new NumFluxUpwind();
}

void MHDSolver::setup_system()
{
    dofHandler.distribute_dofs(feSystem);

    dealii::DoFRenumbering::component_wise(dofHandler);

    DynamicSparsityPattern dsp(dofHandler.n_dofs());
    DoFTools::make_sparsity_pattern(dofHandler, dsp);
    sparsityPattern.copy_from(dsp);

    systemMatrix.reinit(sparsityPattern);
    rightHandSide.reinit(dofHandler.n_dofs());
    solution.reinit(dofHandler.n_dofs());
    slnPrev.reinit(dofHandler.n_dofs());
}

void MHDSolver::assemble_system()
{
    MeshWorker::IntegrationInfoBox<COMPONENT_COUNT, DIM> info_box;

    // \todo This is wrong probably.
    info_box.initialize_gauss_quadrature(3 * DG_ORDER, 3 * DG_ORDER, 3 * DG_ORDER);

    AnyData solution_data;
    solution_data.add(&slnPrev, "solution");
    info_box.cell_selector.add("solution", true, false, false);
    info_box.boundary_selector.add("solution", true, false, false);
    info_box.face_selector.add("solution", true, false, false);

    info_box.initialize_update_flags();
    UpdateFlags update_flags = update_quadrature_points | update_values | update_gradients;
    info_box.add_update_flags(update_flags, true, true, true, true);

    // \todo What about multiple FEs in feCollection?
    info_box.initialize(feSystem, mapping, solution_data, rightHandSide);

    // \todo This has to be done properly for hpDoFHandler (varying number of DOFs per cell)
    MeshWorker::DoFInfo<COMPONENT_COUNT, DIM> dof_info(dofHandler);

    MeshWorker::Assembler::SystemSimple < SparseMatrix<d>, Vector<d> > assembler;
    assembler.initialize(systemMatrix, rightHandSide);

    // \todo This comes from tutorial, it may need some adjustment.
    MeshWorker::loop<COMPONENT_COUNT, DIM, MeshWorker::DoFInfo<COMPONENT_COUNT, DIM>, MeshWorker::IntegrationInfoBox<COMPONENT_COUNT, DIM> >
        (dofHandler.begin_active(), dofHandler.end(), dof_info, info_box, &MHDSolver::assembleVolumetric,
        &MHDSolver::assembleBoundaryEdge, &MHDSolver::assembleInternalEdge, assembler);
}

void MHDSolver::assembleVolumetric(DoFInfo &dinfo,
    CellInfo &info)
{
    const FEValuesBase<COMPONENT_COUNT, DIM> &fe_v = info.fe_values();
    FullMatrix<d> &local_matrix = dinfo.matrix(0).matrix;
    Vector<d> &local_vector = dinfo.vector(0).block(0);
    const std::vector<d> &JxW = fe_v.get_JxW_values();

    const ui dofs_per_cell = info.finite_element().dofs_per_cell;

    std::vector<dealii::Vector<double> > prev_values(fe_v.n_quadrature_points, dealii::Vector<double>(COMPONENT_COUNT));
    for (ui i = 0; i < COMPONENT_COUNT; i++)
        for (ui point = 0; point < fe_v.n_quadrature_points; ++point)
            prev_values[point][i] = info.values[0][i][point];

    // Components
    std::vector<int> components(dofs_per_cell);
    for (ui i = 0; i < dofs_per_cell; ++i)
        components[i] = info.finite_element().system_to_component_index(i).first;

    for (ui point = 0; point < fe_v.n_quadrature_points; ++point)
    {
        for (ui i = 0; i < fe_v.dofs_per_cell; ++i)
        {
            for (ui j = 0; j < fe_v.dofs_per_cell; ++j)
            {
                local_matrix(i, j) += JxW[point] * Eq::matrixVolValue(components[j], components[i],
                    fe_v.shape_value(j, point), fe_v.shape_value(i, point),
                    prev_values[point], fe_v.shape_grad(j, point), fe_v.shape_grad(i, point),
                    vecDimVec(), fe_v.quadrature_point(point));
            }
            local_vector(i) += JxW[point] * Eq::rhsVolValue(components[i],
                fe_v.shape_value(i, point),
                prev_values[point], fe_v.shape_grad(i, point),
                vecDimVec(), fe_v.quadrature_point(point));
        }
    }
}

void MHDSolver::assembleBoundaryEdge(DoFInfo &dinfo,
    CellInfo &info)
{
    const FEValuesBase<COMPONENT_COUNT, DIM> &fe_v = info.fe_values();
    FullMatrix<d> &local_matrix = dinfo.matrix(0).matrix;
    Vector<d> &local_vector = dinfo.vector(0).block(0);

    const std::vector<d> &JxW = fe_v.get_JxW_values();
    const std::vector<Point<DIM> > &normals = fe_v.get_normal_vectors();

    const ui dofs_per_cell = info.finite_element().dofs_per_cell;

    // Previous values.
    std::vector<dealii::Vector<double> > prev_values(fe_v.n_quadrature_points, dealii::Vector<double>(COMPONENT_COUNT));
    for (ui i = 0; i < COMPONENT_COUNT; i++)
        for (ui point = 0; point < fe_v.n_quadrature_points; ++point)
            prev_values[point][i] = info.values[0][i][point];

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
                    vecDimVec(), vec(), fe_v.quadrature_point(point), normals[point], numFlux);
            }
            local_vector(i) += JxW[point] * Eq::rhsBoundaryEdgeValue(components[i],
                fe_v.shape_value(i, point),
                prev_values[point], fe_v.shape_grad(i, point),
                vecDimVec(), vec(), fe_v.quadrature_point(point), normals[point], numFlux, &bc);
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
    const FEValuesBase<COMPONENT_COUNT, DIM> &fe_v = info1.fe_values();

    // For additional shape functions, we have to ask the neighbors
    // FEValuesBase.
    const FEValuesBase<COMPONENT_COUNT, DIM> &fe_v_neighbor = info2.fe_values();

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

    // Previous values.
    std::vector<dealii::Vector<double> > prev_values1(fe_v.n_quadrature_points, dealii::Vector<double>(COMPONENT_COUNT));
    for (ui i = 0; i < COMPONENT_COUNT; i++)
        for (ui point = 0; point < fe_v.n_quadrature_points; ++point)
            prev_values1[point][i] = info1.values[0][i][point];

    std::vector<dealii::Vector<double> > prev_values2(fe_v.n_quadrature_points, dealii::Vector<double>(COMPONENT_COUNT));
    for (ui i = 0; i < COMPONENT_COUNT; i++)
        for (ui point = 0; point < fe_v_neighbor.n_quadrature_points; ++point)
            prev_values2[point][i] = info2.values[0][i][point];

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
                    fe_v.shape_value(j, point), fe_v.shape_value(i, point),
                    prev_values1[point], prev_values2[point], fe_v.shape_grad(j, point), fe_v.shape_grad(i, point),
                    vecDimVec(), vecDimVec(), false, false, fe_v.quadrature_point(point), normals[point], numFlux);
            }
        }

        for (ui k = 0; k < fe_v_neighbor.dofs_per_cell; ++k)
        {
            for (ui j = 0; j < fe_v.dofs_per_cell; ++j)
            {
                u1_v2_matrix(k, j) += JxW[point] * Eq::matrixInternalEdgeValue(components1[j], components1[k],
                    fe_v.shape_value(j, point), fe_v_neighbor.shape_value(k, point),
                    prev_values1[point], prev_values2[point], fe_v.shape_grad(j, point), fe_v_neighbor.shape_grad(k, point),
                    vecDimVec(), vecDimVec(), false, true, fe_v.quadrature_point(point), normals[point], numFlux);
            }
        }

        for (ui i = 0; i < fe_v.dofs_per_cell; ++i)
        {
            for (ui l = 0; l < fe_v_neighbor.dofs_per_cell; ++l)
            {
                u2_v1_matrix(i, l) += JxW[point] * Eq::matrixInternalEdgeValue(components1[l], components1[i],
                    fe_v_neighbor.shape_value(l, point), fe_v.shape_value(i, point),
                    prev_values1[point], prev_values2[point], fe_v_neighbor.shape_grad(l, point), fe_v.shape_grad(i, point),
                    vecDimVec(), vecDimVec(), true, false, fe_v.quadrature_point(point), normals[point], numFlux);
            }
        }

        for (ui k = 0; k < fe_v_neighbor.dofs_per_cell; ++k)
        {
            for (ui l = 0; l < fe_v_neighbor.dofs_per_cell; ++l)
            {
                u2_v2_matrix(k, l) += JxW[point] * Eq::matrixInternalEdgeValue(components1[l], components1[k],
                    fe_v_neighbor.shape_value(l, point), fe_v_neighbor.shape_value(k, point),
                    prev_values1[point], prev_values2[point], fe_v_neighbor.shape_grad(l, point), fe_v_neighbor.shape_grad(k, point),
                    vecDimVec(), vecDimVec(), true, true, fe_v.quadrature_point(point), normals[point], numFlux);
            }
        }
    }
}

void MHDSolver::solve(Vector<d> &solution)
{
    if (PRINT_ALGEBRA)
    {
        std::cout << "  Printing system... " << std::endl;

        std::string matrix_file = "MagMatrix_";
        std::string rhs_file = "MagRhs_";

        std::ofstream matrix_out(matrix_file);
        std::ofstream rhs_out(rhs_file);

        systemMatrix.print(matrix_out);
        rightHandSide.print(rhs_out, 3, true, false);

        matrix_out.close();
        rhs_out.close();
    }
    SolverControl           solver_control(1000, 1e-12);
    SolverRichardson<>      solver(solver_control);

    // Here we create the preconditioner,
    PreconditionBlockSSOR<SparseMatrix<d> > preconditioner;

    // then assign the matrix to it and set the right block size:
    preconditioner.initialize(systemMatrix, feSystem.dofs_per_cell);

    // After these preparations we are ready to start the linear solver.
    solver.solve(systemMatrix, solution, rightHandSide, preconditioner);
}

void MHDSolver::outputResults(ui timeStep, d currentTime) const
{
    Postprocessor postprocessor;
    DataOut<COMPONENT_COUNT, DoFHandler<COMPONENT_COUNT, DIM> > data_out;
    data_out.attach_dof_handler(dofHandler);
    const DataOut<COMPONENT_COUNT, DoFHandler<COMPONENT_COUNT, DIM> >::DataVectorType data_vector_type = DataOut<COMPONENT_COUNT, DoFHandler<COMPONENT_COUNT, DIM> >::type_dof_data;
    data_out.add_data_vector(slnPrev, postprocessor);
    data_out.build_patches(mapping);
    std::stringstream ss;
    ss << "solution-";
    ss << timeStep;
    ss << ".vtk";
    std::ofstream output(ss.str());
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

    // Initial sln.
    VectorFunctionFromScalarFunctionObject<DIM> initialSln(InitialSln::value, 0, COMPONENT_COUNT);
    VectorTools::interpolate(this->dofHandler, initialSln, this->slnPrev);

    d currentTime = 0.;
    for (ui timeStep = 0; currentTime < T_FINAL; timeStep++, currentTime += DELTA_T)
    {
        Timer timer;
        timer.start();
        assemble_system();
        solve(solution);
        this->slnPrev = solution;
        outputResults(timeStep, currentTime);
        timer.stop();
        std::cout << "Time step #" << timeStep << " : " << timer.wall_time() << " s.";
    }

}
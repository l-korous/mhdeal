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
  void outputResults(ui timeStep, d currentTime) const;
  void add_markers(Triangulation<COMPONENT_COUNT, DIM>::cell_iterator cell);

  Triangulation<COMPONENT_COUNT, DIM>   triangulation;
  FESystem<COMPONENT_COUNT, DIM> feSystem;
  const MappingQ1<COMPONENT_COUNT, DIM> mapping;
  const QGauss<DIM> quad;
  const QGauss<DIM - 1> quadFace;
  DoFHandler<COMPONENT_COUNT, DIM>      dofHandler;
  ConstraintMatrix     hangingNodeConstraints;

  SparsityPattern      sparsityPattern;
  SparseMatrix<d> systemMatrix;

  Vector<d>       solution;
  static Vector<d>       slnPrev;
  Vector<d>       rightHandSide;

  typedef MeshWorker::DoFInfo<COMPONENT_COUNT, DIM> DoFInfo;
  typedef MeshWorker::IntegrationInfo<COMPONENT_COUNT, DIM> CellInfo;

  static void assembleVolumetric(DoFInfo &dinfo, CellInfo &info);
  static void assembleBoundaryEdge(DoFInfo &dinfo, CellInfo &info);
  static void assembleInternalEdge(DoFInfo &dinfo1, DoFInfo &dinfo2, CellInfo &info1, CellInfo &info2);
};
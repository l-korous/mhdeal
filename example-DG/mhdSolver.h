#ifndef MHD_SOLVER_H
#define MHD_SOLVER_H

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
  void add_markers(Triangulation<DIM>::cell_iterator cell);

  Triangulation<DIM>   triangulation;
  FESystem<DIM> feSystem;
  const MappingQ1<DIM> mapping;
  const QGauss<DIM> quad;
  const QGauss<DIM - 1> quadFace;
  DoFHandler<DIM>      dofHandler;
  ConstraintMatrix     hangingNodeConstraints;

  SparsityPattern      sparsityPattern;
  SparseMatrix<d> systemMatrix;

  Vector<d>       solution;
  static Vector<d>       slnPrev;
  Vector<d>       rightHandSide;

  typedef MeshWorker::DoFInfo<DIM> DoFInfo;
  typedef MeshWorker::IntegrationInfo<DIM> CellInfo;

  static void assembleVolumetric(DoFInfo &dinfo, CellInfo &info);
  static void assembleBoundaryEdge(DoFInfo &dinfo, CellInfo &info);
  static void assembleInternalEdge(DoFInfo &dinfo1, DoFInfo &dinfo2, CellInfo &info1, CellInfo &info2);
};

#endif
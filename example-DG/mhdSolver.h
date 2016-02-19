#ifndef MHD_SOLVER_H
#define MHD_SOLVER_H

#include <deal.II/fe/fe_system.h>
#include <deal.II/meshworker/simple.h>
#include <deal.II/meshworker/loop.h>


#include "definitions.h"


typedef dealii::MeshWorker::DoFInfo<DIM> DoFInfo;
typedef dealii::MeshWorker::IntegrationInfo<DIM> CellInfo;

class PeriodicBdrInfo
{
public:
  dealii::Point<DIM> center;
  std::vector<dealii::Vector<double> > prev_values;

  PeriodicBdrInfo(dealii::Point<DIM> center, std::vector<dealii::Vector<double> > prev_values) :
    center(center), prev_values(prev_values)
  {}

};

class MHDSolver
{
public:
  MHDSolver();
  void run();

private:
  void setup_system();
  void assemble_system(bool firstIteration);
  void solve(dealii::Vector<d> &solution);
  void solveOneStep(dealii::Vector<d> &solution);
  void outputResults(ui timeStep, d currentTime, int linStep = -1) const;
  void add_markers(dealii::Triangulation<DIM>::cell_iterator cell);
  static void JacobiM(double A[3][COMPONENT_COUNT][COMPONENT_COUNT], 
               dealii::Vector<double> lv);

  dealii::Triangulation<DIM>   triangulation;
  dealii::FESystem<DIM> feSystem;
  const dealii::MappingQ1<DIM> mapping;
  const dealii::QGauss<DIM> quad;
  const dealii::QGauss<DIM - 1> quadFace;
  dealii::DoFHandler<DIM>      dofHandler;
  dealii::ConstraintMatrix     hangingNodeConstraints;

  dealii::SparsityPattern      sparsityPattern;
  dealii::SparseMatrix<d> systemMatrix;

  dealii::Vector<d>       solution;
  static dealii::Vector<d> slnPrev;
  static dealii::Vector<d> slnLin;
  static dealii::Vector<d> slnUtil;

  dealii::Vector<d>       rightHandSide;
  static d A[3][COMPONENT_COUNT][COMPONENT_COUNT];         // Jacobi matrixes of the fluxes

  void JacobiM(double A[][COMPONENT_COUNT][COMPONENT_COUNT], std::vector<dealii::Vector<double> > lv, const unsigned int qp);

  static std::vector<PeriodicBdrInfo> periodicBdr;

  // todo: this is now implemented assuming that periodic BC is on left and right boundary only
  static std::vector<dealii::Vector<double> > findCorrespondingInfo(dealii::Point<DIM> myCenter);

  static void assembleVolumetric(DoFInfo &dinfo, CellInfo &info);
  static void assembleVolumetricEmpty(DoFInfo &dinfo, CellInfo &info) {}
  static void saveInfoBoundaryEdge(DoFInfo &dinfo, CellInfo &info);
  static void assembleBoundaryEdge(DoFInfo &dinfo, CellInfo &info);
  static void assembleInternalEdge(DoFInfo &dinfo1, DoFInfo &dinfo2, CellInfo &info1, CellInfo &info2);
  static void assembleInternalEdgeEmpty(DoFInfo &dinfo1, DoFInfo &dinfo2, CellInfo &info1, CellInfo &info2) {}
};

#endif

#ifndef EQUATIONS_H
#define EQUATIONS_H

class DirichletBoundaryCondition;
class NumFlux;
namespace dealii
{
template <typename number> class FullMatrix;
}

#include "definitions.h"

class EquationImplementation
{
public:

  static d matrixVolValue(ui component_i, ui component_j,
    d u_val, d v_val, vec Un_val, dimVec u_grad, dimVec v_grad,
    vecDimVec Un_grad, dealii::Point<DIM> quadPoint);

  static d matrixBoundaryEdgeValue(ui component_i, ui component_j,
    d u_val, d v_val, vec Un_val, dimVec u_grad, dimVec v_grad,
    vecDimVec Un_grad, vec U_bnd_val, dealii::Point<DIM> quadPoint, dealii::Point<DIM> normal, NumFlux* num_flux, DirichletBoundaryCondition* bc, dealii::types::boundary_id);

  static d matrixInternalEdgeValue(ui component_i, ui component_j,
    d u_val, d v_val, vec Un_val, vec Un_valN, dimVec u_grad, dimVec v_grad,
    vecDimVec Un_grad, vecDimVec Un_gradN, bool u_N, bool v_N,
    dealii::Point<DIM> quadPoint, dealii::Point<DIM> normal, NumFlux* num_flux);


  static d rhsVolValue(ui component_i,
    d v_val, vec Un_val, dimVec v_grad,
    vecDimVec Un_grad, dealii::Point<DIM> quadPoint);

  static d rhsBoundaryEdgeValue(ui component_i,
    d v_val, vec Un_val, dimVec v_grad,
    vecDimVec Un_grad, vec U_bnd_val, dealii::Point<DIM> quadPoint, dealii::Point<DIM> normal, NumFlux* num_flux, DirichletBoundaryCondition* bc, dealii::types::boundary_id);

  static void Jacobians(dealii::FullMatrix<double> *J,
    std::vector<dealii::Vector<double> > lv,
    const unsigned int qp);

  static d rhsInternalEdgeValue(ui comp_i,
      d v_val, dimVec v_grad, bool v_N, vec Un_val,
      vecDimVec Un_grad, vec Un_valN, vecDimVec Un_gradN, dealii::Point<DIM> quadPoint, dealii::Point<DIM> normal, NumFlux* num_flux);

  static d currentTime;
};

#endif

#include "common.h"


d EquationImplementation::matrixVolValue(ui comp_i, ui comp_j,
  d u_val, d v_val, vec Un_val, dimVec u_grad, dimVec v_grad,
  vecDimVec Un_grad, Point<DIM> quadPoint)
{
  d result = 0.;
  // Time derivative.
  if (comp_i == comp_j)
    result += u_val * v_val;

  return result;
}


d EquationImplementation::matrixBoundaryEdgeValue(ui comp_i, ui comp_j,
  d u_val, d v_val, vec Un_val, dimVec u_grad, dimVec v_grad,
  vecDimVec Un_grad, vec U_bnd_val, Point<DIM> quadPoint, Point<DIM> normal)
{
  d result = 0.;

  return result;
}


d EquationImplementation::matrixInternalEdgeValue(ui comp_i, ui comp_j,
  d u_val, d v_val, vec Un_val, vec Un_valN, dimVec u_grad, dimVec v_grad,
  vecDimVec Un_grad, vecDimVec Un_gradN, bool u_N, bool v_N, 
  Point<DIM> quadPoint, Point<DIM> normal)
{
  d result = 0.;

  return result;
}



d EquationImplementation::rhsVolValue(ui comp_i,
  d v_val, vec Un_val, dimVec v_grad,
  vecDimVec Un_grad, Point<DIM> quadPoint)
{
  d result = 0.;

  // Time derivative.
  result += Un_val[comp_i] * v_val;

  return result;
}


d EquationImplementation::rhsBoundaryEdgeValue(ui comp_i,
  d v_val, vec Un_val, dimVec v_grad,
  vecDimVec Un_grad, vec U_bnd_val, Point<DIM> quadPoint, Point<DIM> normal)
{
  d result = 0.;

  return result;
}


d EquationImplementation::rhsInternalEdgeValue(ui comp_i,
  d v_val, vec Un_val, dimVec v_grad,
  vecDimVec Un_grad,
  d v_valN, vec Un_valN, dimVec v_gradN,
  vecDimVec Un_gradN, bool v_N, Point<DIM> quadPoint, Point<DIM> normal)
{
  d result = 0.;

  return result;
}

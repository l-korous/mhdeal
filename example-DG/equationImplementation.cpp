#include "common.h"


d EquationImplementation::matrixVolValue(ui comp_i, ui comp_j,
  d u_val, d v_val, vec Un_val, dimVec u_grad, dimVec v_grad,
  vecDimVec Un_grad, Point<DIM> quadPoint)
{
  if (comp_i != comp_j)
    return 0.;

  Point<DIM> beta;
  beta(0) = -quadPoint(1);
  beta(1) = quadPoint(0);
  beta /= beta.norm();

  return -beta * v_grad * u_val;
}


d EquationImplementation::matrixBoundaryEdgeValue(ui comp_i, ui comp_j,
  d u_val, d v_val, vec Un_val, dimVec u_grad, dimVec v_grad,
  vecDimVec Un_grad, vec U_bnd_val, Point<DIM> quadPoint, Point<DIM> normal)
{
  if (comp_i != comp_j)
    return 0.;

  Point<DIM> beta;
  beta(0) = -quadPoint(1);
  beta(1) = quadPoint(0);
  beta /= beta.norm();

  const d beta_n = beta * normal;
  if (beta_n > 0)
    return beta_n * u_val * v_val;
  else
    return 0.;
}


d EquationImplementation::matrixInternalEdgeValue(ui comp_i, ui comp_j,
  d u_val, d v_val, vec Un_val, vec Un_valN, dimVec u_grad, dimVec v_grad,
  vecDimVec Un_grad, vecDimVec Un_gradN, bool u_N, bool v_N, 
  Point<DIM> quadPoint, Point<DIM> normal)
{
  if (comp_i != comp_j)
    return 0.;

  Point<DIM> beta;
  beta(0) = -quadPoint(1);
  beta(1) = quadPoint(0);
  beta /= beta.norm();

  const d beta_n = beta * normal;
  if (beta_n > 0)
  {
    if (!u_N && !v_N)
      return beta_n * u_val * v_val;
    else if (!u_N && v_N)
      return -beta_n * u_val * v_val;
  }
  else
  {
    if (u_N && v_N)
      return -beta_n * u_val * v_val;
    else if (u_N && !v_N)
      return beta_n * u_val * v_val;
  }

  return 0.;
}



d EquationImplementation::rhsVolValue(ui comp_i,
  d v_val, vec Un_val, dimVec v_grad,
  vecDimVec Un_grad, Point<DIM> quadPoint)
{
  return 0.0;
}


d EquationImplementation::rhsBoundaryEdgeValue(ui comp_i,
  d v_val, vec Un_val, dimVec v_grad,
  vecDimVec Un_grad, vec U_bnd_val, Point<DIM> quadPoint, Point<DIM> normal)
{
  Point<DIM> beta;
  beta(0) = -quadPoint(1);
  beta(1) = quadPoint(0);
  beta /= beta.norm();

  const d beta_n = beta * normal;
  if (beta_n <= 0)
    return -beta_n * v_val;
  else return 0.;
}


d EquationImplementation::rhsInternalEdgeValue(ui comp_i,
  d v_val, vec Un_val, dimVec v_grad,
  vecDimVec Un_grad,
  d v_valN, vec Un_valN, dimVec v_gradN,
  vecDimVec Un_gradN, bool v_N, Point<DIM> quadPoint, Point<DIM> normal)
{
  return 0.0;
}

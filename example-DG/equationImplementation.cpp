#include "common.h"

d calculate_flux(double x, double y, double vx, double vy)
{
    double norm = std::max<double>(1e-12, std::sqrt(std::pow(x, 2.) + std::pow(y, 2.)));
    return -y / norm*vx + x / norm*vy;
}

d EquationImplementation::matrixVolValue(ui comp_i, ui comp_j,
  d u_val, d v_val, vec Un_val, dimVec u_grad, dimVec v_grad,
  vecDimVec Un_grad, Point<DIM> quadPoint)
{
  d result = 0.;
  d flux;

  // Time derivative.
  if (comp_i == comp_j)
  {
      result += u_val * v_val;
      flux = calculate_flux(quadPoint(0), quadPoint(1), v_grad[0], v_grad[1]);
      result -= u_val * flux;
  }

  return result;
}


d EquationImplementation::matrixBoundaryEdgeValue(ui comp_i, ui comp_j,
  d u_val, d v_val, vec Un_val, dimVec u_grad, dimVec v_grad,
  vecDimVec Un_grad, vec U_bnd_val, Point<DIM> quadPoint, Point<DIM> normal, NumFlux* num_flux)
{
  d result = 0.;

  d a_dot_n = calculate_flux(quadPoint(0), quadPoint(1), normal(0), normal(1));
  result = a_dot_n * v_val;

  return result;
}


d EquationImplementation::matrixInternalEdgeValue(ui comp_i, ui comp_j,
  d u_val, d v_val, vec Un_val, vec Un_valN, dimVec u_grad, dimVec v_grad,
  vecDimVec Un_grad, vecDimVec Un_gradN, bool u_N, bool v_N, 
  Point<DIM> quadPoint, Point<DIM> normal, NumFlux* num_flux)
{
  d result = 0.;

  d jump_v = v_N ? -v_val : v_val;

  vec numFlux(COMPONENT_COUNT);

  num_flux->calculate(Un_val, Un_valN, quadPoint, normal, numFlux);

  result = numFlux[comp_i] * jump_v;

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
  vecDimVec Un_grad, vec U_bnd_val, Point<DIM> quadPoint, Point<DIM> normal, NumFlux* num_flux, DirichletBoundaryCondition* bc)
{
  d result = 0.;

  d a_dot_n = calculate_flux(quadPoint(0), quadPoint(1), normal(0), normal(1));

  vec bc_state(COMPONENT_COUNT);
  for (ui i = 0; i < COMPONENT_COUNT; i++)
      bc_state = bc->calculate(i, quadPoint);

  vec numFlux(COMPONENT_COUNT);

  num_flux->calculate(Un_val, bc_state, quadPoint, normal, numFlux);

  result = -numFlux[0] * v_val;

  return result;
}


d EquationImplementation::rhsInternalEdgeValue(ui comp_i,
  d v_val, vec Un_val, dimVec v_grad,
  vecDimVec Un_grad,
  d v_valN, vec Un_valN, dimVec v_gradN,
  vecDimVec Un_gradN, bool v_N, Point<DIM> quadPoint, Point<DIM> normal, NumFlux* num_flux)
{
  d result = 0.;

  return result;
}

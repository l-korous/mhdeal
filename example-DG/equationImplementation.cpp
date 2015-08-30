#include "common.h"

d calculate_flux(double x, double y, double nx, double ny)
{
    return FLUX;
}

d EquationImplementation::matrixVolValue(ui comp_i, ui comp_j,
    d u_val, d v_val, vec Un_val, dimVec u_grad, dimVec v_grad,
    vecDimVec Un_grad, Point<DIM> quadPoint)
{
    d result = 0.;

    // Time derivative.
    if (comp_i == comp_j)
    {
        result += u_val * v_val / DELTA_T;
        result -= calculate_flux(quadPoint(0), quadPoint(1), v_grad[0], v_grad[1]) * u_val;
    }

    return result;
}


d EquationImplementation::matrixBoundaryEdgeValue(ui comp_i, ui comp_j,
    d u_val, d v_val, vec Un_val, dimVec u_grad, dimVec v_grad,
    vecDimVec Un_grad, vec U_bnd_val, Point<DIM> quadPoint, Point<DIM> normal, NumFlux* num_flux)
{
    d result = 0.;

    vec numFlux(COMPONENT_COUNT);

    vec bc_state(COMPONENT_COUNT);
    vec u_state(COMPONENT_COUNT);
    for (ui i = 0; i < COMPONENT_COUNT; i++)
    {
      bc_state = 0.;
      u_state = u_val;
    }

    num_flux->calculate(u_state, bc_state, quadPoint, normal, numFlux);

    result = numFlux[0] * v_val;

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

  vec u_state(COMPONENT_COUNT);
  vec u_stateN(COMPONENT_COUNT);
  for (ui i = 0; i < COMPONENT_COUNT; i++)
  {
    if (u_N)
    {
      u_stateN = u_val;
      u_state = 0.;
    }
    else
    {
      u_state = u_val;
      u_stateN = 0.;
    }
  }

  num_flux->calculate(u_state, u_stateN, quadPoint, normal, numFlux);

  result = numFlux[comp_i] * jump_v;

  return result;
}



d EquationImplementation::rhsVolValue(ui comp_i,
    d v_val, vec Un_val, dimVec v_grad,
    vecDimVec Un_grad, Point<DIM> quadPoint)
{
    d result = 0.;
    d flux;

    // Time derivative.
    result += Un_val[comp_i] * v_val / DELTA_T;

    return result;
}


d EquationImplementation::rhsBoundaryEdgeValue(ui comp_i,
    d v_val, vec Un_val, dimVec v_grad,
    vecDimVec Un_grad, vec U_bnd_val, Point<DIM> quadPoint, Point<DIM> normal, NumFlux* num_flux, DirichletBoundaryCondition* bc)
{
    d result = 0.;

    vec bc_state(COMPONENT_COUNT);
    vec u_state(COMPONENT_COUNT);
    for (ui i = 0; i < COMPONENT_COUNT; i++)
    {
      u_state = 0.;
      bc_state = bc->calculate(i, quadPoint);
    }
    vec numFlux(COMPONENT_COUNT);

    num_flux->calculate(u_state, bc_state, quadPoint, normal, numFlux);

    result = -numFlux[0] * v_val;

    return result;
}


d EquationImplementation::rhsInternalEdgeValue(ui comp_i,
    d v_val, dimVec v_grad, bool v_N, vec Un_val,
    vecDimVec Un_grad, vec Un_valN, vecDimVec Un_gradN, Point<DIM> quadPoint, Point<DIM> normal, NumFlux* num_flux)
{
  d result = 0.;

  d jump_v = v_N ? -v_val : v_val;

  vec numFlux(COMPONENT_COUNT);

  num_flux->calculate(Un_val, Un_valN, quadPoint, normal, numFlux);

  result = numFlux[comp_i] * jump_v;

  return -result;
}
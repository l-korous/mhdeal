#include <deal.II/lac/full_matrix.h>

#include "equationImplementation.h"
#include "boundaryConditions.h"
#include "numericalFlux.h"

using namespace dealii;

EquationImplementation::EquationImplementation()
{
  //   for(int i=0;i<Ne;i++){
  //     F[0][i]=F[1][i]=F[2][i]=0.;
  //     for(int j=0;j<Ne;j++){
  //       A[0][i][j]=0.;
  //       A[1][i][j]=0.;
  //       A[2][i][j]=0.;
  //     }
  //   }
}

d EquationImplementation::matrixVolValue(ui comp_i, ui comp_j,
                                         d u_val, d v_val, vec Un_val, dimVec u_grad, dimVec v_grad,
                                         vecDimVec Un_grad, Point<DIM> quadPoint)
{
  d result = 0.;

  // Time derivative.
  if (comp_i == comp_j)
    result += u_val * v_val / DELTA_T;

  d r = Un_val[0];
  d p1 = Un_val[1];
  d p2 = Un_val[2];
  d p3 = Un_val[3];
  d E = Un_val[4];
  d rum;
  switch (comp_i)
  {
  case 0:
    switch (comp_j)
    {
    case 0:
      result -= u_val * (A_1_1_1 * v_grad[0] + A_2_1_1 * v_grad[1] + A_3_1_1 * v_grad[2]);
      break;
    case 1:
      result -= u_val * (A_1_1_2 * v_grad[0] + A_2_1_2 * v_grad[1] + A_3_1_2 * v_grad[2]);
      break;
    case 2:
      result -= u_val * (A_1_1_3 * v_grad[0] + A_2_1_3 * v_grad[1] + A_3_1_3 * v_grad[2]);
      break;
    case 3:
      result -= u_val * (A_1_1_4 * v_grad[0] + A_2_1_4 * v_grad[1] + A_3_1_4 * v_grad[2]);
      break;
    case 4:
      result -= u_val * (A_1_1_5 * v_grad[0] + A_2_1_5 * v_grad[1] + A_3_1_5 * v_grad[2]);
      break;
    }
    break;
  case 1:
    switch (comp_j)
    {
    case 0:
      result -= u_val * (A_1_2_1 * v_grad[0] + A_2_2_1 * v_grad[1] + A_3_2_1 * v_grad[2]);
      break;
    case 1:
      result -= u_val * (A_1_2_2 * v_grad[0] + A_2_2_2 * v_grad[1] + A_3_2_2 * v_grad[2]);
      break;
    case 2:
      result -= u_val * (A_1_2_3 * v_grad[0] + A_2_2_3 * v_grad[1] + A_3_2_3 * v_grad[2]);
      break;
    case 3:
      result -= u_val * (A_1_2_4 * v_grad[0] + A_2_2_4 * v_grad[1] + A_3_2_4 * v_grad[2]);
      break;
    case 4:
      result -= u_val * (A_1_2_5 * v_grad[0] + A_2_2_5 * v_grad[1] + A_3_2_5 * v_grad[2]);
      break;
    }
    break;
  case 2:
    switch (comp_j)
    {
    case 0:
      result -= u_val * (A_1_3_1 * v_grad[0] + A_2_3_1 * v_grad[1] + A_3_3_1 * v_grad[2]);
      break;
    case 1:
      result -= u_val * (A_1_3_2 * v_grad[0] + A_2_3_2 * v_grad[1] + A_3_3_2 * v_grad[2]);
      break;
    case 2:
      result -= u_val * (A_1_3_3 * v_grad[0] + A_2_3_3 * v_grad[1] + A_3_3_3 * v_grad[2]);
      break;
    case 3:
      result -= u_val * (A_1_3_4 * v_grad[0] + A_2_3_4 * v_grad[1] + A_3_3_4 * v_grad[2]);
      break;
    case 4:
      result -= u_val * (A_1_3_5 * v_grad[0] + A_2_3_5 * v_grad[1] + A_3_3_5 * v_grad[2]);
      break;
    }
    break;
  case 3:
    switch (comp_j)
    {
    case 0:
      result -= u_val * (A_1_4_1 * v_grad[0] + A_2_4_1 * v_grad[1] + A_3_4_1 * v_grad[2]);
      break;
    case 1:
      result -= u_val * (A_1_4_2 * v_grad[0] + A_2_4_2 * v_grad[1] + A_3_4_2 * v_grad[2]);
      break;
    case 2:
      result -= u_val * (A_1_4_3 * v_grad[0] + A_2_4_3 * v_grad[1] + A_3_4_3 * v_grad[2]);
      break;
    case 3:
      result -= u_val * (A_1_4_4 * v_grad[0] + A_2_4_4 * v_grad[1] + A_3_4_4 * v_grad[2]);
      break;
    case 4:
      result -= u_val * (A_1_4_5 * v_grad[0] + A_2_4_5 * v_grad[1] + A_3_4_5 * v_grad[2]);
      break;
    }
    break;
  case 4:
    switch (comp_j)
    {
    case 0:
      result -= u_val * (A_1_5_1 * v_grad[0] + A_2_5_1 * v_grad[1] + A_3_5_1 * v_grad[2]);
      break;
    case 1:
      rum = -A_1_5_2;
      result -= u_val * (A_1_5_2 * v_grad[0] + A_2_5_2 * v_grad[1] + A_3_5_2 * v_grad[2]);
      break;
    case 2:
      result -= u_val * (A_1_5_3 * v_grad[0] + A_2_5_3 * v_grad[1] + A_3_5_3 * v_grad[2]);
      break;
    case 3:
      result -= u_val * (A_1_5_4 * v_grad[0] + A_2_5_4 * v_grad[1] + A_3_5_4 * v_grad[2]);
      break;
    case 4:
      result -= u_val * (A_1_5_5 * v_grad[0] + A_2_5_5 * v_grad[1] + A_3_5_5 * v_grad[2]);
      break;
    }
    break;
  }

  return result;
}

d EquationImplementation::rhsVolValue(ui comp_i,
                                      d v_val, vec Un_val, dimVec v_grad,
                                      vecDimVec Un_grad, Point<DIM> quadPoint)
{
  d result = 0.;

  // Time derivative.
  result += Un_val[comp_i] * v_val / DELTA_T;

  return result;
}


d EquationImplementation::rhsBoundaryEdgeValue(ui comp_i,
                                               d v_val, vec Un_val, dimVec v_grad,
                                               vecDimVec Un_grad, vec U_bnd_val, Point<DIM> quadPoint, Point<DIM> normal,
                                               NumFlux* num_flux, DirichletBoundaryCondition* bc, dealii::types::boundary_id bnd_id)
{
  d result = 0.;

  if(std::abs(normal(2))>1e-6)
    return 0.;

  // Pro proudovou hustotu se nic na hranicich nedela.
  if (comp_i >= COMPONENT_COUNT_T)
    return result;

  vec bc_state(COMPONENT_COUNT_T);

  for (ui i = 0; i < COMPONENT_COUNT_T; i++)
    bc_state[i] = bc->calculate(i, quadPoint);

  vec numFlux(COMPONENT_COUNT_T);


  num_flux->calculate(Un_val, bc_state, quadPoint, normal, numFlux);

  result -= DELTA_T*numFlux[comp_i] * v_val;

  return result;
}

d EquationImplementation::rhsInternalEdgeValue(ui comp_i,
                                               d v_val, dimVec v_grad, bool v_N, vec Un_val,
                                               vecDimVec Un_grad, vec Un_valN, vecDimVec Un_gradN, Point<DIM> quadPoint, Point<DIM> normal, NumFlux* num_flux)
{
  // Pro proudovou hustotu se nic na hranicich nedela.
  d result = 0.;
  if (comp_i >= COMPONENT_COUNT_T)
    return result;

  d jump_v = v_N ? -v_val : v_val;

  vec numFlux(COMPONENT_COUNT_T);

  num_flux->calculate(Un_val, Un_valN, quadPoint, normal, numFlux);

  result -= DELTA_T*numFlux[comp_i] * jump_v;

  return result;
}

d EquationImplementation::currentTime;

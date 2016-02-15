#include <deal.II/lac/full_matrix.h>

#include "equationImplementation.h"
#include "boundaryConditions.h"
#include "numericalFlux.h"

using namespace dealii;

d calculate_flux(double x, double y, double nx, double ny)
{
  return nx + ny;
}

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


d EquationImplementation::matrixBoundaryEdgeValue(ui comp_i, ui comp_j,
  d u_val, d v_val, vec Un_val, dimVec u_grad, dimVec v_grad,
  vecDimVec Un_grad, vec U_bnd_val, Point<DIM> quadPoint, Point<DIM> normal, NumFlux* num_flux, DirichletBoundaryCondition* bc, dealii::types::boundary_id bnd_id)
{
    d result = 0.;

    // Normal
    d nx = normal(0), ny = normal(1), nz = normal(2);

    if (BC_INFLOW_OUTFLOW(bnd_id))
    {
      // Use the right kind of num flux
      NumFluxVijayasundaram* num_flux_Vijayasundaram = (NumFluxVijayasundaram*)num_flux;

      // BC state
      d bc_state[COMPONENT_COUNT];
      for (ui i = 0; i < COMPONENT_COUNT; i++)
        bc_state[i] = bc->calculate(i, quadPoint);

      // Utility declaration & initialization
      double w_L[COMPONENT_COUNT], eigenvalues[COMPONENT_COUNT], alpha[COMPONENT_COUNT], q_ji_star[COMPONENT_COUNT], beta[COMPONENT_COUNT], q_ji[COMPONENT_COUNT], w_ji[COMPONENT_COUNT];
      double T[COMPONENT_COUNT][COMPONENT_COUNT];
      double T_inv[COMPONENT_COUNT][COMPONENT_COUNT];

      for (ui i = 0; i < COMPONENT_COUNT; i++)
      {
        w_L[i] = Un_val[i];
        alpha[i] = 0.;
        beta[i] = 0.;
        q_ji[i] = 0.;
        w_ji[i] = 0.;
        eigenvalues[i] = 0.;
        for (ui j = 0; j < COMPONENT_COUNT; j++)
        {
          T[i][j] = 0.;
          T_inv[i][j] = 0.;
        }
      }

      // Transformation of the inner state to the local coordinates.
      num_flux_Vijayasundaram->Q(num_flux_Vijayasundaram->q, w_L, nx, ny, nz);

      // Calculate Lambda^-.
      num_flux_Vijayasundaram->Lambda(eigenvalues);
      num_flux_Vijayasundaram->T(T);
      num_flux_Vijayasundaram->T_inv(T_inv);

      // Rotate the BC state
      num_flux_Vijayasundaram->Q(q_ji_star, bc_state, nx, ny, nz);

      for (unsigned int ai = 0; ai < COMPONENT_COUNT; ai++)
        for (unsigned int aj = 0; aj < COMPONENT_COUNT; aj++)
          alpha[ai] += T_inv[ai][aj] * num_flux_Vijayasundaram->q[aj];

      for (unsigned int bi = 0; bi < COMPONENT_COUNT; bi++)
        for (unsigned int bj = 0; bj < COMPONENT_COUNT; bj++)
          beta[bi] += T_inv[bi][bj] * q_ji_star[bj];

      for (unsigned int si = 0; si < COMPONENT_COUNT; si++)
        for (unsigned int sj = 0; sj < COMPONENT_COUNT; sj++)
          if (eigenvalues[sj] < 0)
            q_ji[si] += beta[sj] * T[si][sj];
          else
            q_ji[si] += alpha[sj] * T[si][sj];

      num_flux_Vijayasundaram->Q_inv(w_ji, q_ji, nx, ny, nz);

      d e[5][5] =
      {
        { 1., 0, 0, 0, 0 },
        { 0, 1., 0, 0, 0 },
        { 0, 0, 1., 0, 0 },
        { 0, 0, 0, 1., 0 },
        { 0, 0, 0, 0, 1. }
      };

      double P_plus[COMPONENT_COUNT];
      d w_temp[COMPONENT_COUNT];

      for (ui i = 0; i < COMPONENT_COUNT; i++)
        w_temp[i] = (w_ji[i] + w_L[i]) / 2.;

      num_flux_Vijayasundaram->P_plus(P_plus, w_temp, e[comp_j], nx, ny, nz);

      result = P_plus[comp_i] * u_val * v_val;
    }

    // Solid wall
    if (BC_SOLID_WALL(bnd_id))
    {
        d rho = Un_val[0];
        d v_1 = Un_val[1] / rho;
        d v_2 = Un_val[2] / rho;
        d v_3 = Un_val[3] / rho;

        double P[5][5];
        for (unsigned int P_i = 0; P_i < 5; P_i++)
            for (unsigned int P_j = 0; P_j < 5; P_j++)
                P[P_i][P_j] = 0.;

        P[1][0] = (KAPPA - 1.) * (v_1 * v_1 + v_2 * v_2 + v_3 * v_3) * nx / 2.;
        P[1][1] = (KAPPA - 1.) * (-v_1) * nx;
        P[1][2] = (KAPPA - 1.) * (-v_2) * nx;
        P[1][3] = (KAPPA - 1.) * (-v_3) * nx;
        P[1][4] = (KAPPA - 1.) * nx;

        P[2][0] = (KAPPA - 1.) * (v_1 * v_1 + v_2 * v_2 + v_3 * v_3) * ny / 2.;
        P[2][1] = (KAPPA - 1.) * (-v_1) * ny;
        P[2][2] = (KAPPA - 1.) * (-v_2) * ny;
        P[2][3] = (KAPPA - 1.) * (-v_3) * ny;
        P[2][4] = (KAPPA - 1.) * ny;

        P[3][0] = (KAPPA - 1.) * (v_1 * v_1 + v_2 * v_2 + v_3 * v_3) * nz / 2.;
        P[3][1] = (KAPPA - 1.) * (-v_1) * nz;
        P[3][2] = (KAPPA - 1.) * (-v_2) * nz;
        P[3][3] = (KAPPA - 1.) * (-v_3) * nz;
        P[3][4] = (KAPPA - 1.) * nz;

        result = P[comp_i][comp_j] * u_val * v_val;
    }

    return result;
}


d EquationImplementation::matrixInternalEdgeValue(ui comp_i, ui comp_j,
  d u_val, d v_val, vec Un_val, vec Un_valN, dimVec u_grad, dimVec v_grad,
  vecDimVec Un_grad, vecDimVec Un_gradN, bool u_N, bool v_N,
  Point<DIM> quadPoint, Point<DIM> normal, NumFlux* num_flux)
{
    d result = 0.;

    d jump_v = v_N ? -v_val : v_val;

    result = num_flux->calculate(Un_val, Un_valN, quadPoint, normal, comp_i, comp_j, (u_N ? 2 : 1)) * u_val * jump_v;

    // std::cout << "coord [" << comp_i << ", " << comp_j << "], point [" << quadPoint(0) << ", " << quadPoint(1) << "], normal [" << normal(0) << ", " << normal(1) << "]: " << result << std::endl;
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

  if (BC_INFLOW_OUTFLOW(bnd_id))
  {
    // Normal
    d nx = normal(0), ny = normal(1), nz = normal(2);

    // Use the right kind of num flux
    NumFluxVijayasundaram* num_flux_Vijayasundaram = (NumFluxVijayasundaram*)num_flux;

    // BC state
    d bc_state[COMPONENT_COUNT];
    for (ui i = 0; i < COMPONENT_COUNT; i++)
      bc_state[i] = bc->calculate(i, quadPoint);

    // Utility declaration & initialization
    double w_L[COMPONENT_COUNT], eigenvalues[COMPONENT_COUNT], alpha[COMPONENT_COUNT], q_ji_star[COMPONENT_COUNT], beta[COMPONENT_COUNT], q_ji[COMPONENT_COUNT], w_ji[COMPONENT_COUNT];
    double T[COMPONENT_COUNT][COMPONENT_COUNT];
    double T_inv[COMPONENT_COUNT][COMPONENT_COUNT];

    for (ui i = 0; i < COMPONENT_COUNT; i++)
    {
      w_L[i] = Un_val[i];
      alpha[i] = 0.;
      beta[i] = 0.;
      q_ji[i] = 0.;
      w_ji[i] = 0.;
      eigenvalues[i] = 0.;
      for (ui j = 0; j < COMPONENT_COUNT; j++)
      {
        T[i][j] = 0.;
        T_inv[i][j] = 0.;
      }
    }

    // Transformation of the inner state to the local coordinates.
    num_flux_Vijayasundaram->Q(num_flux_Vijayasundaram->q, w_L, nx, ny, nz);

    // Calculate Lambda^-.
    num_flux_Vijayasundaram->Lambda(eigenvalues);
    num_flux_Vijayasundaram->T(T);
    num_flux_Vijayasundaram->T_inv(T_inv);

    // Rotate the BC state
    num_flux_Vijayasundaram->Q(q_ji_star, bc_state, nx, ny, nz);

    for (unsigned int ai = 0; ai < COMPONENT_COUNT; ai++)
      for (unsigned int aj = 0; aj < COMPONENT_COUNT; aj++)
        alpha[ai] += T_inv[ai][aj] * num_flux_Vijayasundaram->q[aj];

    for (unsigned int bi = 0; bi < COMPONENT_COUNT; bi++)
      for (unsigned int bj = 0; bj < COMPONENT_COUNT; bj++)
        beta[bi] += T_inv[bi][bj] * q_ji_star[bj];

    for (unsigned int si = 0; si < COMPONENT_COUNT; si++)
      for (unsigned int sj = 0; sj < COMPONENT_COUNT; sj++)
        if (eigenvalues[sj] < 0)
          q_ji[si] += beta[sj] * T[si][sj];
        else
          q_ji[si] += alpha[sj] * T[si][sj];

    num_flux_Vijayasundaram->Q_inv(w_ji, q_ji, nx, ny, nz);

    d P_minus[COMPONENT_COUNT];
    d w_temp[COMPONENT_COUNT];

    for (ui i = 0; i < COMPONENT_COUNT; i++)
      w_temp[i] = (w_ji[i] + w_L[i]) / 2.;

    num_flux_Vijayasundaram->P_minus(P_minus, w_temp, w_ji, nx, ny, nz);

    result -= P_minus[comp_i] * v_val;
  }
  else if (BC_IS_SYMMETRIC(bnd_id))
  {
    assert(false);
  }

  return result;
}

d EquationImplementation::rhsInternalEdgeValue(ui comp_i,
  d v_val, dimVec v_grad, bool v_N, vec Un_val,
  vecDimVec Un_grad, vec Un_valN, vecDimVec Un_gradN, Point<DIM> quadPoint, Point<DIM> normal, NumFlux* num_flux)
{
  d result = 0.;

  return result;
}


void EquationImplementation::JacobiM(double A[][COMPONENT_COUNT][COMPONENT_COUNT], 
                              std::vector<Vector<double> > lv,
                              const unsigned int qp)
{
  double v[COMPONENT_COUNT], iRh, iRh2, Uk, p, gmmo, Ec1, Ec2, Ec3, E1, E2, E3;

  // using shorter notation for old solution
  // order of the variables is following: rho, v(3), B(3), U, J(3)
  for (unsigned int i = 0; i < COMPONENT_COUNT; i++)
    v[i] = lv[qp](i);

  iRh = 1.0 / v[0];
  iRh2 = iRh*iRh;
  Uk = iRh*(v[1] * v[1] + v[2] * v[2] + v[3] * v[3]);
  p = GAMMA*(v[7] - (v[4] * v[4] + v[5] * v[5] + v[6] * v[6]) - Uk);
  gmmo = GAMMA - 1.0;
  Ec1 = (v[3] * v[5] - v[2] * v[6])*iRh;
  Ec2 = (v[1] * v[6] - v[3] * v[4])*iRh;
  Ec3 = (v[2] * v[4] - v[1] * v[5])*iRh;
  E1 = Ec1 + ETA*v[8];
  E2 = Ec2 + ETA*v[9];
  E3 = Ec3 + ETA*v[10];

  //A[0][0][0] = 0;
  A[0][0][1] = 1;
  //A[0][0][2] = 0;
  //A[0][0][3] = 0;
  //A[0][0][4] = 0;
  //A[0][0][5] = 0;
  //A[0][0][6] = 0;
  //A[0][0][7] = 0;
  //A[0][0][8] = 0;
  //A[0][0][9] = 0;
  //A[0][0][10] = 0;

  A[0][1][0] = -(v[1] * v[1] * iRh2) + gmmo*Uk*.5*iRh;
  A[0][1][1] = (2  - gmmo)*v[1]*iRh;
  A[0][1][2] = -gmmo*v[2]*iRh;
  A[0][1][3] = -gmmo*v[3]*iRh;
  A[0][1][4] = -GAMMA*v[4];
  A[0][1][5] = (1 - gmmo)*v[5];
  A[0][1][6] = (1 - gmmo)*v[6];
  A[0][1][7] = 0.5*gmmo;
  //A[0][1][8] = 0;
  //A[0][1][9] = 0;
  //A[0][1][10] = 0;

  A[0][2][0] = -v[1] * v[2]*iRh2;
  A[0][2][1] = v[2] * iRh;
  A[0][2][2] = v[1] * iRh;
  //A[0][2][3] = 0;
  A[0][2][4] = -v[5];
  A[0][2][5] = -v[4];
  //A[0][2][6] = 0;
  //A[0][2][7] = 0;
  //A[0][2][8] = 0;
  //A[0][2][9] = 0;
  //A[0][2][10] = 0;

  A[0][3][0] = -v[1] * v[3]*iRh2;
  A[0][3][1] = v[3] * iRh;
  //A[0][3][2] = 0;
  A[0][3][3] = v[1] * iRh;
  A[0][3][4] = -v[6];
  //A[0][3][5] = 0;
  A[0][3][6] = -v[4];
  //A[0][3][7] = 0;
  //A[0][3][8] = 0;
  //A[0][3][9] = 0;
  //A[0][3][10] = 0;

  //A[0][4][0] = 0;
  //A[0][4][1] = 0;
  //A[0][4][2] = 0;
  //A[0][4][3] = 0;
  //A[0][4][4] = 0;
  //A[0][4][5] = 0;
  //A[0][4][6] = 0;
  //A[0][4][7] = 0;
  //A[0][4][8] = 0;
  //A[0][4][9] = 0;
  //A[0][4][10] = 0;

  A[0][5][0] = Ec3*iRh;
  A[0][5][1] = v[5] * iRh;
  A[0][5][2] =-v[4] * iRh;
  //A[0][5][3] = 0;
  A[0][5][4] =-v[2] * iRh;
  A[0][5][5] = v[1] * iRh;
  //A[0][5][6] = 0;
  //A[0][5][7] = 0;
  //A[0][5][8] = 0;
  //A[0][5][9] = 0;
  A[0][5][10] = -ETA;

  A[0][6][0] = -Ec2*iRh;
  A[0][6][1] = v[6] * iRh;
  //A[0][6][2] = 0;
  A[0][6][3] = -v[4] * iRh;
  A[0][6][4] = -v[3] * iRh;
  //A[0][6][5] = 0;
  A[0][6][6] = v[1] * iRh;
  //A[0][6][7] = 0;
  //A[0][6][8] = 0;
  A[0][6][9] = ETA;
  //A[0][6][10] = 0;

  A[0][7][0] = 2 * iRh*(v[5]*Ec3 - v[6]*Ec2) + v[1] * gmmo*Uk*iRh2 - v[1]*(Uk + p)*iRh2;
  A[0][7][1] = 2 * (v[5] * v[5] + v[6] * v[6])*iRh - 2*v[1] * gmmo*v[1]*iRh2 + (Uk + p)*iRh;
  A[0][7][2] = -2 * v[4] * v[5] * iRh - v[1] * 2 * gmmo*v[2]*iRh2;
  A[0][7][3] = -2 * v[4] * v[6] * iRh - v[1] * 2 * gmmo*v[3]*iRh2;
  A[0][7][4] = -2 * GAMMA*v[4] * v[1] * iRh + 2 * (-v[5] * v[2] - v[6] * v[3])*iRh;
  A[0][7][5] = -2 * GAMMA*v[5] * v[1] * iRh + 2 * (v[5] * v[1]*iRh - E3);
  A[0][7][6] = -2 * GAMMA*v[6] * v[1] * iRh + 2 * (v[6] * v[1]*iRh + E2);
  A[0][7][7] = GAMMA*v[1]*iRh;
  //A[0][7][8] = 0;
  A[0][7][9] = 2 * ETA*v[6];
  A[0][7][10] = -2 * ETA*v[5];

  //A[0][8][0] = 0;
  //A[0][8][1] = 0;
  //A[0][8][2] = 0;
  //A[0][8][3] = 0;
  //A[0][8][4] = 0;
  //A[0][8][5] = 0;
  //A[0][8][6] = 0;
  //A[0][8][7] = 0;
  //A[0][8][8] = 0;
  //A[0][8][9] = 0;
  //A[0][8][10] = 0;

  //A[0][9][0] = 0;
  //A[0][9][1] = 0;
  //A[0][9][2] = 0;
  //A[0][9][3] = 0;
  //A[0][9][4] = 0;
  //A[0][9][5] = 0;
  A[0][9][6] = 1;
  //A[0][9][7] = 0;
  //A[0][9][8] = 0;
  //A[0][9][9] = 0;
  //A[0][9][10] = 0;

  //A[0][10][0] = 0;
  //A[0][10][1] = 0;
  //A[0][10][2] = 0;
  //A[0][10][3] = 0;
  //A[0][10][4] = 0;
  A[0][10][5] = -1;
  //A[0][10][6] = 0;
  //A[0][10][7] = 0;
  //A[0][10][8] = 0;
  //A[0][10][9] = 0;
  //A[0][10][10] = 0;

  if (DIM>1){
  //A[1][0][0] = 0;
  //A[1][0][1] = 0;
  A[1][0][2] = 1;
  //A[1][0][3] = 0;
  //A[1][0][4] = 0;
  //A[1][0][5] = 0;
  //A[1][0][6] = 0;
  //A[1][0][7] = 0;
  //A[1][0][8] = 0;
  //A[1][0][9] = 0;
  //A[1][0][10] = 0;

  A[1][1][0] = -v[1] * v[2]*iRh2;
  A[1][1][1] = v[2] * iRh;
  A[1][1][2] = v[1] * iRh;
  //A[1][1][3] = 0;
  A[1][1][4] = -v[5];
  A[1][1][5] = -v[4];
  //A[1][1][6] = 0;
  //A[1][1][7] = 0;
  //A[1][1][8] = 0;
  //A[1][1][9] = 0;
  //A[1][1][10] = 0;

  A[1][2][0] = -v[2] * v[2] * iRh2 + gmmo*Uk*.5*iRh;
  A[1][2][1] = -gmmo*v[1]*iRh;
  A[1][2][2] = (2 - gmmo)*v[2]*iRh;
  A[1][2][3] = -gmmo*v[3]*iRh;
  A[1][2][4] = (1 - gmmo)*v[4];
  A[1][2][5] = -GAMMA*v[5];
  A[1][2][6] = (1 - gmmo)*v[6];
  A[1][2][7] = 0.5*gmmo;
  //A[1][2][8] = 0;
  //A[1][2][9] = 0;
  //A[1][2][10] = 0;

  A[1][3][0] = -v[2] * v[3]*iRh2;
  //A[1][3][1] = 0;
  A[1][3][2] = v[3] * iRh;
  A[1][3][3] = v[2] * iRh;
  //A[1][3][4] = 0;
  A[1][3][5] = -v[6];
  A[1][3][6] = -v[5];
  //A[1][3][7] = 0;
  //A[1][3][8] = 0;
  //A[1][3][9] = 0;
  //A[1][3][10] = 0;

  A[1][4][0] = -Ec3*iRh;
  A[1][4][1] =-v[5] * iRh;
  A[1][4][2] = v[4] * iRh;
  //A[1][4][3] = 0;
  A[1][4][4] = v[2] * iRh;
  A[1][4][5] =-v[1] * iRh;
  //A[1][4][6] = 0;
  //A[1][4][7] = 0;
  //A[1][4][8] = 0;
  //A[1][4][9] = 0;
  A[1][4][10] = ETA;

  //A[1][5][0] = 0;
  //A[1][5][1] = 0;
  //A[1][5][2] = 0;
  //A[1][5][3] = 0;
  //A[1][5][4] = 0;
  //A[1][5][5] = 0;
  //A[1][5][6] = 0;
  //A[1][5][7] = 0;
  //A[1][5][8] = 0;
  //A[1][5][9] = 0;
  //A[1][5][10] = 0;

  A[1][6][0] = Ec1*iRh;
  //A[1][6][1] = 0;
  A[1][6][2] = v[6] * iRh;
  A[1][6][3] =-v[5] * iRh;
  //A[1][6][4] = 0;
  A[1][6][5] =-v[3] * iRh;
  A[1][6][6] = v[2] * iRh;
  //A[1][6][7] = 0;
  A[1][6][8] = -ETA;
  //A[1][6][9] = 0;
  //A[1][6][10] = 0;

  A[1][7][0] = 2 * iRh*(-v[4]*Ec3 + v[6]*Ec1) + v[2] * gmmo*Uk*iRh2 - v[2]*(Uk + p)*iRh2;
  A[1][7][1] = -2 * v[4] * v[5] * iRh - 2 * gmmo*v[1] * v[2] * iRh2;
  A[1][7][2] = 2 * (v[4] * v[4] + v[6] * v[6])*iRh - 2*v[2] * gmmo*v[2]*iRh2 + (Uk + p)*iRh;
  A[1][7][3] = -2 * v[5] * v[6] * iRh - 2 * gmmo*v[2] * v[3] * iRh2;
  A[1][7][4] = -2 * GAMMA*v[4] * v[2] * iRh + 2 * (v[4] * v[2]*iRh + E3);
  A[1][7][5] = -2 * GAMMA*v[5] * v[2] * iRh + 2 * (-v[4] * v[1] - v[6] * v[3])*iRh;
  A[1][7][6] = -2 * GAMMA*v[6] * v[2] * iRh + 2 * (v[6] * v[2]*iRh - E1);
  A[1][7][7] = GAMMA*v[2]*iRh;
  A[1][7][8] = -2 * ETA*v[6];
  //A[1][7][9] = 0;
  A[1][7][10] = 2 * ETA*v[4];

  //A[1][8][0] 0;
  //A[1][8][1] 0;
  //A[1][8][2] 0;
  //A[1][8][3] 0;
  //A[1][8][4] 0;
  //A[1][8][5] 0;
  A[1][8][6] = -1;
  //A[1][8][7] 0;
  //A[1][8][8] 0;
  //A[1][8][9] 0;
  //A[1][8][10] 0;

  //A[1][9][0] = 0;
  //A[1][9][1] = 0;
  //A[1][9][2] = 0;
  //A[1][9][3] = 0;
  //A[1][9][4] = 0;
  //A[1][9][5] = 0;
  //A[1][9][6] = 0;
  //A[1][9][7] = 0;
  //A[1][9][8] = 0;
  //A[1][9][9] = 0;
  //A[1][9][10] = 0;

  //A[1][10][0] - 0;
  //A[1][10][1] - 0;
  //A[1][10][2] - 0;
  //A[1][10][3] - 0;
  A[1][10][4] = 1;
  //A[1][10][5] - 0;
  //A[1][10][6] - 0;
  //A[1][10][7] - 0;
  //A[1][10][8] - 0;
  //A[1][10][9] - 0;
  //A[1][10][10] - 0;
  }

  if (DIM>2){
  //A[2][0][0] = 0;
  //A[2][0][1] = 0;
  //A[2][0][2] = 0;
  A[2][0][3] = 1;
  //A[2][0][4] = 0;
  //A[2][0][5] = 0;
  //A[2][0][6] = 0;
  //A[2][0][7] = 0;
  //A[2][0][8] = 0;
  //A[2][0][9] = 0;
  //A[2][0][10] = 0;

  A[2][1][0] = -v[1] * v[3]*iRh2;
  A[2][1][1] = v[3] * iRh;
  //A[2][1][2] = 0;
  A[2][1][3] = v[1] * iRh;
  A[2][1][4] = -v[6];
  //A[2][1][5] = 0;
  A[2][1][6] = -v[4];
  //A[2][1][7] = 0;
  //A[2][1][8] = 0;
  //A[2][1][9] = 0;
  //A[2][1][10] = 0;

  A[2][2][0] = -v[2] * v[3]*iRh2;
  //A[2][2][1] = 0;
  A[2][2][2] = v[3] * iRh;
  A[2][2][3] = v[2] * iRh;
  //A[2][2][4] = 0;
  A[2][2][5] = -v[6];
  A[2][2][6] = -v[5];
  //A[2][2][7] = 0;
  //A[2][2][8] = 0;
  //A[2][2][9] = 0;
  //A[2][2][10] = 0;

  A[2][3][0] = -v[3] * v[3] * iRh2 + (gmmo*Uk)*.5*iRh;
  A[2][3][1] = -gmmo*v[1]*iRh;
  A[2][3][2] = -gmmo*v[2]*iRh;
  A[2][3][3] = (2 - gmmo)*v[3]*iRh;
  A[2][3][4] = (1 - gmmo)*v[4];
  A[2][3][5] = (1 - gmmo)*v[5];
  A[2][3][6] = -GAMMA*v[6];
  A[2][3][7] = 0.5*gmmo;
  //A[2][3][8] = 0;
  //A[2][3][9] = 0;
  //A[2][3][10] = 0;

  A[2][4][0] = Ec2*iRh;
  A[2][4][1] =-v[6] * iRh;
  //A[2][4][2] = 0;
  A[2][4][3] = v[4] * iRh;
  A[2][4][4] = v[3] * iRh;
  //A[2][4][5] = 0;
  A[2][4][6] =-v[1] * iRh;
  //A[2][4][7] = 0;
  //A[2][4][8] = 0;
  A[2][4][9] = -ETA;
  //A[2][4][10] = 0;

  A[2][5][0] = -Ec1*iRh;
  //A[2][5][1] = 0;
  A[2][5][2] =-v[6] * iRh;
  A[2][5][3] = v[5] * iRh;
  //A[2][5][4] = 0;
  A[2][5][5] = v[3] * iRh;
  A[2][5][6] =-v[2] * iRh;
  //A[2][5][7] = 0;
  A[2][5][8] = ETA;
  //A[2][5][9] = 0;
  //A[2][5][10] = 0;

  //A[2][6][0] = 0;
  //A[2][6][1] = 0;
  //A[2][6][2] = 0;
  //A[2][6][3] = 0;
  //A[2][6][4] = 0;
  //A[2][6][5] = 0;
  //A[2][6][6] = 0;
  //A[2][6][7] = 0;
  //A[2][6][8] = 0;
  //A[2][6][9] = 0;
  //A[2][6][10] = 0;

  A[2][7][0] = 2 * iRh*(v[4]*Ec2 - v[5]*Ec1) + v[3] * gmmo*Uk*iRh2 - v[3]*(Uk + p)*iRh2;
  A[2][7][1] = -2*v[4]*v[6]*iRh - 2*gmmo*v[1]*v[3]*iRh2;
  A[2][7][2] = -2*v[5]*v[6]*iRh - 2*gmmo*v[2]*v[3]*iRh2;
  A[2][7][3] = 2 * (v[4] * v[4] + v[5] * v[5])*iRh + 2*v[3] * gmmo*v[3]*iRh2 + (Uk + p)*iRh;
  A[2][7][4] = -2 * GAMMA*v[4] * v[3] * iRh + 2 * (v[4] * v[3]*iRh - E2);
  A[2][7][5] = -2 * GAMMA*v[5] * v[3] * iRh + 2 * (v[5] * v[3]*iRh + E1);
  A[2][7][6] = -2 * GAMMA*v[6] * v[3] * iRh + 2 * (-v[4] * v[1] - v[5] * v[2])*iRh;
  A[2][7][7] = GAMMA*v[3]*iRh;
  A[2][7][8] = 2 * ETA*v[5];
  A[2][7][9] = -2 * ETA*v[4];
  //A[2][7][10] = 0;

  //A[2][8][0] = 0;
  //A[2][8][1] = 0;
  //A[2][8][2] = 0;
  //A[2][8][3] = 0;
  //A[2][8][4] = 0;
  A[2][8][5] = 1;
  //A[2][8][6] = 0;
  //A[2][8][7] = 0;
  //A[2][8][8] = 0;
  //A[2][8][9] = 0;
  //A[2][8][10] = 0;

  //A[2][9][0] = 0;
  //A[2][9][1] = 0;
  //A[2][9][2] = 0;
  //A[2][9][3] = 0;
  A[2][9][4] = -1;
  //A[2][9][5] = 0;
  //A[2][9][6] = 0;
  //A[2][9][7] = 0;
  //A[2][9][8] = 0;
  //A[2][9][9] = 0;
  //A[2][9][10] = 0;

  //A[2][10][0] = 0;
  //A[2][10][1] = 0;
  //A[2][10][2] = 0;
  //A[2][10][3] = 0;
  //A[2][10][4] = 0;
  //A[2][10][5] = 0;
  //A[2][10][6] = 0;
  //A[2][10][7] = 0;
  //A[2][10][8] = 0;
  //A[2][10][9] = 0;
  //A[2][10][10] = 0;
    }
}

d EquationImplementation::currentTime;
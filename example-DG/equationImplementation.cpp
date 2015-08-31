#include "common.h"

// Not used so far
#pragma region FluxJacobians

#define Rcv (R / c_v)
#define v1 (p1 / r)
#define v2 (p2 / r)
#define v3 (p3 / r)
#define e (E / r)
#define Uk ((p1 * p1 + p2 * p2 + p3 * p3) / (2 * r * r))

// A_x
#define A_1_1_1 0.
#define A_1_1_2 1.
#define A_1_1_3 0.
#define A_1_1_4 0.
#define A_1_1_5 0.

#define A_1_2_1 -v1 * v1 + Rcv * Uk
#define A_1_2_2 2 * v1 - Rcv * v1
#define A_1_2_3 -Rcv * v2
#define A_1_2_4 -Rcv * v3
#define A_1_2_5 Rcv

#define A_1_3_1 -v1 * v2
#define A_1_3_2 v2
#define A_1_3_3 v1
#define A_1_3_4 0.
#define A_1_3_5 0.

#define A_1_4_1 -v1 * v3
#define A_1_4_2 v3
#define A_1_4_3 0.
#define A_1_4_4 v1
#define A_1_4_5 0.

#define A_1_5_1 -(v1 * E) - (v1 / r) * Rcv * (E - Uk) + v1 * Rcv * Uk
#define A_1_5_2 e + (1. / r) * Rcv * (E - Uk) - Rcv * Uk * v1 * v1
#define A_1_5_3 -Rcv * v1 * v2
#define A_1_5_4 -Rcv * v1 * v3
#define A_1_5_5 v1 + Rcv * v1

// A_y
#define A_2_1_1 0.
#define A_2_1_2 0.
#define A_2_1_3 1.
#define A_2_1_4 0.
#define A_2_1_5 0.

#define A_2_2_1 -v2 * v1
#define A_2_2_2 v2
#define A_2_2_3 v1
#define A_2_2_4 0.
#define A_2_2_5 0.

#define A_2_3_1 -v2 * v2 + Rcv * Uk
#define A_2_3_2 -Rcv * v1
#define A_2_3_3 2 * v2 - Rcv * v2
#define A_2_3_4 -Rcv * v3
#define A_2_3_5 Rcv

#define A_2_4_1 -v2 * v3
#define A_2_4_2 0.
#define A_2_4_3 v3
#define A_2_4_4 v2
#define A_2_4_5 0.

#define A_2_5_1 -(v2 * E) - (v2 / r) * Rcv * (E - Uk) + v2 * Rcv * Uk
#define A_2_5_2 -Rcv * v1 * v2
#define A_2_5_3 e + (1. / r) * Rcv * (E - Uk) - Rcv * Uk * v2 * v2
#define A_2_5_4 -Rcv * v2 * v3
#define A_2_5_5 v2 + Rcv * v2

// A_z
#define A_3_1_1 0.
#define A_3_1_2 0.
#define A_3_1_3 0.
#define A_3_1_4 1.
#define A_3_1_5 0.

#define A_3_2_1 -v3 * v1
#define A_3_2_2 v3
#define A_3_2_3 0.
#define A_3_2_4 v1
#define A_3_2_5 0.

#define A_3_3_1 -v3 * v2
#define A_3_3_2 0.
#define A_3_3_3 v3
#define A_3_3_4 v2
#define A_3_3_5 0.

#define A_3_4_1 -v3 * v3 + Rcv * Uk
#define A_3_4_2 -Rcv * v1
#define A_3_4_3 -Rcv * v2
#define A_3_4_4 2. * v3 - Rcv * v3
#define A_3_4_5 Rcv

#define A_3_5_1 -(v3 * E) - (v3 / r) * Rcv * (E - Uk) + v3 * Rcv * Uk
#define A_3_5_2 -Rcv * v1 * v3
#define A_3_5_3 -Rcv * v2 * v3
#define A_3_5_4 e + (1. / r) * Rcv * (E - Uk) - Rcv * Uk * v3 * v3
#define A_3_5_5 v3 + Rcv * v3

#pragma endregion

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

#if TIME_DISCRETIZATION_SEMI_IMPLICIT
        result -= calculate_flux(quadPoint(0), quadPoint(1), v_grad[0], v_grad[1]) * u_val;
#endif

    }

    return result;
}


d EquationImplementation::matrixBoundaryEdgeValue(ui comp_i, ui comp_j,
    d u_val, d v_val, vec Un_val, dimVec u_grad, dimVec v_grad,
    vecDimVec Un_grad, vec U_bnd_val, Point<DIM> quadPoint, Point<DIM> normal, NumFlux* num_flux)
{
    d result = 0.;

#if TIME_DISCRETIZATION_SEMI_IMPLICIT

    vec numFlux(COMPONENT_COUNT);

    vec bc_state(COMPONENT_COUNT);
    vec u_state(COMPONENT_COUNT);
    for (ui i = 0; i < COMPONENT_COUNT; i++)
    {
        bc_state = 0.;
        u_state = u_val;
    }

    // Jelikoz kdyz tok vyhodnoti ze je treba pouzit hodnotu nikoliv z elementu, ale ze souseda a jsme na hranici,
    // tak hodnoty z hranice musi jit na pravou stranu (nepouzije se bazova fce). Proto je tady boundary hodnota 0,
    // u prave strany (rhsBoundaryEdgeValue) to bude naopak.
    num_flux->calculate(u_state, bc_state, quadPoint, normal, numFlux);

    result += numFlux[0] * v_val;
#endif

    return result;
}


d EquationImplementation::matrixInternalEdgeValue(ui comp_i, ui comp_j,
    d u_val, d v_val, vec Un_val, vec Un_valN, dimVec u_grad, dimVec v_grad,
    vecDimVec Un_grad, vecDimVec Un_gradN, bool u_N, bool v_N,
    Point<DIM> quadPoint, Point<DIM> normal, NumFlux* num_flux)
{
    d result = 0.;

#if TIME_DISCRETIZATION_SEMI_IMPLICIT
    d jump_v = v_N ? -v_val : v_val;

    vec numFlux(COMPONENT_COUNT);

    vec u_state(COMPONENT_COUNT);
    vec u_stateN(COMPONENT_COUNT);

    // Tohle je takove asi matouci, ale takto muzeme mit pouze jednu funkci na vsechny pripady.
    // Vlastne vyuzivame toho, ze kdyz je u_N (== "u ma support na Neighborovi"), tak je u na aktualnim elementu 0,
    // a tok pocitam takto. Kdyz neni u_N, tak zase jak to odpovida.
    for (ui i = 0; i < COMPONENT_COUNT; i++)
    {
        u_state[i] = 0.;
        u_stateN[i] = 0.;

        if (u_N)
            u_stateN[comp_i] = u_val;
        else
            u_state[comp_i] = u_val;
    }

    num_flux->calculate(u_state, u_stateN, quadPoint, normal, numFlux);

    result += numFlux[comp_i] * jump_v;
#endif

    return result;
}



d EquationImplementation::rhsVolValue(ui comp_i,
    d v_val, vec Un_val, dimVec v_grad,
    vecDimVec Un_grad, Point<DIM> quadPoint)
{
    d result = 0.;

    // Time derivative.
    result += Un_val[comp_i] * v_val / DELTA_T;

#if !TIME_DISCRETIZATION_SEMI_IMPLICIT
    result += calculate_flux(quadPoint(0), quadPoint(1), v_grad[0], v_grad[1]) * Un_val[0];
#endif

    return result;
}


d EquationImplementation::rhsBoundaryEdgeValue(ui comp_i,
    d v_val, vec Un_val, dimVec v_grad,
    vecDimVec Un_grad, vec U_bnd_val, Point<DIM> quadPoint, Point<DIM> normal, NumFlux* num_flux, DirichletBoundaryCondition* bc)
{
    d result = 0.;

    vec bc_state(COMPONENT_COUNT);
    vec u_state(COMPONENT_COUNT);

    // Viz poznamka u matrixBoundaryEdgeValue.
    for (ui i = 0; i < COMPONENT_COUNT; i++)
    {
        u_state = 0.;
        bc_state = bc->calculate(i, quadPoint);
    }
    vec numFlux(COMPONENT_COUNT);

    num_flux->calculate(Un_val, bc_state, quadPoint, normal, numFlux);

    result -= numFlux[0] * v_val;

    return result;
}


d EquationImplementation::rhsInternalEdgeValue(ui comp_i,
    d v_val, dimVec v_grad, bool v_N, vec Un_val,
    vecDimVec Un_grad, vec Un_valN, vecDimVec Un_gradN, Point<DIM> quadPoint, Point<DIM> normal, NumFlux* num_flux)
{
    d result = 0.;

#if !TIME_DISCRETIZATION_SEMI_IMPLICIT

    d jump_v = v_N ? -v_val : v_val;

    vec numFlux(COMPONENT_COUNT);

    num_flux->calculate(Un_val, Un_valN, quadPoint, normal, numFlux);

    result -= numFlux[comp_i] * jump_v;

#endif

    return result;
}
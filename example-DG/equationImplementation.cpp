#include "common.h"

// Not used so far
#pragma region FluxJacobians

#define Rcv (R / c_v)
#define v1 (p1 / r)
#define v2 (p2 / r)
#define v3 (p3 / r)
#define e (E / r)
#define U_k ((p1 * p1 + p2 * p2 + p3 * p3) / (2 * r * r))

// A_x
#define A_1_1_1 0.
#define A_1_1_2 1.
#define A_1_1_3 0.
#define A_1_1_4 0.
#define A_1_1_5 0.

#define A_1_2_1 -v1 * v1 + Rcv * U_k
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

#define A_1_5_1 -(v1 * E) - (v1 / r) * Rcv * (E - U_k) + v1 * Rcv * U_k
#define A_1_5_2 e + (1. / r) * Rcv * (E - U_k) - Rcv * U_k * v1 * v1
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

#define A_2_3_1 -v2 * v2 + Rcv * U_k
#define A_2_3_2 -Rcv * v1
#define A_2_3_3 2 * v2 - Rcv * v2
#define A_2_3_4 -Rcv * v3
#define A_2_3_5 Rcv

#define A_2_4_1 -v2 * v3
#define A_2_4_2 0.
#define A_2_4_3 v3
#define A_2_4_4 v2
#define A_2_4_5 0.

#define A_2_5_1 -(v2 * E) - (v2 / r) * Rcv * (E - U_k) + v2 * Rcv * U_k
#define A_2_5_2 -Rcv * v1 * v2
#define A_2_5_3 e + (1. / r) * Rcv * (E - U_k) - Rcv * U_k * v2 * v2
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

#define A_3_4_1 -v3 * v3 + Rcv * U_k
#define A_3_4_2 -Rcv * v1
#define A_3_4_3 -Rcv * v2
#define A_3_4_4 2. * v3 - Rcv * v3
#define A_3_4_5 Rcv

#define A_3_5_1 -(v3 * E) - (v3 / r) * Rcv * (E - U_k) + v3 * Rcv * U_k
#define A_3_5_2 -Rcv * v1 * v3
#define A_3_5_3 -Rcv * v2 * v3
#define A_3_5_4 e + (1. / r) * Rcv * (E - U_k) - Rcv * U_k * v3 * v3
#define A_3_5_5 v3 + Rcv * v3

#pragma endregion

d calculate_flux(double x, double y, double nx, double ny)
{
    return nx + ny;
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

void EquationImplementation::Jacobians(FullMatrix<double> *J,
                                    std::vector<Vector<double> > lv,
                                    const unsigned int qp)
{
    double v[11], iRh, iRh2, Uk, p, gmmo, Ec1, Ec2, Ec3, E1, E2, E3;
  
  // using shorter notation for old solution
  // order of the variables is following: rho, v(3), B(3), U, J(3)
  for(unsigned int i=0;i<11;i++)
    v[i]=lv[qp](i);
  
  iRh=1.0/v[0];
  iRh2=iRh*iRh;
  Uk=iRh*(v[1]*v[1]+v[2]*v[2]+v[3]*v[3]);
  p=GAMMA*(v[7]-v[4]*v[4]+v[5]*v[5]+v[6]*v[6]-Uk);
  gmmo=GAMMA-1.0;
  Ec1=(v[3]*v[5]-v[2]*v[6])*iRh;
  Ec1=(v[1]*v[6]-v[3]*v[4])*iRh;
  Ec1=(v[2]*v[4]-v[1]*v[5])*iRh;
  E1=Ec1+ETA*v[8];
  E2=Ec2+ETA*v[9];
  E3=Ec3+ETA*v[10];
  
  J[0](0,0)=0;
  J[0](1,0)=1;
  J[0](2,0)=0;
  J[0](3,0)=0;
  J[0](4,0)=0;
  J[0](5,0)=0;
  J[0](6,0)=0;
  J[0](7,0)=0;
  J[0](8,0)=0;
  J[0](9,0)=0;
  J[0](10,0)=0;

  J[0](0,1)=-(v[1]*v[1]*iRh2)+(gmmo*Uk)*.5*iRh;
  J[0](1,1)=(2*v[1])*iRh-(gmmo*v[1])*iRh;
  J[0](2,1)=-((gmmo*v[2])*iRh);
  J[0](3,1)=-((gmmo*v[3])*iRh);
  J[0](4,1)=-2*v[4]+0.5*(2*v[4]-2*gmmo*v[4]);
  J[0](5,1)=0.5*(2*v[5]-2*gmmo*v[5]);
  J[0](6,1)=0.5*(2*v[6]-2*gmmo*v[6]);
  J[0](7,1)=0.5*gmmo;
  J[0](8,1)=0;
  J[0](9,1)=0;
  J[0](10,1)=0;

  J[0](0,2)=-((v[1]*v[2])*iRh2);
  J[0](1,2)=v[2]*iRh;
  J[0](2,2)=v[1]*iRh;
  J[0](3,2)=0;
  J[0](4,2)=-v[5];
  J[0](5,2)=-v[4];
  J[0](6,2)=0;
  J[0](7,2)=0;
  J[0](8,2)=0;
  J[0](9,2)=0;
  J[0](10,2)=0;

  J[0](0,3)=-((v[1]*v[3])*iRh2);
  J[0](1,3)=v[3]*iRh;
  J[0](2,3)=0;
  J[0](3,3)=v[1]*iRh;
  J[0](4,3)=-v[6];
  J[0](5,3)=0;
  J[0](6,3)=-v[4];
  J[0](7,3)=0;
  J[0](8,3)=0;
  J[0](9,3)=0;
  J[0](10,3)=0;

  J[0](0,4)=0;
  J[0](1,4)=0;
  J[0](2,4)=0;
  J[0](3,4)=0;
  J[0](4,4)=0;
  J[0](5,4)=0;
  J[0](6,4)=0;
  J[0](7,4)=0;
  J[0](8,4)=0;
  J[0](9,4)=0;
  J[0](10,4)=0;

  J[0](0,5)=Ec3*iRh;
  J[0](1,5)=v[5]*iRh;
  J[0](2,5)=-(v[4]*iRh);
  J[0](3,5)=0;
  J[0](4,5)=-(v[2]*iRh);
  J[0](5,5)=v[1]*iRh;
  J[0](6,5)=0;
  J[0](7,5)=0;
  J[0](8,5)=0;
  J[0](9,5)=0;
  J[0](10,5)=-ETA;

  J[0](0,6)=-Ec2*iRh;
  J[0](1,6)=v[6]*iRh;
  J[0](2,6)=0;
  J[0](3,6)=-(v[4]*iRh);
  J[0](4,6)=-(v[3]*iRh);
  J[0](5,6)=0;
  J[0](6,6)=v[1]*iRh;
  J[0](7,6)=0;
  J[0](8,6)=0;
  J[0](9,6)=ETA;
  J[0](10,6)=0;

  J[0](0,7)=2*iRh*(v[5]*Ec3-v[6]*Ec2)+v[1]*gmmo*Uk*iRh2-(v[1]*(Uk+p))*iRh2;
  J[0](1,7)=2*(v[5]*v[5]*iRh+v[6]*v[6]*iRh)+(v[1]*((2*v[1])*iRh-(2*GAMMA*v[1])*iRh))*iRh+(Uk+p)*iRh;
  J[0](2,7)=-(2*v[4]*v[5]*iRh)-(v[1]*2*gmmo*v[2])*iRh2;
  J[0](3,7)=-(2*v[4]*v[6]*iRh)-(v[1]*2*gmmo*v[3])*iRh2;
  J[0](4,7)=-(2*GAMMA*v[4]*v[1]*iRh)+2*(-((v[5]*v[2])*iRh)-(v[6]*v[3])*iRh);
  J[0](5,7)=-(2*GAMMA*v[5]*v[1]*iRh)+2*((v[5]*v[1])*iRh-E3);
  J[0](6,7)=-(2*GAMMA*v[6]*v[1]*iRh)+2*((v[6]*v[1])*iRh+E2);
  J[0](7,7)=(GAMMA*v[1])*iRh;
  J[0](8,7)=0;
  J[0](9,7)=2*ETA*v[6];
  J[0](10,7)=-2*ETA*v[5];

  J[0](0,8)=0;
  J[0](1,8)=0;
  J[0](2,8)=0;
  J[0](3,8)=0;
  J[0](4,8)=0;
  J[0](5,8)=0;
  J[0](6,8)=0;
  J[0](7,8)=0;
  J[0](8,8)=0;
  J[0](9,8)=0;
  J[0](10,8)=0;

  J[0](0,9)=0;
  J[0](1,9)=0;
  J[0](2,9)=0;
  J[0](3,9)=0;
  J[0](4,9)=0;
  J[0](5,9)=0;
  J[0](6,9)=1;
  J[0](7,9)=0;
  J[0](8,9)=0;
  J[0](9,9)=0;
  J[0](10,9)=0;

  J[0](0,10)=0;
  J[0](1,10)=0;
  J[0](2,10)=0;
  J[0](3,10)=0;
  J[0](4,10)=0;
  J[0](5,10)=-1;
  J[0](6,10)=0;
  J[0](7,10)=0;
  J[0](8,10)=0;
  J[0](9,10)=0;
  J[0](10,10)=0;

  
  
  J[1](0,0)=0;
  J[1](1,0)=0;
  J[1](2,0)=1;
  J[1](3,0)=0;
  J[1](4,0)=0;
  J[1](5,0)=0;
  J[1](6,0)=0;
  J[1](7,0)=0;
  J[1](8,0)=0;
  J[1](9,0)=0;
  J[1](10,0)=0;

  J[1](0,1)=-((v[1]*v[2])*iRh2);
  J[1](1,1)=v[2]*iRh;
  J[1](2,1)=v[1]*iRh;
  J[1](3,1)=0;
  J[1](4,1)=-v[5];
  J[1](5,1)=-v[4];
  J[1](6,1)=0;
  J[1](7,1)=0;
  J[1](8,1)=0;
  J[1](9,1)=0;
  J[1](10,1)=0;

  J[1](0,2)=-(v[2]*v[2]*iRh2)+(gmmo*Uk)*.5*iRh;
  J[1](1,2)=-((gmmo*v[1])*iRh);
  J[1](2,2)=(2*v[2])*iRh-(gmmo*v[2])*iRh;
  J[1](3,2)=-((gmmo*v[3])*iRh);
  J[1](4,2)=0.5*(2*v[4]-2*gmmo*v[4]);
  J[1](5,2)=-2*v[5]+0.5*(2*v[5]-2*gmmo*v[5]);
  J[1](6,2)=0.5*(2*v[6]-2*gmmo*v[6]);
  J[1](7,2)=0.5*gmmo;
  J[1](8,2)=0;
  J[1](9,2)=0;
  J[1](10,2)=0;

  J[1](0,3)=-((v[2]*v[3])*iRh2);
  J[1](1,3)=0;
  J[1](2,3)=v[3]*iRh;
  J[1](3,3)=v[2]*iRh;
  J[1](4,3)=0;
  J[1](5,3)=-v[6];
  J[1](6,3)=-v[5];
  J[1](7,3)=0;
  J[1](8,3)=0;
  J[1](9,3)=0;
  J[1](10,3)=0;

  J[1](0,4)=-(Ec3*iRh);
  J[1](1,4)=-(v[5]*iRh);
  J[1](2,4)=v[4]*iRh;
  J[1](3,4)=0;
  J[1](4,4)=v[2]*iRh;
  J[1](5,4)=-(v[1]*iRh);
  J[1](6,4)=0;
  J[1](7,4)=0;
  J[1](8,4)=0;
  J[1](9,4)=0;
  J[1](10,4)=ETA;

  J[1](0,5)=0;
  J[1](1,5)=0;
  J[1](2,5)=0;
  J[1](3,5)=0;
  J[1](4,5)=0;
  J[1](5,5)=0;
  J[1](6,5)=0;
  J[1](7,5)=0;
  J[1](8,5)=0;
  J[1](9,5)=0;
  J[1](10,5)=0;

  J[1](0,6)=Ec1*iRh;
  J[1](1,6)=0;
  J[1](2,6)=v[6]*iRh;
  J[1](3,6)=-(v[5]*iRh);
  J[1](4,6)=0;
  J[1](5,6)=-(v[3]*iRh);
  J[1](6,6)=v[2]*iRh;
  J[1](7,6)=0;
  J[1](8,6)=-ETA;
  J[1](9,6)=0;
  J[1](10,6)=0;

  J[1](0,7)=2*iRh*(-(v[4]*Ec3)+v[6]*Ec1)+v[2]*gmmo*Uk*iRh2-(v[2]*(Uk+p))*iRh2;
  J[1](1,7)=-(2*v[4]*v[5]*iRh)-2*gmmo*v[1]*v[2]*iRh2;
  J[1](2,7)=2*(v[4]*v[4]*iRh+v[6]*v[6]*iRh)+(v[2]*((2*v[2])*iRh-(2*GAMMA*v[2])*iRh))*iRh+(Uk+p)*iRh;
  J[1](3,7)=-(2*v[5]*v[6]*iRh)-2*gmmo*v[2]*v[3]*iRh2;
  J[1](4,7)=-(2*GAMMA*v[4]*v[2]*iRh)+2*((v[4]*v[2])*iRh+E3);
  J[1](5,7)=-(2*GAMMA*v[5]*v[2]*iRh)+2*(-((v[4]*v[1])*iRh)-(v[6]*v[3])*iRh);
  J[1](6,7)=-(2*GAMMA*v[6]*v[2]*iRh)+2*((v[6]*v[2])*iRh-E1);
  J[1](7,7)=(GAMMA*v[2])*iRh;
  J[1](8,7)=-2*ETA*v[6];
  J[1](9,7)=0;
  J[1](10,7)=2*ETA*v[4];

  J[1](0,8)=0;
  J[1](1,8)=0;
  J[1](2,8)=0;
  J[1](3,8)=0;
  J[1](4,8)=0;
  J[1](5,8)=0;
  J[1](6,8)=-1;
  J[1](7,8)=0;
  J[1](8,8)=0;
  J[1](9,8)=0;
  J[1](10,8)=0;

  J[1](0,9)=0;
  J[1](1,9)=0;
  J[1](2,9)=0;
  J[1](3,9)=0;
  J[1](4,9)=0;
  J[1](5,9)=0;
  J[1](6,9)=0;
  J[1](7,9)=0;
  J[1](8,9)=0;
  J[1](9,9)=0;
  J[1](10,9)=0;

  J[1](0,10)=0;
  J[1](1,10)=0;
  J[1](2,10)=0;
  J[1](3,10)=0;
  J[1](4,10)=1;
  J[1](5,10)=0;
  J[1](6,10)=0;
  J[1](7,10)=0;
  J[1](8,10)=0;
  J[1](9,10)=0;
  J[1](10,10)=0;

  
  
  J[2](0,0)=0;
  J[2](1,0)=0;
  J[2](2,0)=0;
  J[2](3,0)=1;
  J[2](4,0)=0;
  J[2](5,0)=0;
  J[2](6,0)=0;
  J[2](7,0)=0;
  J[2](8,0)=0;
  J[2](9,0)=0;
  J[2](10,0)=0;

  J[2](0,1)=-((v[1]*v[3])*iRh2);
  J[2](1,1)=v[3]*iRh;
  J[2](2,1)=0;
  J[2](3,1)=v[1]*iRh;
  J[2](4,1)=-v[6];
  J[2](5,1)=0;
  J[2](6,1)=-v[4];
  J[2](7,1)=0;
  J[2](8,1)=0;
  J[2](9,1)=0;
  J[2](10,1)=0;

  J[2](0,2)=-((v[2]*v[3])*iRh2);
  J[2](1,2)=0;
  J[2](2,2)=v[3]*iRh;
  J[2](3,2)=v[2]*iRh;
  J[2](4,2)=0;
  J[2](5,2)=-v[6];
  J[2](6,2)=-v[5];
  J[2](7,2)=0;
  J[2](8,2)=0;
  J[2](9,2)=0;
  J[2](10,2)=0;

  J[2](0,3)=-(v[3]*v[3]*iRh2)+(gmmo*Uk)*.5*iRh;
  J[2](1,3)=-((gmmo*v[1])*iRh);
  J[2](2,3)=-((gmmo*v[2])*iRh);
  J[2](3,3)=(2*v[3])*iRh-(gmmo*v[3])*iRh;
  J[2](4,3)=0.5*(2*v[4]-2*gmmo*v[4]);
  J[2](5,3)=0.5*(2*v[5]-2*gmmo*v[5]);
  J[2](6,3)=-2*v[6]+0.5*(2*v[6]-2*gmmo*v[6]);
  J[2](7,3)=0.5*gmmo;
  J[2](8,3)=0;
  J[2](9,3)=0;
  J[2](10,3)=0;

  J[2](0,4)=Ec2*iRh;
  J[2](1,4)=-(v[6]*iRh);
  J[2](2,4)=0;
  J[2](3,4)=v[4]*iRh;
  J[2](4,4)=v[3]*iRh;
  J[2](5,4)=0;
  J[2](6,4)=-(v[1]*iRh);
  J[2](7,4)=0;
  J[2](8,4)=0;
  J[2](9,4)=-ETA;
  J[2](10,4)=0;

  J[2](0,5)=-(Ec1*iRh);
  J[2](1,5)=0;
  J[2](2,5)=-(v[6]*iRh);
  J[2](3,5)=v[5]*iRh;
  J[2](4,5)=0;
  J[2](5,5)=v[3]*iRh;
  J[2](6,5)=-(v[2]*iRh);
  J[2](7,5)=0;
  J[2](8,5)=ETA;
  J[2](9,5)=0;
  J[2](10,5)=0;

  J[2](0,6)=0;
  J[2](1,6)=0;
  J[2](2,6)=0;
  J[2](3,6)=0;
  J[2](4,6)=0;
  J[2](5,6)=0;
  J[2](6,6)=0;
  J[2](7,6)=0;
  J[2](8,6)=0;
  J[2](9,6)=0;
  J[2](10,6)=0;

  J[2](0,7)=2*iRh*(v[4]*Ec2-v[5]*Ec1)+v[3]*gmmo*Uk*iRh2-(v[3]*(Uk+p))*iRh2;
  // J[2](1,7)=-((2*v[4]*v[6])*iRh)-2*gmmo*v[1]*v[3])*iRh2;
  // J[2](2,7)=-((2*v[5]*v[6])*iRh)-2*gmmo*v[2]*v[3])*iRh2;
  J[2](3,7)=2*(v[4]*v[4]*iRh+v[5]*v[5]*iRh)+(v[3]*((2*v[3])*iRh-(2*GAMMA*v[3])*iRh))*iRh+(Uk+p)*iRh;
  J[2](4,7)=-(2*GAMMA*v[4]*v[3]*iRh)+2*((v[4]*v[3])*iRh-E2);
  J[2](5,7)=-(2*GAMMA*v[5]*v[3]*iRh)+2*((v[5]*v[3])*iRh+E1);
  J[2](6,7)=2*(-((v[4]*v[1])*iRh)-(v[5]*v[2])*iRh)-(2*GAMMA*v[6]*v[3])*iRh;
  J[2](7,7)=(GAMMA*v[3])*iRh;
  J[2](8,7)=2*ETA*v[5];
  J[2](9,7)=-2*ETA*v[4];
  J[2](10,7)=0;

  J[2](0,8)=0;
  J[2](1,8)=0;
  J[2](2,8)=0;
  J[2](3,8)=0;
  J[2](4,8)=0;
  J[2](5,8)=1;
  J[2](6,8)=0;
  J[2](7,8)=0;
  J[2](8,8)=0;
  J[2](9,8)=0;
  J[2](10,8)=0;

  J[2](0,9)=0;
  J[2](1,9)=0;
  J[2](2,9)=0;
  J[2](3,9)=0;
  J[2](4,9)=-1;
  J[2](5,9)=0;
  J[2](6,9)=0;
  J[2](7,9)=0;
  J[2](8,9)=0;
  J[2](9,9)=0;
  J[2](10,9)=0;

  J[2](0,10)=0;
  J[2](1,10)=0;
  J[2](2,10)=0;
  J[2](3,10)=0;
  J[2](4,10)=0;
  J[2](5,10)=0;
  J[2](6,10)=0;
  J[2](7,10)=0;
  J[2](8,10)=0;
  J[2](9,10)=0;
  J[2](10,10)=0;

}
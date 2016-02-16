#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <deal.II/lac/vector.h>
#include <deal.II/base/point.h>

typedef unsigned int ui;
typedef double d;

// These are discretization definitions.
// This is preprocessor because it is used in templates.
#define DIM 3

extern const ui DG_ORDER;
extern const ui INIT_REF_NUM_X, INIT_REF_NUM_Y, INIT_REF_NUM_Z;
#define COMPONENT_COUNT 11
#define COMPONENT_COUNT_T 8
extern const ui TIME_DISCRETIZATION_SEMI_IMPLICIT;

// boundary id
extern const unsigned int BOUNDARY_FRONT;
extern const unsigned int BOUNDARY_RIGHT;
extern const unsigned int BOUNDARY_BACK;
extern const unsigned int BOUNDARY_LEFT;
extern const unsigned int BOUNDARY_BOTTOM;
extern const unsigned int BOUNDARY_TOP;

bool BC_INFLOW_OUTFLOW(const unsigned int bnd_marker);
bool BC_SOLID_WALL(const unsigned int bnd_marker);
bool BC_IS_SYMMETRIC(const unsigned int bnd_marker);

// Points defining geometry
extern const dealii::Point<DIM> p1;
extern const dealii::Point<DIM> p2;
extern const dealii::Point<DIM> p3;
extern const dealii::Point<DIM> p4;

extern const d T_FINAL;
extern const d DELTA_T;
extern const bool PRINT_ALGEBRA;
extern const bool DELETE_VTK;

typedef dealii::Tensor<1, DIM> dimVec;
typedef std::vector<dimVec> vecDimVec;
typedef dealii::Vector<d> vec;

// Problem definitions
extern const d GAMMA;
extern const d ETA;
extern const d R;
extern const d C_V;
extern const d KAPPA;

// Left
extern const double RHO_EXT;
extern const double V1_EXT;
extern const double V2_EXT;
extern const double V3_EXT;
extern const double P_EXT;
extern const double E_EXT;
extern const double B1_EXT;
extern const double B2_EXT;
extern const double B3_EXT;
extern const double J1_EXT;
extern const double J2_EXT;
extern const double J3_EXT;

// Not used so far
#pragma region FluxJacobians
// http://www.theoretical-physics.net/dev/fluid-dynamics/euler.html

#define Rcv (R / C_V)
#define v1_ (p1 / r)
#define v2_ (p2 / r)
#define v3_ (p3 / r)
#define en (E / r)
#define U_k ((p1 * p1 + p2 * p2 + p3 * p3) / (2. * r * r))

// A_x
#define A_1_1_1 (0.)
#define A_1_1_2 (1.)
#define A_1_1_3 (0.)
#define A_1_1_4 (0.)
#define A_1_1_5 (0.)

#define A_1_2_1 (-v1_ * v1_ + Rcv * U_k)
#define A_1_2_2 (2 * v1_ - Rcv * v1_)
#define A_1_2_3 (-Rcv * v2_)
#define A_1_2_4 (-Rcv * v3_)
#define A_1_2_5 (Rcv)

#define A_1_3_1 (-v1_ * v2_)
#define A_1_3_2 (v2_)
#define A_1_3_3 (v1_)
#define A_1_3_4 (0.)
#define A_1_3_5 (0.)

#define A_1_4_1 (-v1_ * v3_)
#define A_1_4_2 (v3_)
#define A_1_4_3 (0.)
#define A_1_4_4 (v1_)
#define A_1_4_5 (0.)

#define A_1_5_1 (-(v1_ * en) - (v1_ / r) * Rcv * (E - (r*U_k)) + v1_ * Rcv * U_k)
#define A_1_5_2 (en + (1. / r) * Rcv * (E - (r*U_k)) - (Rcv * v1_ * v1_))
#define A_1_5_3 (-Rcv * v1_ * v2_)
#define A_1_5_4 (-Rcv * v1_ * v3_)
#define A_1_5_5 (v1_ + (Rcv * v1_))

// A_y)
#define A_2_1_1 (0.)
#define A_2_1_2 (0.)
#define A_2_1_3 (1.)
#define A_2_1_4 (0.)
#define A_2_1_5 (0.)

#define A_2_2_1 (-v2_ * v1_)
#define A_2_2_2 (v2_)
#define A_2_2_3 (v1_)
#define A_2_2_4 (0.)
#define A_2_2_5 (0.)

#define A_2_3_1 (-v2_ * v2_ + Rcv * U_k)
#define A_2_3_2 (-Rcv * v1_)
#define A_2_3_3 (2 * v2_ - Rcv * v2_)
#define A_2_3_4 (-Rcv * v3_)
#define A_2_3_5 (Rcv)

#define A_2_4_1 (-v2_ * v3_)
#define A_2_4_2 (0.)
#define A_2_4_3 (v3_)
#define A_2_4_4 (v2_)
#define A_2_4_5 (0.)

#define A_2_5_1 (-(v2_ * en) - (v2_ / r) * Rcv * (E - (r*U_k)) + v2_ * Rcv * U_k)
#define A_2_5_2 (-Rcv * v1_ * v2_)
#define A_2_5_3 (en + (1. / r) * Rcv * (E - (r*U_k)) - (Rcv * v2_ * v2_))
#define A_2_5_4 (-Rcv * v2_ * v3_)
#define A_2_5_5 (v2_ + Rcv * v2_)

// A_z
#define A_3_1_1 (0.)
#define A_3_1_2 (0.)
#define A_3_1_3 (0.)
#define A_3_1_4 (1.)
#define A_3_1_5 (0.)

#define A_3_2_1 (-v3_ * v1_)
#define A_3_2_2 (v3_)
#define A_3_2_3 (0.)
#define A_3_2_4 (v1_)
#define A_3_2_5 (0.)

#define A_3_3_1 (-v3_ * v2_)
#define A_3_3_2 (0.)
#define A_3_3_3 (v3_)
#define A_3_3_4 (v2_)
#define A_3_3_5 (0.)

#define A_3_4_1 (-v3_ * v3_ + Rcv * U_k)
#define A_3_4_2 (-Rcv * v1_)
#define A_3_4_3 (-Rcv * v2_)
#define A_3_4_4 (2. * v3_ - Rcv * v3_)
#define A_3_4_5 (Rcv)

#define A_3_5_1 (-(v3_ * en) - (v3_ / r) * Rcv * (E - (r*U_k)) + v3_ * Rcv * U_k)
#define A_3_5_2 (-Rcv * v1_ * v3_)
#define A_3_5_3 (-Rcv * v2_ * v3_)
#define A_3_5_4 (en + (1. / r) * Rcv * (E - (r*U_k)) - (Rcv * v3_ * v3_))
#define A_3_5_5 (v3_ + Rcv * v3_)

#define f_1_1 (A_1_1_1 * r + A_1_1_2 * p1 + A_1_1_3 * p2 + A_1_1_4 * p3 + A_1_1_5 * E)
#define f_1_2 (A_1_2_1 * r + A_1_2_2 * p1 + A_1_2_3 * p2 + A_1_2_4 * p3 + A_1_2_5 * E)
#define f_1_3 (A_1_3_1 * r + A_1_3_2 * p1 + A_1_3_3 * p2 + A_1_3_4 * p3 + A_1_3_5 * E)
#define f_1_4 (A_1_4_1 * r + A_1_4_2 * p1 + A_1_4_3 * p2 + A_1_4_4 * p3 + A_1_4_5 * E)
#define f_1_5 (A_1_5_1 * r + A_1_5_2 * p1 + A_1_5_3 * p2 + A_1_5_4 * p3 + A_1_5_5 * E)

#define f_2_1 (A_2_1_1 * r + A_2_1_2 * p1 + A_2_1_3 * p2 + A_2_1_4 * p3 + A_2_1_5 * E)
#define f_2_2 (A_2_2_1 * r + A_2_2_2 * p1 + A_2_2_3 * p2 + A_2_2_4 * p3 + A_2_2_5 * E)
#define f_2_3 (A_2_3_1 * r + A_2_3_2 * p1 + A_2_3_3 * p2 + A_2_3_4 * p3 + A_2_3_5 * E)
#define f_2_4 (A_2_4_1 * r + A_2_4_2 * p1 + A_2_4_3 * p2 + A_2_4_4 * p3 + A_2_4_5 * E)
#define f_2_5 (A_2_5_1 * r + A_2_5_2 * p1 + A_2_5_3 * p2 + A_2_5_4 * p3 + A_2_5_5 * E)

#define f_3_1 (A_3_1_1 * r + A_3_1_2 * p1 + A_3_1_3 * p2 + A_3_1_4 * p3 + A_3_1_5 * E)
#define f_3_2 (A_3_2_1 * r + A_3_2_2 * p1 + A_3_2_3 * p2 + A_3_2_4 * p3 + A_3_2_5 * E)
#define f_3_3 (A_3_3_1 * r + A_3_3_2 * p1 + A_3_3_3 * p2 + A_3_3_4 * p3 + A_3_3_5 * E)
#define f_3_4 (A_3_4_1 * r + A_3_4_2 * p1 + A_3_4_3 * p2 + A_3_4_4 * p3 + A_3_4_5 * E)
#define f_3_5 (A_3_5_1 * r + A_3_5_2 * p1 + A_3_5_3 * p2 + A_3_5_4 * p3 + A_3_5_5 * E)

#pragma endregion

#endif

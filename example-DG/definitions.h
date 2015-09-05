#ifndef DEFINITIONS_H
#define DEFINITIONS_H

typedef unsigned int ui;
typedef double d;

// These are discretization definitions.
// This is preprocessor because it is used in templates.
#define DIM 2

extern const ui DG_ORDER;
extern const ui INIT_REF_NUM;
extern const ui COMPONENT_COUNT;
extern const ui TIME_DISCRETIZATION_SEMI_IMPLICIT;

extern const d T_FINAL;
extern const d DELTA_T;
extern const bool PRINT_ALGEBRA;

typedef dealii::Tensor<1, DIM> dimVec;
typedef std::vector<dimVec> vecDimVec;
typedef dealii::Vector<d> vec;

// Problem definitions
extern const d GAMMA;
extern const d ETA;

#endif
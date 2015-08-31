#ifndef DEFINITIONS_H
#define DEFINITIONS_H

// These are problem definitions.

#define DIM 2
#define DG_ORDER 1
#define INIT_REF_NUM 4
#define COMPONENT_COUNT 1

// 1 - SEMI-IMPLICIT
// 2 - EXPLICIT
#define TIME_DISCRETIZATION_SEMI_IMPLICIT 1

#define T_FINAL 1000.0
#define DELTA_T 0.1

const bool PRINT_ALGEBRA = false;

#define FLUX nx + ny

#endif
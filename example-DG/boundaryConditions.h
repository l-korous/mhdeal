#ifndef DIRICHLET_H
#define DIRICHLET_H

#include "definitions.h"

class DirichletBoundaryCondition
{
public:
    d calculate(ui component, dealii::Point<DIM> point);
};

#endif

#include "common.h"

d DirichletBoundaryCondition::calculate(ui component, Point<DIM> point)
{
    if (point(0) <= 0.501 && point(1) < 0.0005 && point(0) >=-0.001)
        return 1.0;
    else
        return 0.;
}

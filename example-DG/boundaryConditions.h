#ifndef DIRICHLET_H
#define DIRICHLET_H

class DirichletBoundaryCondition
{
public:
    d calculate(ui component, Point<DIM> point)
    {
        if (point(0) < 0.5 && point(1) < 0.0005)
            return 1.0;
        else
            return 0.;
    }
};

#endif
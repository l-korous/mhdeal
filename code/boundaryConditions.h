#ifndef _BOUNDARY_CONDITIONS_H
#define _BOUNDARY_CONDITIONS_H

#include "util.h"
#include "equations.h"

// Boundary conditions
template <EquationsType equationsType, int dim>
class BoundaryConditions
{
public:
  BoundaryConditions();
  static const int max_n_boundaries = 10;
  typename Equations<equationsType, dim>::BoundaryKind kind[max_n_boundaries][Equations<equationsType, dim>::n_components];
  void bc_vector_value(int boundary_no, const std::vector<Point<dim> > &points, std::vector<Vector<double> >&) const;
};

#endif
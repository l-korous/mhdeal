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

  // Maximum allowed number of boundaries.
  static const int max_n_boundaries = 10;

  // Pairing between boundary parts and the boundary types (BoundaryKinds) for the used equations.
  typename Equations<equationsType, dim>::BoundaryKind kind[max_n_boundaries][Equations<equationsType, dim>::n_components];
  
  // Values for this boundary identifier.
  void bc_vector_value(int boundary_no, const std::vector<Point<dim> > &points, std::vector<Vector<double> >&) const;
};

#endif
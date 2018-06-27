#ifndef _BOUNDARY_CONDITIONS_H
#define _BOUNDARY_CONDITIONS_H

#include "util.h"
#include "equations.h"
#include "parameters.h"

// Boundary conditions
template <EquationsType equationsType, int dim>
class BoundaryCondition
{
public:
  typedef std::array<double, Equations<equationsType, dim>::n_components> values_vector;
  typedef std::array<std::array<double, dim>, Equations<equationsType, dim>::n_components> grad_vector;

  BoundaryCondition(Parameters<dim>& parameters);

  // Values for this boundary identifier.
  virtual void bc_vector_value(int boundary_no, const Point<dim> &point, const Tensor<1, dim> &normal, values_vector &result,
    const grad_vector &grads, const values_vector &W_plus, double time, typename DoFHandler<dim>::active_cell_iterator&) const;

  // Passed as a constructor parameter
  Parameters<dim>& parameters;
};

#endif
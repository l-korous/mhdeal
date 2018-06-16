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
  typedef std::array<double, Equations<equationsType, dim>::n_components> InputVector;

  BoundaryCondition(Parameters<dim>& parameters);

  // Values for this boundary identifier.
  virtual void bc_vector_value(int boundary_no, const Point<dim> &point, InputVector &result, const InputVector &W_plus, double time) const;

  double compute_energy_from_given_pressure(const InputVector &W, double pressure) const;

  // Passed as a constructor parameter
  Parameters<dim>& parameters;
};

#endif
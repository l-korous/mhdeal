#ifndef _INITIAL_CONDITION_H
#define _INITIAL_CONDITION_H

#include "util.h"
#include "parameters.h"
#include "equations.h"

// Initial condition
template <EquationsType equationsType, int dim>
class InitialCondition
{
public:
  // InitialCondition constructor takes parameters as an attribute - to set up whatever necessary and use it in the initial values definition
  InitialCondition(Parameters<dim>&);

  // To be overwritten.
  virtual void vector_value(const std::vector<Point<dim> >& , std::vector<std::array<double, Equations<equationsType, dim>::n_components> >&) const = 0;

protected:
  Parameters<dim>& parameters;
};

#endif

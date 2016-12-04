#ifndef _INITIAL_CONDITION_H
#define _INITIAL_CONDITION_H

#include "util.h"
#include "parameters.h"
#include "equations.h"

// Initial condition
template <EquationsType equationsType, int dim>
class InitialCondition : public Function<dim>
{
public:
  InitialCondition(Parameters<dim>& parameters);
  double value(const Point<dim> &p, const unsigned int  component = 0) const;

private:
  Parameters<dim>& parameters;
};

#endif
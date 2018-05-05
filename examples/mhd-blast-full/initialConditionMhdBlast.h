#ifndef _INITIAL_CONDITION_MHD_BLAST_H
#define _INITIAL_CONDITION_MHD_BLAST_H

#include "util.h"
#include "parameters.h"
#include "initialCondition.h"
#include "equations.h"

// Initial condition
template <EquationsType equationsType, int dim>
class InitialConditionMhdBlast : public InitialCondition<equationsType, dim>
{
public:
  InitialConditionMhdBlast(Parameters<dim>&);
  void vector_value(const std::vector<Point<dim> >&, std::vector<std::array<double, Equations<equationsType, dim>::n_components> >&) const;
};

#endif

#include "completeEllipticIntegrals.h"
#include "initialCondition.h"
#include "equationsMhd.h"

template <EquationsType equationsType, int dim>
InitialCondition<equationsType, dim>::InitialCondition(Parameters<dim>& parameters) : parameters(parameters)
{
};

template class InitialCondition<EquationsTypeMhd, 3>;

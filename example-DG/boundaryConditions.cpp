#include "common.h"

d DirichletBoundaryCondition::calculate(ui component, Point<DIM> point)
{
  switch (component)
  {
  case 0:
    return RHO_EXT;
  case 1:
    return RHO_EXT * V1_EXT;
  case 2:
    return RHO_EXT * V2_EXT;
  case 3:
    return RHO_EXT * V3_EXT;
  case 4:
    return E_EXT;
  }

  abort();
  return 0.;
}

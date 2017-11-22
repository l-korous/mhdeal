#include "boundaryConditions.h"

using namespace dealii;

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
    return B1_EXT;
  case 5:
    return B2_EXT;
  case 6:
    return B3_EXT;
  case 7:
    return E_EXT;
  case 8:
    return J1_EXT;
  case 9:
    return J2_EXT;
  case 10:
    return J3_EXT;
  }

  abort();
  return 0.;
}

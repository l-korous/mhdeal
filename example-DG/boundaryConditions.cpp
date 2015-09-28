#include "common.h"

d DirichletBoundaryCondition::calculate(ui component, Point<DIM> point)
{
  switch (component)
  {
  case 0:
    if (point(0) > 1e-8)
      return RHO_IN_TOP;
    else
      return RHO_IN_LEFT;
  case 1:
    if(point(0) > 1e-8)
      return RHO_IN_TOP * V1_IN_TOP;
    else
      return RHO_IN_LEFT * V1_IN_LEFT;
  case 2:
    if (point(0) > 1e-8)
      return RHO_IN_TOP * V2_IN_TOP;
    else
      return RHO_IN_LEFT * V2_IN_LEFT;
  case 3:
    if (point(0) > 1e-8)
      return RHO_IN_TOP * V3_IN_TOP;
    else
      return RHO_IN_LEFT * V3_IN_LEFT;
  case 4:
    if (point(0) > 1e-8)
      return E_IN_TOP;
    else
      return E_IN_LEFT;
  }

  abort();
  return 0.;
}
